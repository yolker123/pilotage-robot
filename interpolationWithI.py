import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import mu_0
from scipy.integrate import quad
import math

class MagneticFieldSimulation:
    def __init__(self, R=0.1, I=1.0, N=10, d=0.01, resolution=2, k=10, S=1):
        self.R = R
        self.I = I
        self.N = N
        self.d = d
        self.resolution = resolution
        self.k = k  # Propagation constant
        self.S = S  # Scaling factor
        self.points_haute_resolution = []
        self.resultats = []

    def magnetic_field_helix(self, x, y, z):
        def dB_component(phi, z0, component):
            x0 = self.R * np.cos(phi)
            y0 = self.R * np.sin(phi)

            r_vec = np.array([x - x0, y - y0, z - z0])
            r_mag = np.linalg.norm(r_vec)
            if r_mag == 0:
                return 0

            dL = np.array([-self.R * np.sin(phi), self.R * np.cos(phi), 0])
            dB_vec = (mu_0 * self.I / (4 * np.pi)) * np.cross(dL, r_vec) / (r_mag ** 3)
            return dB_vec[component]

        B_total = np.array([0.0, 0.0, 0.0])

        for i in range(self.N):
            z0 = i * self.d
            Bx, _ = quad(lambda phi: dB_component(phi, z0, 0), 0, 2 * np.pi, limit=100)
            By, _ = quad(lambda phi: dB_component(phi, z0, 1), 0, 2 * np.pi, limit=100)
            Bz, _ = quad(lambda phi: dB_component(phi, z0, 2), 0, 2 * np.pi, limit=100)
            B_total += np.array([Bx, By, Bz])

        return B_total

    def generate_simulated_points(self):
        """Simulate grid points for initial magnetic field calculation."""
        self.resultats = [{'x': x, 'y': y, 'z': z} for x in np.linspace(-0.1, 0.1, 5)
                          for y in np.linspace(-0.1, 0.1, 5)
                          for z in np.linspace(-0.1, 0.1, 5)]
        for point in self.resultats:
            Bx, By, Bz = self.magnetic_field_helix(point['x'], point['y'], point['z'])
            Hx, Hy, Hz = Bx / mu_0, By / mu_0, Bz / mu_0
            point.update({'Hx': Hx, 'Hy': Hy, 'Hz': Hz})
        return self.resultats

    def calcul_r_teta(self, x, y, z):
        r = math.sqrt(x ** 2 + y ** 2 + z ** 2)
        teta = math.acos(z / r) if r != 0 else 0.0
        return r, teta

    def selectionner_points_proches(self):
        resultats_sans_origine = [r for r in self.resultats if not (r['x'] == 0 and r['y'] == 0 and r['z'] == 0)]
        points_selectionnes = []

        for axis in ['x', 'y', 'z']:
            axis_points = sorted(resultats_sans_origine,
                                 key=lambda r: abs(r[axis]) if all(r[a] == 0 for a in 'xyz' if a != axis) else float('inf'))[:2]
            if len(axis_points) > 1 and (axis_points[0][axis] * axis_points[1][axis] > 0):
                axis_points = [axis_points[0]]
            points_selectionnes.extend(axis_points)
        return points_selectionnes

    def moyenne_I(self, points_proches):
        I_valeurs = []
        for point in points_proches:
            r, teta = self.calcul_r_teta(point['x'], point['y'], point['z'])
            H_r, H_theta = point['Hx'], point['Hy']  # Simplified for now
            I_r = self.calculer_I(H_r, r, teta)
            I_theta = self.calculer_I_Htetha(H_theta, r, teta)
            I_valeurs.append((I_r + I_theta) / 2 if I_r and I_theta else (I_r or I_theta))
        return sum(I_valeurs) / len(I_valeurs) if I_valeurs else None

    def calculer_I(self, H_r, r, theta):
        facteur1 = 1j / (self.k**2 * r**2)
        facteur2 = 1 / (self.k**3 * r**3)
        denominateur = self.S * (self.k**3) * (facteur1 + facteur2) * math.cos(math.radians(theta))
        return abs(2 * math.pi * H_r / denominateur) if denominateur else None

    def calculer_I_Htetha(self, H_theta, r, theta):
        sin_theta = np.sin(np.radians(theta))
        facteur = -1 / (self.k * r) + 1j / (self.k**2 * r**2) + 1 / (self.k**3 * r**3)
        return abs(4 * np.pi * H_theta / (self.S * self.k**3 * facteur * sin_theta)) if sin_theta else None
    
    def calculer_H_total(self, x, y, z, I):
        Hr, Htheta = self.calculer_Hr(x, y, z, I), self.calculer_Htheta(x, y, z, I)
        return math.sqrt(Hr**2 + Htheta**2)

    def calculer_Hr(self, x, y, z, I):
        r, theta = self.calcul_r_teta(x, y, z)
        facteur = (1j / (self.k**2 * r**2)) + (1 / (self.k**3 * r**3))
        Hr = (I * self.S * self.k**3 / (2 * math.pi)) * facteur * math.cos(math.radians(theta)) * np.exp(-1j * self.k * r)
        return abs(Hr)

    def calculer_Htheta(self, x, y, z, I):
        r, theta = self.calcul_r_teta(x, y, z)
        sin_theta = np.sin(np.radians(theta))
        facteur = (-1 / (self.k * r) + 1j / (self.k**2 * r**2) + 1 / (self.k**3 * r**3))
        Htheta = (I * self.S * self.k**3 / (4 * math.pi)) * facteur * sin_theta * np.exp(-1j * self.k * r)
        return abs(Htheta)
    
    def augmenter_resolution(self, points, I_moyen):
        interpolated_points = []
        for i in range(len(points) - 1):
            p1, p2 = points[i], points[i + 1]
            interpolated_points.extend(self.interpoler_points(p1, p2, I_moyen, self.resolution))
        self.points_haute_resolution = interpolated_points

    def interpoler_points(self, p1, p2, I, resolution):
        interpolated = []
        for i in range(1, resolution):
            fraction = i / resolution
            x = p1['x'] + (p2['x'] - p1['x']) * fraction
            y = p1['y'] + (p2['y'] - p1['y']) * fraction
            z = p1['z'] + (p2['z'] - p1['z']) * fraction
            
            Hr = self.calculer_Hr(x, y, z, I)
            print(f"H_r: {Hr}")
            Htheta = self.calculer_Htheta(x, y, z, I)
            print(f"H_theta: {Htheta}")
            Hphi = 0  # Pas d'influence azimutale dans ce modèle simplifié

            # Conversion des composantes cylindriques en cartésiennes

            Hx = Hr * np.sin(Htheta) * np.cos(Hphi) + Htheta * np.cos(Htheta) * np.cos(Hphi) - Hphi * np.sin(Hphi)
            Hy = Hr * np.sin(Htheta) * np.sin(Hphi) + Htheta * np.cos(Htheta) * np.sin(Hphi) + Hphi * np.cos(Hphi)
            Hz = Hr * np.cos(Htheta) - Htheta * np.sin(Htheta)
                
            print(f"Hx: {Hx}, Hy: {Hy}, Hz: {Hz}")
            interpolated.append({'x': x, 'y': y, 'z': z, 'Hx': Hx, 'Hy': Hy, 'Hz': Hz})
        return interpolated


    def afficher_vecteurs_3D(self, points=None):
        points_base = self.resultats
        points_interpolés = points or self.points_haute_resolution
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        # Affichage des points de base (bleu)
        for point in points_base:
            x, y, z = point['x'], point['y'], point['z']
            Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
            ax.quiver(x, y, z, Hx, Hy, Hz, color='b', length=0.01, normalize=True, label="Base" if 'base' not in locals() else "")
            base = True
        
        # Affichage des points interpolés (rouge)
        for point in points_interpolés:
            x, y, z = point['x'], point['y'], point['z']
            Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
            ax.quiver(x, y, z, Hx, Hy, Hz, color='r', length=0.01, normalize=True, label="Interpolé" if 'interpolé' not in locals() else "")
            interpolé = True
        
        # Gestion des labels de la légende
        handles, labels = ax.get_legend_handles_labels()
        unique_labels = dict(zip(labels, handles))
        ax.legend(unique_labels.values(), unique_labels.keys())

        plt.show()


    def execute_pipeline(self):
        simulated_points = self.generate_simulated_points()
        points_proches = self.selectionner_points_proches()
        I_moyen = self.moyenne_I(points_proches)
        print(I_moyen)
        self.augmenter_resolution(simulated_points, I_moyen)
        self.afficher_vecteurs_3D()

# Execution
simulation = MagneticFieldSimulation()
simulation.execute_pipeline()
