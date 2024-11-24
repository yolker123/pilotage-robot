
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import mu_0
from scipy.integrate import quad
import math


class MagneticFieldSimulation:
    def __init__(self, R=0.05, I=1.0, N=10, a=0.001, resolution=2):
            self.R = R
            self.I = I
            self.N = N
            self.a = a
            self.resolution = resolution
            self.c = 3e8  # Vitesse de la lumière en m/s
            self.mu0 = 4 * math.pi * 1e-7  # Perméabilité magnétique du vide en T·m/A
            self.F = 13.56e6  # Fréquence en Hz
            self.omega = 2 * math.pi * self.F  # Pulsation angulaire en rad/s
            self.S = math.pi * R**2  # Surface de l'antenne en m² (1 cm²)
            self.k = self.omega / self.c  # Nombre d'onde en rad/m
            self.points_haute_resolution = []
            self.resultats = []
            

    def magnetic_field_helix(self, x, y, z, R=0.1, I=1.0):
        def dB_element(phi):
            # Position of the current element
            x0 = R * np.cos(phi)
            y0 = R * np.sin(phi)
            z0 = 0  # Loop lies in the xy-plane, so z0 = 0 for all elements
            # Distance vector from the current element to the observation point (x, y, z)
            r_vec = np.array([x - x0, y - y0, z - z0])
            r_mag = np.linalg.norm(r_vec)
            if r_mag == 0:  # Avoid singularity at the source location
                return np.array([0, 0, 0])
            # Direction of current element (tangent to the loop)
            dL = np.array([-R * np.sin(phi), R * np.cos(phi), 0])
            # Magnetic field contribution from this current element
            dB = (mu_0 * I / (4 * np.pi)) * np.cross(dL, r_vec) / (r_mag ** 3)
            return dB

        Bx, _ = quad(lambda phi: dB_element(phi)[0], 0, 2 * np.pi)
        By, _ = quad(lambda phi: dB_element(phi)[1], 0, 2 * np.pi)
        Bz, _ = quad(lambda phi: dB_element(phi)[2], 0, 2 * np.pi)
        B_total = np.array([Bx, By, Bz])
        return B_total
    
    def generate_simulated_points(self):
        """Simulate grid points for initial magnetic field calculation."""
        x_min, x_max = -0.1, 0.1
        y_min, y_max = -0.1, 0.1
        z_min, z_max = -0.1, 0.1
        
        self.resultats = [{'x': x, 'y': y, 'z': z} for x in np.linspace(x_min, x_max, 5)
                          for y in np.linspace(y_min, y_max, 10)
                          for z in np.linspace(z_min, z_max, 10)]
        
        for point in self.resultats:
            Bx, By, Bz = self.magnetic_field_helix(point['x'], point['y'], point['z'])
            Hx, Hy, Hz = Bx / self.mu0, By / self.mu0, Bz / self.mu0
            point.update({'Hx': Hx, 'Hy': Hy, 'Hz': Hz})
            
            r, theta, phi = self.calcul_r_teta_phi(point['x'], point['y'], point['z'])
            Hr, Htheta, Hphi = self.convertir_cartesien_to_spherique(Hx, Hy, Hz, theta, phi)
            # print(f"Point: {point['x']}, {point['y']}, {point['z']}, Hr: {Hr}, Htheta: {Htheta}, Hphi: {Hphi}")
        return self.resultats

    def calcul_r_teta_phi(self, x, y, z):
        r = math.sqrt(x ** 2 + y ** 2 + z ** 2)
        teta = math.acos(z / r) if r != 0 else 0.0
        phi = math.atan2(y, x)
        return r, teta, phi
    
    def selectionner_points_proches(self):
        for point in self.resultats:
            point['r'] = math.sqrt(point['x'] ** 2 + point['y'] ** 2 + point['z'] ** 2)

        points_tries = sorted(self.resultats, key=lambda point: point['r'])

        points_proches = points_tries[:6]

        print("Les 6 points proches les plus proches :")
        for idx, point in enumerate(points_proches, start=1):
            print(f"Point {idx}: {point}, Distance r: {point['r']}")

        return points_proches

    def moyenne_I(self, points_proches):
        I_valeurs = []
        print("Points proches" , points_proches)
        for point in points_proches:
            print(point["x"], point["y"], point["z"])
            r, teta, phi = self.calcul_r_teta_phi(point['x'], point['y'], point['z'])
            H_r, H_theta, H_phi = self.convertir_cartesien_to_spherique(point['Hx'], point['Hy'], point['Hz'], teta, phi)
            print(teta)
            I_r = self.calculer_I(H_r, r, teta)
            I_theta = self.calculer_I_Htetha(H_theta, r, teta)
            print(f"I_r: {I_r}, I_theta: {I_theta}")
            I_valeurs.append((I_r + I_theta) / 2 if I_r and I_theta else (I_r or I_theta))
        return sum(I_valeurs) / len(I_valeurs) if I_valeurs else None

    def calculer_I(self, H_r, r, theta):
        facteur1 = 1j / (self.k**2 * r**2)
        facteur2 = 1 / (self.k**3 * r**3)
        denominateur = self.S * (self.k**3) * (facteur1 + facteur2) * math.cos(theta)
        res =  (2 * math.pi * H_r / denominateur) if denominateur else None
        return math.sqrt(res.real**2 + abs(res.imag**2)) if res else None

    def calculer_I_Htetha(self, H_theta, r, theta):
        sin_theta = np.sin(theta)
        facteur = -1 / (self.k * r) + 1j / (self.k**2 * r**2) + 1 / (self.k**3 * r**3)
        res = (4 * np.pi * H_theta / (self.S * self.k**3 * facteur * sin_theta)) if sin_theta else None
        return math.sqrt(res.real**2 + abs(res.imag**2))  if res else None
    
    def calculer_Hr(self, x, y, z, I):
        if x == 0 and y == 0 and z == 0:
            return 0
        r, theta, phi = self.calcul_r_teta_phi(x, y, z)
        
        facteur = (1j / (self.k**2 * r**2) + 1 / (self.k**3 * r**3)) 
        
        Hr = (I * self.S * self.k**3 / (2 * np.pi)) * facteur * np.cos(theta)
        modHr = math.sqrt(Hr.real**2 + Hr.imag**2)
        return modHr

    def calculer_Htheta(self, x, y, z, I):
        if x == 0 and y == 0 and z == 0:
            return 0
        r, theta, phi = self.calcul_r_teta_phi(x, y, z)
        sin_theta = np.sin(theta)
        
        facteur = (-1 / (self.k * r) + 1j / (self.k**2 * r**2) + 1 / (self.k**3 * r**3))
        
        Htheta = (I * self.S * self.k**3 / (4 * np.pi)) * facteur * sin_theta
        modHtheta = math.sqrt(Htheta.real**2 + Htheta.imag**2)
        return modHtheta
    
    def convertir_spherique_to_cartesien(self, Hr, Htheta, Hphi, theta, phi):
        Hx = Hr * np.sin(theta) * np.cos(phi) + Htheta * np.cos(theta) * np.cos(phi) - Hphi * np.sin(phi)
        Hy = Hr * np.sin(theta) * np.sin(phi) + Htheta * np.cos(theta) * np.sin(phi) +  Hphi * np.cos(phi)
        Hz = Hr * np.cos(theta) - Htheta * np.sin(theta)
        return Hx, Hy, Hz
    
    def convertir_cartesien_to_spherique(self, Hx, Hy, Hz, theta, phi):
        Hr = Hx * np.sin(theta) * np.cos(phi) + Hy * np.sin(theta) * np.sin(phi) + Hz * np.cos(theta)
        Htheta = Hx * np.cos(theta) * np.cos(phi) + Hy * np.cos(theta) * np.sin(phi) - Hz * np.sin(theta)
        Hphi = -Hx * np.sin(phi) + Hy * np.cos(phi)
        return Hr, Htheta, Hphi
    

    
    def recalcul_hx_hy_hz_point_de_base(self, points, I_moyen):
        new = []
        
        for point in points:
            x, y, z = point['x'], point['y'], point['z']
        
            Hr = self.calculer_Hr(x, y, z, I_moyen)
            Htheta = self.calculer_Htheta(x, y, z, I_moyen)
            Hphi = 0  
            
            r, theta, phi = self.calcul_r_teta_phi(x, y, z)
            Hx, Hy, Hz = self.convertir_spherique_to_cartesien(Hr, Htheta, Hphi, theta, phi)
            
            new.append({'x': x, 'y': y, 'z': z, 'Hx': Hx, 'Hy': Hy, 'Hz': Hz})
            
        self.points_haute_resolution = new

    def afficher_vecteurs_3D(self, points=None):
        points_base = self.resultats
        points_interpolés = points or self.points_haute_resolution
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        # Affichage des points de base (bleu)
        for point in points_base:
            x, y, z = point['x'], point['y'], point['z']
            Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
            if x == 0:
                # print(x, y, z)
                # print(Hx, Hy, Hz)
                ax.quiver(x, y, z, Hx, Hy, Hz, color='b', length=0.01, normalize=True, label="Base" if 'base' not in locals() else "")
            base = True
        
        # Affichage des points interpolés (rouge)
        for point in points_interpolés:
            x, y, z = point['x'], point['y'], point['z']
            Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
            if x == 0:
                # print(x, y, z)
                # print(Hx, Hy, Hz)
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
        print(f"I_moyen:  {I_moyen}")
        self.recalcul_hx_hy_hz_point_de_base(simulated_points, I_moyen)
        self.afficher_vecteurs_3D()
        
        for point in self.points_haute_resolution:
            print(f"Point recaculé: {point['Hx']}, {point['Hy']}, {point['Hz']}")
            for p in self.resultats:
                if (point["x"], point["y"], point["z"]) == (p["x"], p["y"], p["z"]):
                    print(f"Point de base: {p['Hx']}, {p['Hy']}, {p['Hz']}")
                    break
# Execution
simulation = MagneticFieldSimulation()
simulation.execute_pipeline()