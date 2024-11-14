
from matplotlib.collections import PolyCollection
from skimage import measure # pip install scikit-image
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import mu_0
from scipy.integrate import quad
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

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
            self.S = 1e-4  # Surface de l'antenne en m² (1 cm²)
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
        
        return self.resultats

    def calcul_r_teta_phi(self, x, y, z):
        r = math.sqrt(x ** 2 + y ** 2 + z ** 2)
        teta = math.acos(z / r) if r != 0 else 0.0
        phi = math.atan2(y, x)
        return r, teta, phi

    import math

    def selectionner_points_proches(self):
        # Calculer la distance r pour chaque point généré
        for point in self.resultats:
            point['r'] = math.sqrt(point['x'] ** 2 + point['y'] ** 2 + point['z'] ** 2)

        # Trier les points par leur distance r (distance la plus petite en premier)
        points_tries = sorted(self.resultats, key=lambda point: point['r'])

        # Sélectionner les 6 points avec les distances les plus petites
        points_proches = points_tries[:6]

        # Afficher les résultats
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
            H_r, H_theta, H_phi = self.calcul_r_teta_phi(point['Hx'], point['Hy'], point['Hz'])
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
        return abs(2 * math.pi * H_r / denominateur) if denominateur else None

    def calculer_I_Htetha(self, H_theta, r, theta):
        sin_theta = np.sin(theta)
        facteur = -1 / (self.k * r) + 1j / (self.k**2 * r**2) + 1 / (self.k**3 * r**3)
        return abs(4 * np.pi * H_theta / (self.S * self.k**3 * facteur * sin_theta)) if sin_theta else None
    
    # def calculer_H_total(self, x, y, z, I):
    #     Hr, Htheta = self.calculer_Hr(x, y, z, I), self.calculer_Htheta(x, y, z, I)
    #     return math.sqrt(Hr**2 + Htheta**2)
    
    def calculer_Hr(self, x, y, z, I):
        if x == 0 and y == 0 and z == 0:
            return 0
        r, theta, phi = self.calcul_r_teta_phi(x, y, z)
        
        facteur = (1j / (self.k**2 * r**2) + 1 / (self.k**3 * r**3)) * np.exp(-1j * self.k * r)
        
        Hr = (I * self.S * self.k**3 / (2 * np.pi)) * facteur * np.cos(theta)
        return abs(Hr)

    def calculer_Htheta(self, x, y, z, I):
        if x == 0 and y == 0 and z == 0:
            return 0
        r, theta, phi = self.calcul_r_teta_phi(x, y, z)
        sin_theta = np.sin(theta)
        
        facteur = (-1 / (self.k * r) + 1j / (self.k**2 * r**2) + 1 / (self.k**3 * r**3)) * np.exp(-1j * self.k * r)
        
        Htheta = (I * self.S * self.k**3 / (4 * np.pi)) * facteur * sin_theta
        return abs(Htheta)
    
    def convertir_spherique_to_cartesien(self, Hr, Htheta, Hphi, r, theta, phi):
        x = Hr * np.sin(theta) * np.cos(phi) + Htheta * np.cos(theta) * np.cos(phi) - Hphi * np.sin(phi)
        y = Hr * np.sin(theta) * np.sin(phi) + Htheta * np.cos(theta) * np.sin(phi) +  Hphi * np.cos(phi)
        z = Hr * np.cos(theta) - Htheta * np.sin(theta)
        return x, y, z

    
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
            
            # Calculer les composantes en coordonnées cylindriques
            Hr = self.calculer_Hr(x, y, z, I)
            Htheta = self.calculer_Htheta(x, y, z, I)
            Hphi = 0  # D'après l'équation donnée
            
            r, theta, phi = self.calcul_r_teta_phi(x, y, z)
            # Convertir en coordonnées cartésiennes
            Hx, Hy, Hz = self.convertir_spherique_to_cartesien(Hr, Htheta, Hphi, r, theta, phi)
            
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
        
    def afficher_surface_2D(self):
        # Créer une grille 2D
        x = np.linspace(-0.1, 0.1, 10)
        y = np.linspace(-0.1, 0.1, 10)
        X, Y = np.meshgrid(x, y)
        Z = np.zeros_like(X)
        
        # Calculer H pour chaque point de la grille
        for i in range(len(x)):
            for j in range(len(y)):
                r = np.sqrt(X[i,j]**2 + Y[i,j]**2)
                if r > 0:
                    theta = np.arccos(0/r)  # z=0 pour la vue 2D
                    phi = np.arctan2(Y[i,j], X[i,j])
                    
                    Hr = self.calculer_Hr(X[i,j], Y[i,j], 0, self.I)
                    Htheta = self.calculer_Htheta(X[i,j], Y[i,j], 0, self.I)
                    Hx, Hy, _ = self.convertir_spherique_to_cartesien(Hr, Htheta, 0, theta, phi)
                    Z[i,j] = np.sqrt(Hx**2 + Hy**2)  # Amplitude du champ H
        
        # Créer une figure avec deux sous-graphiques
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Surface 3D du champ H
        ax1 = fig.add_subplot(121, projection='3d')
        surf = ax1.plot_surface(X, Y, Z, cmap='viridis')
        ax1.set_title('Amplitude du champ H')
        fig.colorbar(surf, ax=ax1, shrink=0.5, aspect=5)
        
        # Champ de vecteurs
        ax2 = fig.add_subplot(122)
        skip = 3  # Pour réduire la densité des vecteurs
        U = np.zeros_like(X)
        V = np.zeros_like(Y)
        
        for i in range(len(x)):
            for j in range(len(y)):
                r = np.sqrt(X[i,j]**2 + Y[i,j]**2)
                if r > 0:
                    theta = np.arccos(0/r)
                    phi = np.arctan2(Y[i,j], X[i,j])
                    Hr = self.calculer_Hr(X[i,j], Y[i,j], 0, self.I)
                    Htheta = self.calculer_Htheta(X[i,j], Y[i,j], 0, self.I)
                    U[i,j], V[i,j], _ = self.convertir_spherique_to_cartesien(Hr, Htheta, 0, theta, phi)

        ax2.quiver(X[::skip, ::skip], Y[::skip, ::skip], 
                U[::skip, ::skip], V[::skip, ::skip],
                scale=50)
        ax2.set_title('Champ de vecteurs H')
        ax2.set_xlim(-0.1, 0.1)
        ax2.set_ylim(-0.1, 0.1)
        
        plt.tight_layout()
        plt.show()

    def afficher_isosurface(self):
        x = np.linspace(-0.1, 0.1, 50)
        y = np.linspace(-0.1, 0.1, 50)
        X, Y = np.meshgrid(x, y)
        Z = np.zeros_like(X)
        U = np.zeros_like(X)
        V = np.zeros_like(Y)
        
        # Calculer l'amplitude du champ H et ses composantes
        for i in range(len(x)):
            for j in range(len(y)):
                r = np.sqrt(X[i,j]**2 + Y[i,j]**2)
                if r > 0:
                    theta = np.arccos(0/r)
                    phi = np.arctan2(Y[i,j], X[i,j])
                    Hr = self.calculer_Hr(X[i,j], Y[i,j], 0, self.I)
                    Htheta = self.calculer_Htheta(X[i,j], Y[i,j], 0, self.I)
                    Hx, Hy, _ = self.convertir_spherique_to_cartesien(Hr, Htheta, 0, theta, phi)
                    Z[i,j] = np.sqrt(Hx**2 + Hy**2)
                    U[i,j], V[i,j] = Hx, Hy
        
        # Créer une figure avec deux sous-graphiques
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Isosurfaces
        levels = np.linspace(np.min(Z), np.max(Z), 20)
        cs = ax1.contour(X, Y, Z, levels=levels)
        ax1.clabel(cs, inline=True, fontsize=8)
        ax1.set_title('Isosurfaces du champ H')
        ax1.set_xlim(-0.1, 0.1)
        ax1.set_ylim(-0.1, 0.1)
        
        # Champ de vecteurs
        skip = 3
        ax2.quiver(X[::skip, ::skip], Y[::skip, ::skip], 
                U[::skip, ::skip], V[::skip, ::skip],
                scale=50)
        ax2.set_title('Champ de vecteurs H')
        ax2.set_xlim(-0.1, 0.1)
        ax2.set_ylim(-0.1, 0.1)
        
        plt.tight_layout()
        plt.show()

    def execute_pipeline(self):
        simulated_points = self.generate_simulated_points()
        points_proches = self.selectionner_points_proches()
        I_moyen = self.moyenne_I(points_proches)
        
        print(I_moyen)
        
        self.augmenter_resolution(simulated_points, I_moyen)
        self.afficher_vecteurs_3D()
        self.afficher_isosurface()
        self.afficher_surface_2D()
# Execution
simulation = MagneticFieldSimulation()
simulation.execute_pipeline()
