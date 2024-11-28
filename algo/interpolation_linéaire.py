import itertools
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import mu_0
from scipy.integrate import quad
import math

class MagneticFieldSimulation:
    def __init__(self, R=0.1, I=1.0, N=10, d=0.01, resolution=1, k=10, S=1):
        self.R = R
        self.I = I
        self.N = N
        self.d = d
        self.resolution = resolution  # Nombre de subdivisions par axe
        self.k = k  # Constante de propagation
        self.S = S  # Facteur d'échelle
        self.mu0 = 4 * math.pi * 1e-7  # Perméabilité du vide en T·m/A
        self.points_haute_resolution = []
        self.resultats = []

    def magnetic_field_helix(self, x, y, z, R=0.1, I=1.0):
        def dB_element(phi):
            x0 = R * np.cos(phi)
            y0 = R * np.sin(phi)
            z0 = 0
            r_vec = np.array([x - x0, y - y0, z - z0])
            r_mag = np.linalg.norm(r_vec)
            if r_mag == 0:
                return np.array([0, 0, 0])
            dL = np.array([-R * np.sin(phi), R * np.cos(phi), 0])
            dB = (self.mu0 * I / (4 * np.pi)) * np.cross(dL, r_vec) / (r_mag ** 3)
            return dB

        Bx, _ = quad(lambda phi: dB_element(phi)[0], 0, 2 * np.pi)
        By, _ = quad(lambda phi: dB_element(phi)[1], 0, 2 * np.pi)
        Bz, _ = quad(lambda phi: dB_element(phi)[2], 0, 2 * np.pi)
        B_total = np.array([Bx, By, Bz])
        return B_total

    def generate_simulated_points(self):
        x_min, x_max = -0.1, 0.1
        y_min, y_max = -0.1, 0.1
        z_min, z_max = -0.1, 0.1

        self.resultats = [{'x': x, 'y': y, 'z': z} for x in np.linspace(x_min, x_max, 5)
                          for y in np.linspace(y_min, y_max, 5)
                          for z in np.linspace(z_min, z_max, 5)]

        for point in self.resultats:
            Bx, By, Bz = self.magnetic_field_helix(point['x'], point['y'], point['z'])
            Hx, Hy, Hz = Bx / self.mu0, By / self.mu0, Bz / self.mu0
            point.update({'Hx': Hx, 'Hy': Hy, 'Hz': Hz})
        return self.resultats

    def interpoler_trilineaire(self, sommets, u, v, w):
        """Interpolation trilineaire entre 8 sommets."""
        # Extraire Hx, Hy, Hz des sommets
        Hx = [s['Hx'] for s in sommets]
        Hy = [s['Hy'] for s in sommets]
        Hz = [s['Hz'] for s in sommets]

        # Formule d'interpolation trilineaire pour chaque composant
        def interp(values):
            return (
                values[0] * (1 - u) * (1 - v) * (1 - w) +
                values[1] * u * (1 - v) * (1 - w) +
                values[2] * (1 - u) * v * (1 - w) +
                values[3] * u * v * (1 - w) +
                values[4] * (1 - u) * (1 - v) * w +
                values[5] * u * (1 - v) * w +
                values[6] * (1 - u) * v * w +
                values[7] * u * v * w
            )

        Hx_interp = interp(Hx)
        Hy_interp = interp(Hy)
        Hz_interp = interp(Hz)

        # Calculer les coordonnées interpolées
        x_coords = [s['x'] for s in sommets]
        y_coords = [s['y'] for s in sommets]
        z_coords = [s['z'] for s in sommets]

        x_interp = (
            x_coords[0] * (1 - u) * (1 - v) * (1 - w) +
            x_coords[1] * u * (1 - v) * (1 - w) +
            x_coords[2] * (1 - u) * v * (1 - w) +
            x_coords[3] * u * v * (1 - w) +
            x_coords[4] * (1 - u) * (1 - v) * w +
            x_coords[5] * u * (1 - v) * w +
            x_coords[6] * (1 - u) * v * w +
            x_coords[7] * u * v * w
        )

        y_interp = (
            y_coords[0] * (1 - u) * (1 - v) * (1 - w) +
            y_coords[1] * u * (1 - v) * (1 - w) +
            y_coords[2] * (1 - u) * v * (1 - w) +
            y_coords[3] * u * v * (1 - w) +
            y_coords[4] * (1 - u) * (1 - v) * w +
            y_coords[5] * u * (1 - v) * w +
            y_coords[6] * (1 - u) * v * w +
            y_coords[7] * u * v * w
        )

        z_interp = (
            z_coords[0] * (1 - u) * (1 - v) * (1 - w) +
            z_coords[1] * u * (1 - v) * (1 - w) +
            z_coords[2] * (1 - u) * v * (1 - w) +
            z_coords[3] * u * v * (1 - w) +
            z_coords[4] * (1 - u) * (1 - v) * w +
            z_coords[5] * u * (1 - v) * w +
            z_coords[6] * (1 - u) * v * w +
            z_coords[7] * u * v * w
        )

        return {
            'x': x_interp,
            'y': y_interp,
            'z': z_interp,
            'Hx': Hx_interp,
            'Hy': Hy_interp,
            'Hz': Hz_interp
        }

    def augmenter_resolution(self, points):
        """Augmente la résolution de la grille avec interpolation trilineaire."""
        interpolated_points = []
        grid_points = np.array([[point['x'], point['y'], point['z']] for point in points])

        # Extraire les coordonnées uniques de basse résolution
        x_unique = np.unique(grid_points[:, 0])
        y_unique = np.unique(grid_points[:, 1])
        z_unique = np.unique(grid_points[:, 2])

        x_unique_sorted = np.sort(x_unique)
        y_unique_sorted = np.sort(y_unique)
        z_unique_sorted = np.sort(z_unique)

        # Indexation des points pour accès rapide
        points_dict = {(round(p['x'], 5), round(p['y'], 5), round(p['z'], 5)): p for p in points}

        # Fonction de recherche tolérante
        def find_point(x, y, z):
            rounded_key = (round(x, 5), round(y, 5), round(z, 5))
            return points_dict.get(rounded_key)

        # Parcours de chaque "cube" défini par les sommets voisins de basse résolution
        for i in range(len(x_unique_sorted) - 1):
            for j in range(len(y_unique_sorted) - 1):
                for k in range(len(z_unique_sorted) - 1):
                    # Récupérer les 8 sommets du cube de basse résolution
                    cube = []
                    for dx, dy, dz in itertools.product([0, 1], repeat=3):
                        x = x_unique_sorted[i + dx]
                        y = y_unique_sorted[j + dy]
                        z = z_unique_sorted[k + dz]
                        point = find_point(x, y, z)
                        if point is None:
                            raise ValueError(f"Point not found for coordinates x={x}, y={y}, z={z}")
                        cube.append(point)

                    # Interpolation pour `resolution + 1` subdivisions
                    steps = self.resolution + 1
                    for u, v, w in itertools.product(np.linspace(0, 1, steps), repeat=3):
                        # Éviter d'ajouter des points qui coïncident avec les sommets du cube
                        if (u == 0 and v == 0 and w == 0) or (u == 1 and v == 1 and w == 1):
                            continue
                        interpolated_point = self.interpoler_trilineaire(cube, u, v, w)
                        interpolated_points.append(interpolated_point)

        self.points_haute_resolution = interpolated_points

    def afficher_vecteurs_3D(self, points=None):
        points_base = self.resultats
        points_interpolés = points or self.points_haute_resolution

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Points bleus : Points initiaux (sommets)
        for point in points_base:
            x, y, z = point['x'], point['y'], point['z']
            Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
            if x == 0:
                ax.quiver(x, y, z, Hx, Hy, Hz, color='b', length=0.01, normalize=True)

        # Points rouges : Points interpolés
        for point in points_interpolés:
            x, y, z = point['x'], point['y'], point['z']
            Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
            if x == 0:
                ax.quiver(x, y, z, Hx, Hy, Hz, color='r', length=0.005, normalize=True)  # Réduire la longueur

        # Légende et affichage
        plt.show()

    def execute_pipeline(self):
        simulated_points = self.generate_simulated_points()
        self.augmenter_resolution(simulated_points)
        self.afficher_vecteurs_3D()

# Exécution avec une résolution plus élevée pour obtenir plus de points par face
simulation = MagneticFieldSimulation(resolution=4)  # Par exemple, résolution 3 pour 4 subdivisions
simulation.execute_pipeline()