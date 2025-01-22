
from scipy.constants import mu_0
from scipy.integrate import quad
import itertools
import numpy as np
import math

class MagneticFieldSimulation:
    def __init__(self, resolution=3, R=0.01):
        self.c = 3e8  # Vitesse de la lumière en m/s
        self.F = 13.56e6  # Fréquence en Hz
        self.omega = 2 * math.pi * self.F  # Pulsation angulaire en rad/s
        self.k = self.omega / self.c  # Nombre d'onde en rad/m
        self.S = math.pi * R ** 2
        self.resolution = resolution
        self.mu_0 = 4 * np.pi * 1e-7  # Perméabilité du vide (T·m/A)
        self.resultats = []
        self.points_haute_resolution = []
        # self.read_file_and_calculate()
        # self.augmenter_resolution(self.resultats)

    def read_file_and_calculate_point(self, line, columns):
        # Définir les indices des colonnes que nous voulons extraire
        try:
            index_Ax = columns.index("CH1_Max_Voltage")
            index_Ay = columns.index("CH2_Max_Voltage")
            index_Az = columns.index("CH3_Max_Voltage")
            index_x = columns.index("x")
            index_y = columns.index("y")
            index_z = columns.index("z")
        except ValueError as e:
            print(f"Erreur : Colonne manquante dans les données : {e}")
            return

        # Extraire les valeurs par index
        values = line.split(',')
        try:
            Ax = float(values[index_Ax])
            Ay = float(values[index_Ay])
            Az = float(values[index_Az])
            x = float(values[index_x])
            y = float(values[index_y])
            z = float(values[index_z])
            print(f"Valeurs extraites : Ax={Ax}, Ay={Ay}, Az={Az}, x={x}, y={y}, z={z}")
        except (IndexError, ValueError) as e:
            # Gérer les erreurs possibles lors de l'extraction et conversion
            print(f"Erreur lors de l'extraction des valeurs : {e}")
            return

        # Calculs
        Bx = Ax / (self.S * self.omega)
        By = Ay / (self.S * self.omega)
        Bz = Az / (self.S * self.omega)

        Hx = Bx / self.mu_0
        Hy = By / self.mu_0
        Hz = Bz / self.mu_0

        # Ajouter le résultat au tableau
        self.resultats.append({
            'x': x, 'y': y, 'z': z,
            'Hx': Hx, 'Hy': Hy, 'Hz': Hz
        })
        H = np.linalg.norm([Hx, Hy, Hz])
        print(f"Résultat ajouté : {H} pour x={x}, y={y}, z={z}")

    def interpoler_trilineaire(self, sommets, u, v, w):
        """Interpolation trilineaire entre 8 sommets."""

        # Function to perform trilinear interpolation for a given set of values
        def interpolate(values):
            c00 = values[0] * (1 - u) + values[1] * u
            c01 = values[2] * (1 - u) + values[3] * u
            c10 = values[4] * (1 - u) + values[5] * u
            c11 = values[6] * (1 - u) + values[7] * u

            c0 = c00 * (1 - v) + c01 * v
            c1 = c10 * (1 - v) + c11 * v

            return c0 * (1 - w) + c1 * w

        # Interpolate Hx, Hy, Hz
        Hx = interpolate([s['Hx'] for s in sommets])
        Hy = interpolate([s['Hy'] for s in sommets])
        Hz = interpolate([s['Hz'] for s in sommets])

        # Interpolate coordinates x, y, z
        x = interpolate([s['x'] for s in sommets])
        y = interpolate([s['y'] for s in sommets])
        z = interpolate([s['z'] for s in sommets])

        return {'x': x, 'y': y, 'z': z, 'Hx': Hx, 'Hy': Hy, 'Hz': Hz}


    def augmenter_resolution(self, points, algorithm, I_moyen=0):
        """Augmente la résolution de la grille avec interpolation trilineaire."""
        print(f"Nombre initial de points : {len(points)}")

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
                            raise ValueError(f"Point not found pour les coordonnées x={x}, y={y}, z={z}")
                        cube.append(point)

                    # Interpolation pour `resolution + 1` subdivisions
                    steps = self.resolution + 1
                    for u, v, w in itertools.product(np.linspace(0, 1, steps), repeat=3):
                        # Éviter d'ajouter des points qui coïncident avec les sommets du cube
                        if (u == 0 and v == 0 and w == 0) or (u == 1 and v == 1 and w == 1):
                            continue
                        interpolated_point = self.interpoler_trilineaire(cube, u, v, w)
                        if algorithm == "linear":
                            Hx, Hy, Hz = interpolated_point['Hx'], interpolated_point['Hy'], interpolated_point['Hz']
                            H_total = np.linalg.norm([Hx, Hy, Hz])
                            interpolated_point.update({'H_total': H_total})
                            interpolated_points.append(interpolated_point)
                        else:
                            x, y, z = interpolated_point['x'], interpolated_point['y'], interpolated_point['z']
                            Hr = self.calculer_Hr(x, y, z, I_moyen)
                            Htheta = self.calculer_Htheta(x, y, z, I_moyen)
                            Hphi = 0  # D'après l'équation donnée

                            r, theta, phi = self.calcul_r_teta_phi(x, y, z)
                            Hx, Hy, Hz = self.convertir_spherique_to_cartesien(Hr, Htheta, Hphi, r, theta, phi)

                            interpolated_points.append({'x': x, 'y': y, 'z': z, 'Hx': Hx, 'Hy': Hy, 'Hz': Hz})
        print(f"Nombre de points interpolés : {len(interpolated_points)}")

        self.points_haute_resolution = interpolated_points

# -----------------NONLINEAIRE
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
        print("Points proches", points_proches)
        for point in points_proches:
            print(point["x"], point["y"], point["z"])
            r, teta, phi = self.calcul_r_teta_phi(point['x'], point['y'], point['z'])
            H_r, H_theta, H_phi = self.convertir_cartesien_to_spherique(point['Hx'], point['Hy'], point['Hz'], teta,
                                                                        phi)
            print(teta)
            I_r = self.calculer_I(H_r, r, teta)
            I_theta = self.calculer_I_Htetha(H_theta, r, teta)
            print(f"I_r: {I_r}, I_theta: {I_theta}")
            I_valeurs.append((I_r + I_theta) / 2 if I_r and I_theta else (I_r or I_theta))
        return sum(I_valeurs) / len(I_valeurs) if I_valeurs else None

    def calculer_I(self, H_r, r, theta):
        facteur1 = 1j / (self.k ** 2 * r ** 2)
        facteur2 = 1 / (self.k ** 3 * r ** 3)

        module_facteurs = math.sqrt(facteur1.real ** 2 + facteur1.imag ** 2) + math.sqrt(
            facteur2.real ** 2 + facteur2.imag ** 2)

        denominateur = self.S * (self.k ** 3) * module_facteurs * math.cos(theta)

        res = (2 * math.pi * H_r / denominateur) if denominateur else None
        return res

    def calculer_I_Htetha(self, H_theta, r, theta):
        sin_theta = np.sin(theta)

        facteur = (-1 / (self.k * r)) + (1j / (self.k ** 2 * r ** 2)) + (1 / (self.k ** 3 * r ** 3))
        module_facteur = math.sqrt(facteur.real ** 2 + facteur.imag ** 2)

        denominateur = self.S * self.k ** 3 * module_facteur * sin_theta

        res = (4 * np.pi * H_theta / denominateur) if sin_theta else None
        return res

    def calculer_Hr(self, x, y, z, I):
        if x == 0 and y == 0 and z == 0:
            return 0
        r, theta, phi = self.calcul_r_teta_phi(x, y, z)

        facteur = (1j / (self.k ** 2 * r ** 2) + 1 / (self.k ** 3 * r ** 3))

        # Calculer le module du facteur
        mod_facteur = math.sqrt(facteur.real ** 2 + facteur.imag ** 2)

        Hr = (self.S * self.k ** 3 / (2 * np.pi)) * mod_facteur * np.cos(theta) * I
        return Hr

    def calculer_Htheta(self, x, y, z, I):
        if x == 0 and y == 0 and z == 0:
            return 0
        r, theta, phi = self.calcul_r_teta_phi(x, y, z)

        facteur = (-1 / (self.k * r) + 1j / (self.k ** 2 * r ** 2) + 1 / (self.k ** 3 * r ** 3))

        # Calculer le module du facteur
        mod_facteur = math.sqrt(facteur.real ** 2 + facteur.imag ** 2)

        Htheta = (self.S * self.k ** 3 / (4 * np.pi)) * mod_facteur * np.sin(theta) * I
        return Htheta

    def convertir_spherique_to_cartesien(self, Hr, Htheta, Hphi, r, theta, phi):
        x = Hr * np.sin(theta) * np.cos(phi) + Htheta * np.cos(theta) * np.cos(phi) - Hphi * np.sin(phi)
        y = Hr * np.sin(theta) * np.sin(phi) + Htheta * np.cos(theta) * np.sin(phi) + Hphi * np.cos(phi)
        z = Hr * np.cos(theta) - Htheta * np.sin(theta)
        return x, y, z

    def convertir_cartesien_to_spherique(self, Hx, Hy, Hz, theta, phi):
        Hr = Hx * np.sin(theta) * np.cos(phi) + Hy * np.sin(theta) * np.sin(phi) + Hz * np.cos(theta)
        Htheta = Hx * np.cos(theta) * np.cos(phi) + Hy * np.cos(theta) * np.sin(phi) - Hz * np.sin(theta)
        Hphi = -Hx * np.sin(phi) + Hy * np.cos(phi)
        return Hr, Htheta, Hphi

