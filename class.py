import os
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class ChampMagnetique:
    def __init__(self, dossier, resolution=5):
        self.dossier = dossier
        self.resolution = resolution
        self.c = 3e8  # Vitesse de la lumière en m/s
        self.mu0 = 4 * math.pi * 1e-7  # Perméabilité magnétique du vide en T·m/A
        self.F = 13.56e6  # Fréquence en Hz
        self.omega = 2 * math.pi * self.F  # Pulsation angulaire en rad/s
        self.S = 1e-4  # Surface de l'antenne en m² (1 cm²)
        self.k = self.omega / self.c  # Nombre d'onde en rad/m
        self.resultats = []
        self.points_haute_resolution = []

    def extraire_amplitude_et_position(self, file_content, filename):
        matches = re.findall(r"C(\d+)_pkpk:(\d+\.\d+)", file_content)
        amplitudes = {f"C{axis}": float(amp) for axis, amp in matches}
        pos_match = re.search(r"logMeasureOscillo_\((-?\d+\.?\d*) (-?\d+\.?\d*) (-?\d+\.?\d*)\)", filename)
        x, y, z = map(float, pos_match.groups()) if pos_match else (0.0, 0.0, 0.0)
        return amplitudes.get('C1', 0), 0, 0, x, y, z

    def calcul_champ_magnetique(self, amp_x, amp_y, amp_z):
        Bx, By, Bz = amp_x / (self.S * self.omega), amp_y / (self.S * self.omega), amp_z / (self.S * self.omega)
        Hx, Hy, Hz = Bx / self.mu0, By / self.mu0, Bz / self.mu0
        H_magnitude = math.sqrt(Hx**2 + Hy**2 + Hz**2)
        Htheta = math.degrees(math.acos(Hz / H_magnitude)) if H_magnitude != 0 else 0.0
        Hphi = math.degrees(math.atan2(Hy, Hx))
        return (Hx, Hy, Hz), H_magnitude, Htheta, Hphi

    def traiter_fichiers(self):
        for filename in os.listdir(self.dossier):
            if filename.startswith("logMeasureOscillo"):
                filepath = os.path.join(self.dossier, filename)
                with open(filepath, 'r') as file:
                    file_content = file.read()
                amp_x, amp_y, amp_z, x, y, z = self.extraire_amplitude_et_position(file_content, filename)
                (Hx, Hy, Hz), H_magnitude, Htheta, Hphi = self.calcul_champ_magnetique(amp_x, amp_y, amp_z)
                self.resultats.append({
                    'fichier': filename, 'x': x, 'y': y, 'z': z, 
                    'Hx': Hx, 'Hy': Hy, 'Hz': Hz, '|H|': H_magnitude, 'Htheta': Htheta, 'Hphi': Hphi
                })
        return self.resultats

    def calcul_r_teta(self, x, y, z):
        r = math.sqrt(x ** 2 + y ** 2 + z ** 2)
        teta = math.degrees(math.acos(z / r)) if r != 0 else 0.0
        phi = math.degrees(math.atan2(y, x))
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
            r, tetha = self.calcul_r_teta(point['x'], point['y'], point['z'])
            H_r, H_theta = point['|H|'], point['Htheta']
            print(f"r: {r}, tetha: {tetha}")
            print(f"H_r: {H_r}, H_theta: {H_theta}")
            point['r'], point['tetha'] = r, tetha
            I_r, I_theta = self.calculer_I(H_r, r, tetha), self.calculer_I_Htetha(H_theta, r, tetha)
            print(f"I_r: {I_r}, I_theta: {I_theta}")
            I_moyenne_point = (I_r + I_theta) / 2 if I_r and I_theta else (I_r or I_theta)
            if I_moyenne_point is not None:
                I_valeurs.append(I_moyenne_point)
        return sum(I_valeurs) / len(I_valeurs) if I_valeurs else None

    def calculer_I(self, H_r, r, theta):
        facteur1, facteur2 = 1j / (self.k**2 * r**2), 1 / (self.k**3 * r**3)
        # denominateur = self.S * (self.k**3) * (facteur1 + facteur2) * math.cos(math.radians(theta))
        denominateur = self.S * (self.k**3) * (facteur1 + facteur2) * 1
        return abs(2 * math.pi * H_r / denominateur) if denominateur else None

    def calculer_I_Htetha(self, H_theta, r, theta):
        # sin_theta = np.sin(np.radians(theta))
        sin_theta = 1
        facteur = -1 / (self.k * r) + 1j / (self.k**2 * r**2) + 1 / (self.k**3 * r**3)
        return abs(4 * np.pi * H_theta / (self.S * self.k**3 * facteur * sin_theta)) if sin_theta else None

    def augmenter_resolution(self, points_proches):
        for i in range(len(points_proches) - 1):
            p1, p2 = points_proches[i], points_proches[i + 1]
            self.points_haute_resolution.append(p1)
            self.points_haute_resolution.extend(self.interpoler_points(p1, p2, self.resolution))
        self.points_haute_resolution.append(points_proches[-1])

    def interpoler_points(self, p1, p2, resolution):
        points_interpolés = []
        for i in range(1, resolution):
            fraction = i / resolution
            x = p1['x'] + (p2['x'] - p1['x']) * fraction
            y = p1['y'] + (p2['y'] - p1['y']) * fraction
            z = p1['z'] + (p2['z'] - p1['z']) * fraction
            if x == 0 and y == 0 and z == 0:
                break
                
            r, teta = self.calcul_r_teta(x, y, z)
            
            # Calculer le champ magnétique pour les points interpolés
            H_magnitude = self.calculer_H_total(x, y, z, self.I_moyenne)
            Hx, Hy, Hz = (H_magnitude * np.cos(np.radians(teta)),
                        H_magnitude * np.sin(np.radians(teta)),
                        H_magnitude)  # Exemple de calcul, à adapter selon le champ
            print(f"|H|: {H_magnitude}, r: {r}, tetha: {teta}")
            points_interpolés.append({
                'x': x, 'y': y, 'z': z, 'Hx': Hx, 'Hy': Hy, 'Hz': Hz, '|H|': H_magnitude, 'r': r, 'tetha': teta
            })
        return points_interpolés


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

    def afficher_vecteurs_3D(self, points=None):
        points = points or self.points_haute_resolution
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for point in points:
            x, y, z = point['x'], point['y'], point['z']
            Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
            ax.quiver(x, y, z, Hx, Hy, Hz, color='b', length=1)
        plt.show()

    def execute_pipeline(self):
        self.traiter_fichiers()
        points_proches = self.selectionner_points_proches()
        
        self.I_moyenne = self.moyenne_I(points_proches)
        self.augmenter_resolution(points_proches)
        self.afficher_vecteurs_3D()

# Définir le chemin du dossier contenant les fichiers de mesures
dossier_mesures = "logMeasure"  # Remplacez par le chemin réel de votre dossier

# Créer une instance de la classe avec la résolution souhaitée pour interpolation
resolution_interpolation = 10  # Définissez la résolution (nombre de points supplémentaires entre chaque paire de mesures)

champ_magnetique = ChampMagnetique(dossier=dossier_mesures, resolution=resolution_interpolation)

# Exécuter la séquence complète
champ_magnetique.execute_pipeline()