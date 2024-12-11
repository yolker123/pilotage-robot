import math

import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from scipy.constants import mu_0
from scipy.integrate import quad


class MagneticFieldSimulation_brouillon(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax = fig.add_subplot(111, projection='3d')
        self.R = 0.1
        self.I = 1.0
        self.N = 10
        self.d = 0.01
        self.k = 10
        self.S = 1
        super().__init__(fig)

        self.resultats = []  # Points à afficher
        self.new_points_buffer = []  # Buffer pour les nouveaux points

    def magnetic_field_helix(self, x, y, z):
        """Calcul du champ magnétique pour un point (x, y, z)."""

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

        def dB_component_Bx(phi, z0):
            return dB_component(phi, z0, 0)

        def dB_component_By(phi, z0):
            return dB_component(phi, z0, 1)

        def dB_component_Bz(phi, z0):
            return dB_component(phi, z0, 2)

        B_total = np.array([0.0, 0.0, 0.0])
        for i in range(self.N):
            z0 = i * self.d

            # Diviser l'intégrale en plusieurs sous-intervalles pour éviter les singularités
            Bx1, _ = quad(dB_component_Bx, 0, np.pi / 2, args=(z0,))
            Bx2, _ = quad(dB_component_Bx, np.pi / 2, np.pi, args=(z0,))
            Bx3, _ = quad(dB_component_Bx, np.pi, 3 * np.pi / 2, args=(z0,))
            Bx4, _ = quad(dB_component_Bx, 3 * np.pi / 2, 2 * np.pi, args=(z0,))

            By1, _ = quad(dB_component_By, 0, np.pi / 2, args=(z0,))
            By2, _ = quad(dB_component_By, np.pi / 2, np.pi, args=(z0,))
            By3, _ = quad(dB_component_By, np.pi, 3 * np.pi / 2, args=(z0,))
            By4, _ = quad(dB_component_By, 3 * np.pi / 2, 2 * np.pi, args=(z0,))

            Bz1, _ = quad(dB_component_Bz, 0, np.pi / 2, args=(z0,))
            Bz2, _ = quad(dB_component_Bz, np.pi / 2, np.pi, args=(z0,))
            Bz3, _ = quad(dB_component_Bz, np.pi, 3 * np.pi / 2, args=(z0,))
            Bz4, _ = quad(dB_component_Bz, 3 * np.pi / 2, 2 * np.pi, args=(z0,))

            # Additionner les résultats des sous-intervalles
            B_total += np.array([Bx1 + Bx2 + Bx3 + Bx4,
                                 By1 + By2 + By3 + By4,
                                 Bz1 + Bz2 + Bz3 + Bz4])

        return B_total

    def generate_simulated_points(self):
        """Simule les points pour le champ magnétique."""
        self.resultats = [{'x': x, 'y': y, 'z': z} for x in np.linspace(-0.1, 0.1, 5)
                          for y in np.linspace(-0.1, 0.1, 5)
                          for z in np.linspace(-0.1, 0.1, 5)]

    def afficher_vecteurs_3D(self):
        """Affiche les nouveaux points périodiquement (toutes les 250 ms)."""
        if not self.new_points_buffer:
            return

        # Ajouter les nouveaux points au graphique
        for point in self.new_points_buffer:
            x, y, z = point['x'], point['y'], point['z']
            Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
            self.ax.quiver(x, y, z, Hx, Hy, Hz, color='b', length=0.01, normalize=True)

        self.new_points_buffer.clear()
        self.draw()

    def calcul_r_teta(self, x, y, z):
        r = math.sqrt(x ** 2 + y ** 2 + z ** 2)
        teta = math.degrees(math.acos(z / r)) if r != 0 else 0.0
        return r, teta

    def selectionner_points_proches(self):
        resultats_sans_origine = [r for r in self.resultats if not (r['x'] == 0 and r['y'] == 0 and r['z'] == 0)]
        points_selectionnes = []

        for axis in ['x', 'y', 'z']:
            axis_points = sorted(resultats_sans_origine,
                                 key=lambda r: abs(r[axis]) if all(r[a] == 0 for a in 'xyz' if a != axis) else float(
                                     'inf'))[:2]
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
        facteur1 = 1j / (self.k ** 2 * r ** 2)
        facteur2 = 1 / (self.k ** 3 * r ** 3)
        denominateur = self.S * (self.k ** 3) * (facteur1 + facteur2) * math.cos(math.radians(theta))
        return abs(2 * math.pi * H_r / denominateur) if denominateur else None

    def calculer_I_Htetha(self, H_theta, r, theta):
        sin_theta = np.sin(np.radians(theta))
        facteur = -1 / (self.k * r) + 1j / (self.k ** 2 * r ** 2) + 1 / (self.k ** 3 * r ** 3)
        return abs(4 * np.pi * H_theta / (self.S * self.k ** 3 * facteur * sin_theta)) if sin_theta else None

    def calculer_H_total(self, x, y, z, I):
        Hr, Htheta = self.calculer_Hr(x, y, z, I), self.calculer_Htheta(x, y, z, I)
        return math.sqrt(Hr ** 2 + Htheta ** 2)

    def calculer_Hr(self, x, y, z, I):
        r, theta = self.calcul_r_teta(x, y, z)
        facteur = (1j / (self.k ** 2 * r ** 2)) + (1 / (self.k ** 3 * r ** 3))
        Hr = (I * self.S * self.k ** 3 / (2 * math.pi)) * facteur * math.cos(math.radians(theta)) * np.exp(
            -1j * self.k * r)
        return abs(Hr)

    def calculer_Htheta(self, x, y, z, I):
        r, theta = self.calcul_r_teta(x, y, z)
        sin_theta = np.sin(np.radians(theta))
        facteur = (-1 / (self.k * r) + 1j / (self.k ** 2 * r ** 2) + 1 / (self.k ** 3 * r ** 3))
        Htheta = (I * self.S * self.k ** 3 / (4 * math.pi)) * facteur * sin_theta * np.exp(-1j * self.k * r)
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
            Htheta = self.calculer_Htheta(x, y, z, I)
            # put in cartesian coordinates
            Hx = Hr * np.cos(np.radians(y))
            Hy = Hr * np.sin(np.radians(y))
            Hz = Htheta
            interpolated.append({'x': x, 'y': y, 'z': z, 'Hx': Hx, 'Hy': Hy, 'Hz': Hz})
        return interpolated

