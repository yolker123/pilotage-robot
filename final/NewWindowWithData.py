import itertools
import sys
from matplotlib.gridspec import GridSpec
import numpy as np
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QTabWidget, QVBoxLayout, QWidget,
    QHBoxLayout, QLabel, QComboBox
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from scipy.constants import mu_0
from scipy.integrate import quad
import math

from scipy.interpolate import griddata


class MagneticFieldSimulation:
    def __init__(self, resolution=3, R=0.01, I=1.0):
        self.R = R
        self.I = I
        self.resolution = resolution
        self.mu0 = 4 * math.pi * 1e-7
        self.resultats = []
        self.points_haute_resolution = []
        self.generate_simulated_points()
        self.augmenter_resolution(self.resultats)

    def magnetic_field_helix(self, x, y, z):
        """Calcule le champ magnétique d'une spire au point (x, y, z)."""
        def dB_element(phi):
            x0 = self.R * np.cos(phi)
            y0 = self.R * np.sin(phi)
            z0 = 0  # Spire dans le plan xy
            r_vec = np.array([x - x0, y - y0, z - z0])
            r_mag = np.linalg.norm(r_vec)
            if r_mag == 0:  # Évite la singularité
                return np.array([0, 0, 0])
            dL = np.array([-self.R * np.sin(phi), self.R * np.cos(phi), 0])  # Direction du courant
            dB = (mu_0 * self.I / (4 * np.pi)) * np.cross(dL, r_vec) / (r_mag ** 3)
            return dB

        Bx, _ = quad(lambda phi: dB_element(phi)[0], 0, 2 * np.pi)
        By, _ = quad(lambda phi: dB_element(phi)[1], 0, 2 * np.pi)
        Bz, _ = quad(lambda phi: dB_element(phi)[2], 0, 2 * np.pi)
        return np.array([Bx, By, Bz])

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
            H_total = np.linalg.norm([Hx, Hy, Hz])
            point.update({'Hx': Hx, 'Hy': Hy, 'Hz': Hz, 'H_total': H_total})
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
                        Hx, Hy, Hz = interpolated_point['Hx'], interpolated_point['Hy'], interpolated_point['Hz']
                        H_total = np.linalg.norm([Hx, Hy, Hz])
                        interpolated_point.update({'H_total': H_total})
                        interpolated_points.append(interpolated_point)
        
        print(f"Nombre de points interpolés : {len(interpolated_points)}")

        self.points_haute_resolution = interpolated_points


class MagneticFieldApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Simulation du Champ Magnétique")
        self.setGeometry(100, 100, 1200, 800)

        # Initialise la simulation
        self.simulation = MagneticFieldSimulation(resolution=3)
        print(f"Nombre de points originaux: {len(self.simulation.resultats)}")
        print(f"Nombre de points interpolés: {len(self.simulation.points_haute_resolution)}")
        
        # Création de la disposition principale
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)

        # Ajouter les onglets
        self.add_tab_3d_vectors()
        self.add_tab_2d_plane()
        self.add_tab_gaussian_and_radial()


    def add_tab_3d_vectors(self):
        """Onglet 1 : Affichage 3D des vecteurs."""
        tab = QWidget()
        layout = QVBoxLayout()
        tab.setLayout(layout)

        figure = plt.figure()
        canvas = FigureCanvas(figure)
        ax = figure.add_subplot(111, projection='3d')

        # Données pour les vecteurs
        points = self.simulation.resultats
        points_interpolés = self.simulation.points_haute_resolution
        for point in points:
            x, y, z = point['x'], point['y'], point['z']
            Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
            if x == 0:
                ax.quiver(x, y, z, Hx, Hy, Hz, color='b', length=0.01, normalize=True)
        print(f"Nombre de points interpolés pour vecteurs 3D: {len(points_interpolés)}")
        for point in points_interpolés:
            x, y, z = point['x'], point['y'], point['z']
            Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
            if x == 0:
                ax.quiver(x, y, z, Hx, Hy, Hz, color='r', length=0.005, normalize=True)  # Réduire la longueur

        ax.set_title("Vecteurs 3D du champ magnétique")
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        layout.addWidget(canvas)
        self.tabs.addTab(tab, "Vecteurs 3D")

    def add_tab_2d_plane(self):
        """Onglet 2 : Affichage 2D dans un plan avec sélecteur de plan."""
        tab = QWidget()
        layout = QVBoxLayout()
        tab.setLayout(layout)

        # Ajout des sélecteurs
        selector_layout = QHBoxLayout()

        # Sélecteur de plan
        plane_label = QLabel("Plan:")
        self.plane_selector_2d = QComboBox()
        self.plane_selector_2d.addItems(['x', 'y', 'z'])
        self.plane_selector_2d.currentTextChanged.connect(self.update_plane_selector_2d_values)

        # Sélecteur de valeur
        value_label = QLabel("Valeur:")
        self.value_selector_2d = QComboBox()
        self.value_selector_2d.currentTextChanged.connect(self.update_tab_2d_plane)

        selector_layout.addWidget(plane_label)
        selector_layout.addWidget(self.plane_selector_2d)
        selector_layout.addWidget(value_label)
        selector_layout.addWidget(self.value_selector_2d)
        layout.addLayout(selector_layout)

        # Figure et axes
        self.figure_2d = plt.figure(figsize=(12, 6))
        self.canvas_2d = FigureCanvas(self.figure_2d)
        self.axes_2d = self.figure_2d.subplots(1, 2)
        layout.addWidget(self.canvas_2d)

        # Initialiser les valeurs du sélecteur de valeurs
        self.update_plane_selector_2d_values()

        # Tracer initial
        self.update_tab_2d_plane()

        self.tabs.addTab(tab, "Champ dans le plan 2D")

    def update_plane_selector_2d_values(self):
        """Met à jour les valeurs disponibles dans le sélecteur de valeurs pour l'onglet 2."""
        plane = self.plane_selector_2d.currentText()
        # Réinitialiser le sélecteur de valeurs
        self.value_selector_2d.blockSignals(True)  # Empêche le déclenchement de signaux lors de la mise à jour
        self.value_selector_2d.clear()
        if plane == 'x':
            unique_values = sorted(np.unique([p['x'] for p in self.simulation.points_haute_resolution]))
        elif plane == 'y':
            unique_values = sorted(np.unique([p['y'] for p in self.simulation.points_haute_resolution]))
        elif plane == 'z':
            unique_values = sorted(np.unique([p['z'] for p in self.simulation.points_haute_resolution]))
        else:
            unique_values = []
        self.value_selector_2d.addItems([f"{v:.5f}" for v in unique_values])
        self.value_selector_2d.blockSignals(False)  # Réactive les signaux

        # Optionnel : définir une valeur par défaut
        if unique_values:
            self.value_selector_2d.setCurrentIndex(len(unique_values) // 2)  # Sélectionne le milieu par défaut

    def update_tab_2d_plane(self):
        """Met à jour les graphiques de l'onglet Champ dans le plan 2D selon le plan et la valeur sélectionnés."""

        plane = self.plane_selector_2d.currentText()
        if self.value_selector_2d.count() == 0:
            return  # Aucun point à tracer
        try:
            value = float(self.value_selector_2d.currentText())
        except ValueError:
            return  # Valeur non valide

        # Filtrer les points selon le plan et la valeur
        epsilon = 1e-5  # Tolérance pour la comparaison
        if plane == 'x':
            filtered_points = [p for p in self.simulation.points_haute_resolution if np.isclose(p['x'], value, atol=epsilon)]
        elif plane == 'y':
            filtered_points = [p for p in self.simulation.points_haute_resolution if np.isclose(p['y'], value, atol=epsilon)]
        elif plane == 'z':
            filtered_points = [p for p in self.simulation.points_haute_resolution if np.isclose(p['z'], value, atol=epsilon)]
        else:
            filtered_points = []

        if not filtered_points:
            print(f"Aucun point trouvé pour le plan {plane}={value}")
            return

        # Extraire les axes non-constantes
        axes_vars = {'x': ['y', 'z'], 'y': ['x', 'z'], 'z': ['x', 'y']}
        axis1, axis2 = axes_vars[plane]

        coord1 = np.array([p[axis1] for p in filtered_points])
        coord2 = np.array([p[axis2] for p in filtered_points])

        H_total = np.sqrt(
            np.array([p['Hx'] for p in filtered_points]) ** 2 +
            np.array([p['Hy'] for p in filtered_points]) ** 2 +
            np.array([p['Hz'] for p in filtered_points]) ** 2
        )
        H_component1 = np.array([p[f'H{axis1}'] for p in filtered_points])
        H_component2 = np.array([p[f'H{axis2}'] for p in filtered_points])

        # Calculer les vecteurs normalisés pour la direction
        with np.errstate(divide='ignore', invalid='ignore'):
            H_component1_normalized = np.where(H_total != 0, H_component1 / H_total, 0)
            H_component2_normalized = np.where(H_total != 0, H_component2 / H_total, 0)

        # Créer une grille régulière
        unique_coord1 = np.unique(coord1)
        unique_coord2 = np.unique(coord2)
        coord1_grid, coord2_grid = np.meshgrid(unique_coord1, unique_coord2)

        # Interpoler les données sur la grille
        H_total_grid = griddata((coord1, coord2), H_total, (coord1_grid, coord2_grid), method='linear', fill_value=0)
        H_component1_norm_grid = griddata((coord1, coord2), H_component1_normalized, (coord1_grid, coord2_grid), method='linear', fill_value=0)
        H_component2_norm_grid = griddata((coord1, coord2), H_component2_normalized, (coord1_grid, coord2_grid), method='linear', fill_value=0)

        # Créer une nouvelle figure avec GridSpec
        self.figure_2d.clear()  # Effacer la figure existante
        gs = GridSpec(1, 3, width_ratios=[4, 0.2, 4], figure=self.figure_2d)  # 3 colonnes : 4:0.2:4 proportions

        # Sous-plot pour le champ scalaire |H|
        ax1 = self.figure_2d.add_subplot(gs[0, 0])  # Colonne de gauche
        contour = ax1.contourf(coord1_grid, coord2_grid, H_total_grid, levels=20, cmap='viridis')
        ax1.set_title(f"Norme du champ |H| ({axis1}, {axis2})")
        ax1.set_xlabel(f'{axis1} (m)')
        ax1.set_ylabel(f'{axis2} (m)')

        # Ajouter une colorbar fine à côté de ax1
        cbar_ax = self.figure_2d.add_subplot(gs[0, 1])  # Colonne étroite au centre
        self.figure_2d.colorbar(contour, cax=cbar_ax, label='|H| (A/m)')

        # Sous-plot pour le champ vectoriel
        ax2 = self.figure_2d.add_subplot(gs[0, 2])  # Colonne de droite
        quiver = ax2.quiver(coord1_grid, coord2_grid, H_component1_norm_grid, H_component2_norm_grid, color='red', scale=20)
        ax2.set_title(f"Direction du champ magnétique (H{axis1}, H{axis2})")
        ax2.set_xlabel(f'{axis1} (m)')
        ax2.set_ylabel(f'{axis2} (m)')
        ax2.set_xlim(coord1_grid.min(), coord1_grid.max())
        ax2.set_ylim(coord2_grid.min(), coord2_grid.max())
        ax2.set_aspect('equal')  # Assure que les axes ont la même échelle

        # Rafraîchir le canvas
        self.canvas_2d.draw()


    def add_tab_gaussian_and_radial(self):
        """Onglet 3 : Affichage du champ radial et de l'amplitude simulée avec sélecteur de plan."""
        tab = QWidget()
        layout = QVBoxLayout()
        tab.setLayout(layout)

        # Ajout des sélecteurs
        selector_layout = QHBoxLayout()

        # Sélecteur de plan
        plane_label = QLabel("Plan:")
        self.plane_selector_3d = QComboBox()
        self.plane_selector_3d.addItems(['x', 'y', 'z'])
        self.plane_selector_3d.currentTextChanged.connect(self.update_plane_selector_3d_values)

        # Sélecteur de valeur
        value_label = QLabel("Valeur:")
        self.value_selector_3d = QComboBox()
        self.value_selector_3d.currentTextChanged.connect(self.update_tab_gaussian_and_radial)

        selector_layout.addWidget(plane_label)
        selector_layout.addWidget(self.plane_selector_3d)
        selector_layout.addWidget(value_label)
        selector_layout.addWidget(self.value_selector_3d)
        layout.addLayout(selector_layout)

        # Figure
        self.figure_gaussian = plt.figure(figsize=(12, 6))
        self.canvas_gaussian = FigureCanvas(self.figure_gaussian)
        layout.addWidget(self.canvas_gaussian)

        # Initialiser les valeurs du sélecteur de valeurs
        self.update_plane_selector_3d_values()

        # Tracer initial
        self.update_tab_gaussian_and_radial()

        self.tabs.addTab(tab, "Champ Amplitude et Vectoriel")

    def update_plane_selector_3d_values(self):
        """Met à jour les valeurs disponibles dans le sélecteur de valeurs pour l'onglet 3."""
        plane = self.plane_selector_3d.currentText()
        # Réinitialiser le sélecteur de valeurs
        self.value_selector_3d.blockSignals(True)  # Empêche le déclenchement de signaux lors de la mise à jour
        self.value_selector_3d.clear()
        if plane == 'x':
            unique_values = sorted(np.unique([p['x'] for p in self.simulation.points_haute_resolution]))
        elif plane == 'y':
            unique_values = sorted(np.unique([p['y'] for p in self.simulation.points_haute_resolution]))
        elif plane == 'z':
            unique_values = sorted(np.unique([p['z'] for p in self.simulation.points_haute_resolution]))
        else:
            unique_values = []
        self.value_selector_3d.addItems([f"{v:.5f}" for v in unique_values])
        self.value_selector_3d.blockSignals(False)  # Réactive les signaux

        # Optionnel : définir une valeur par défaut
        if unique_values:
            self.value_selector_3d.setCurrentIndex(len(unique_values) // 2)  # Sélectionne le milieu par défaut

    def update_tab_gaussian_and_radial(self):
        """Met à jour les graphiques de l'onglet Champ Amplitude et Vectoriel selon le plan et la valeur sélectionnés."""
        plane = self.plane_selector_3d.currentText()
        if self.value_selector_3d.count() == 0:
            return  # Aucun point à tracer
        try:
            value = float(self.value_selector_3d.currentText())
        except ValueError:
            return  # Valeur non valide

        # Filtrer les points selon le plan et la valeur
        epsilon = 1e-5  # Tolérance pour la comparaison
        if plane == 'x':
            filtered_points = [p for p in self.simulation.points_haute_resolution if np.isclose(p['x'], value, atol=epsilon)]
        elif plane == 'y':
            filtered_points = [p for p in self.simulation.points_haute_resolution if np.isclose(p['y'], value, atol=epsilon)]
        elif plane == 'z':
            filtered_points = [p for p in self.simulation.points_haute_resolution if np.isclose(p['z'], value, atol=epsilon)]
        else:
            filtered_points = []

        if not filtered_points:
            print(f"Aucun point trouvé pour le plan {plane}={value}")
            return

        # Extraire les axes non-constantes
        axes_vars = {'x': ['y', 'z'], 'y': ['x', 'z'], 'z': ['x', 'y']}
        axis1, axis2 = axes_vars[plane]

        coord1 = np.array([p[axis1] for p in filtered_points])
        coord2 = np.array([p[axis2] for p in filtered_points])

        H_total = np.sqrt(
            np.array([p['Hx'] for p in filtered_points]) ** 2 +
            np.array([p['Hy'] for p in filtered_points]) ** 2 +
            np.array([p['Hz'] for p in filtered_points]) ** 2
        )
        H_component1 = np.array([p[f'H{axis1}'] for p in filtered_points])
        H_component2 = np.array([p[f'H{axis2}'] for p in filtered_points])

        # Calculer les vecteurs normalisés pour la direction
        with np.errstate(divide='ignore', invalid='ignore'):
            H_component1_normalized = np.where(H_total != 0, H_component1 / H_total, 0)
            H_component2_normalized = np.where(H_total != 0, H_component2 / H_total, 0)

        # Créer une grille régulière
        unique_coord1 = np.unique(coord1)
        unique_coord2 = np.unique(coord2)
        coord1_grid, coord2_grid = np.meshgrid(unique_coord1, unique_coord2)

        # Interpoler les données sur la grille
        H_total_grid = griddata((coord1, coord2), H_total, (coord1_grid, coord2_grid), method='linear', fill_value=0)
        H_component1_norm_grid = griddata((coord1, coord2), H_component1_normalized, (coord1_grid, coord2_grid), method='linear', fill_value=0)
        H_component2_norm_grid = griddata((coord1, coord2), H_component2_normalized, (coord1_grid, coord2_grid), method='linear', fill_value=0)

        # Effacer la figure existante
        self.figure_gaussian.clf()

        # Utiliser GridSpec pour organiser les subplots
        gs = self.figure_gaussian.add_gridspec(1, 2, width_ratios=[1, 1])

        # **Graphique de gauche : Amplitude du vecteur H en 3D**
        ax1 = self.figure_gaussian.add_subplot(gs[0, 0], projection='3d')

        # Tracer la surface 3D
        surf = ax1.plot_surface(coord1_grid, coord2_grid, H_total_grid, cmap='viridis', edgecolor='k', alpha=0.8)
        ax1.set_title(f'Amplitude du champ magnétique |H| ({axis1}, {axis2})')
        ax1.set_xlabel(f'{axis1} (m)')
        ax1.set_ylabel(f'{axis2} (m)')
        ax1.set_zlabel('Amplitude |H| (A/m)')
        self.figure_gaussian.colorbar(surf, ax=ax1, shrink=0.5, aspect=10, label='|H| (A/m)')

        # **Graphique de droite : Champ Vectoriel 2D (Hy, Hz)**
        ax2 = self.figure_gaussian.add_subplot(gs[0, 1])

        # Tracer le champ vectoriel avec coloration basée sur la magnitude
        quiver = ax2.quiver(coord1_grid, coord2_grid, H_component1_norm_grid, H_component2_norm_grid, H_total_grid, cmap='inferno', scale=20, scale_units='width', angles='xy')
        ax2.set_title(f"Champ vectoriel (H{axis1}, H{axis2})")
        ax2.set_xlabel(f'{axis1} (m)')
        ax2.set_ylabel(f'{axis2} (m)')
        ax2.set_xlim(coord1_grid.min(), coord1_grid.max())
        ax2.set_ylim(coord2_grid.min(), coord2_grid.max())
        ax2.set_aspect('equal')  # Assure que les axes ont la même échelle
        self.figure_gaussian.colorbar(quiver, ax=ax2, shrink=0.5, aspect=10, label='Magnitude |H| (A/m)')

        # Rafraîchir le canvas
        self.canvas_gaussian.draw()


# Lancement de l'application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_window = MagneticFieldApp()
    main_window.show()
    sys.exit(app.exec_())