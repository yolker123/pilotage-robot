
import sys
import numpy as np
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QTabWidget, QVBoxLayout, QWidget,
    QHBoxLayout, QLabel, QComboBox, QPushButton
)
from PyQt5.QtWidgets import QDialog, QLineEdit, QPushButton, QVBoxLayout, QFormLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtCore import Qt

from scipy.interpolate import griddata

from MagneticFieldSimulation import MagneticFieldSimulation

class MagneticFieldApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Simulation du Champ Magnétique")
        self.setGeometry(100, 100, 1200, 800)

        # Initialise la simulation
        self.simulation = MagneticFieldSimulation(resolution=1)
        self.initUI()

    def initUI(self):
        # Création de la disposition principale
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)

        # Ajouter les onglets
        self.add_tab_3d_vectors()
        self.add_tab_2d_plane()
        self.add_tab_gaussian_and_radial()
        
    def set_resolution(self, resolution_value, dialog):
        try:
            resolution_value = int(resolution_value)
            self.simulation.resolution = resolution_value  # Modifier la résolution de la simulation
            self.simulation.augmenter_resolution(self.simulation.resultats)  # Augmenter la résolution
            self.update_all_graphs()  # Appeler une méthode pour rafraîchir les graphiques
            dialog.accept()  # Fermer la fenêtre popup
            print(f"Résolution modifiée à : {resolution_value}")
        except ValueError:
            print("Veuillez entrer un nombre entier valide.")
            
    def update_all_graphs(self):
        self.update_tab_2d_plane()  # Mettre à jour les graphiques 2D
        self.update_tab_gaussian_and_radial()  # Mettre à jour les graphiques 3D
        self.plot_3d_vectors()  # Mettre à jour les vecteurs 3D
        
    def add_tab_3d_vectors(self):
        """Onglet 1 : Affichage 3D des vecteurs."""
        self.tab_3d = QWidget()  # Créer un attribut pour l'onglet afin de pouvoir le mettre à jour
        layout = QVBoxLayout()
        self.tab_3d.setLayout(layout)

        # Canvas pour le tracé 3D
        self.figure_3d = plt.figure()
        self.canvas_3d = FigureCanvas(self.figure_3d)
        self.ax_3d = self.figure_3d.add_subplot(111, projection='3d')

        # Tracer les vecteurs
        self.plot_3d_vectors()

        # Ajouter le bouton pour modifier la résolution
        resolution_button = QPushButton("Modifier Résolution")
        resolution_button.clicked.connect(self.open_resolution_dialog)
        layout.addWidget(resolution_button, alignment=Qt.AlignRight)

        layout.addWidget(self.canvas_3d)
        self.tabs.addTab(self.tab_3d, "Vecteurs 3D")

    def plot_3d_vectors(self):
        """Dessine les vecteurs 3D dans l'onglet correspondant."""
        self.ax_3d.clear()  # Effacer les anciens vecteurs

        # Données pour les vecteurs
        points = self.simulation.resultats
        points_interpolés = self.simulation.points_haute_resolution

        # Tracer les vecteurs des points originaux
        for point in points:
            x, y, z = point['x'], point['y'], point['z']
            Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
            self.ax_3d.quiver(x, y, z, Hx, Hy, Hz, color='b', length=0.01, normalize=True)
        # Tracer les vecteurs interpolés
        for point in points_interpolés:
            x, y, z = point['x'], point['y'], point['z']
            Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
            self.ax_3d.quiver(x, y, z, Hx, Hy, Hz, color='r', length=0.005, normalize=True)

        # Configurer les axes
        self.ax_3d.set_title("Vecteurs 3D du champ magnétique")
        self.ax_3d.set_xlabel('x')
        self.ax_3d.set_ylabel('y')
        self.ax_3d.set_zlabel('z')

        # Rafraîchir le canvas
        self.canvas_3d.draw()

    def add_tab_2d_plane(self):
        """Onglet 2 : Affichage 2D dans un plan avec sélecteur de plan."""
        self.add_tab_with_plane_selector(
            title="Champ dans le plan 2D",
            update_selector_function=self.update_plane_selector_2d_values,
            update_plot_function=self.update_tab_2d_plane,
            figure_attribute='figure_2d',
            canvas_attribute='canvas_2d'
        )

    def add_tab_gaussian_and_radial(self):
        """Onglet 3 : Affichage du champ radial et de l'amplitude simulée avec sélecteur de plan."""
        self.add_tab_with_plane_selector(
            title="Champ Amplitude et Vectoriel",
            update_selector_function=self.update_plane_selector_3d_values,
            update_plot_function=self.update_tab_gaussian_and_radial,
            figure_attribute='figure_gaussian',
            canvas_attribute='canvas_gaussian'
        )

    def add_tab_with_plane_selector(self, title, update_selector_function, update_plot_function, figure_attribute, canvas_attribute):
        """Ajoute un onglet avec un sélecteur de plan et une figure."""
        tab = QWidget()
        layout = QVBoxLayout()
        tab.setLayout(layout)

        # Créer les sélecteurs
        selector_layout, plane_selector, value_selector, plane_label = self.create_plane_and_value_selectors(
            update_selector_function, update_plot_function
        )
        selector_layout.addWidget(plane_label, 0, Qt.AlignLeft)
        selector_layout.addWidget(plane_selector, 1, Qt.AlignLeft)
        
        # Ajouter le bouton pour modifier la résolution
        resolution_button = QPushButton("Modifier Résolution")
        resolution_button.clicked.connect(self.open_resolution_dialog)
        selector_layout.addWidget(resolution_button, alignment=Qt.AlignLeft)
        
        layout.addLayout(selector_layout)

        # Initialiser les attributs pour les sélecteurs
        if title == "Champ dans le plan 2D":
            self.plane_selector_2d = plane_selector
            self.value_selector_2d = value_selector  # Initialiser l'attribut de sélecteur 2D
        else:
            self.plane_selector_3d = plane_selector
            self.value_selector_3d = value_selector  # Initialiser l'attribut de sélecteur 3D

        # Ajouter une figure et un canvas
        figure = plt.figure(figsize=(12, 6))
        canvas = FigureCanvas(figure)
        setattr(self, figure_attribute, figure)
        setattr(self, canvas_attribute, canvas)
        layout.addWidget(canvas)

        # Initialiser les valeurs du sélecteur et le tracé
        update_selector_function()
        update_plot_function()

        # Ajouter l'onglet
        self.tabs.addTab(tab, title)

    def open_resolution_dialog(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Modifier la Résolution")

        layout = QVBoxLayout()
        form_layout = QFormLayout()

        resolution_input = QLineEdit()
        resolution_input.setPlaceholderText("Entrez la nouvelle résolution")

        form_layout.addRow("Résolution :", resolution_input)

        ok_button = QPushButton("OK")
        ok_button.clicked.connect(lambda: self.set_resolution(resolution_input.text(), dialog))

        layout.addLayout(form_layout)
        layout.addWidget(ok_button)

        dialog.setLayout(layout)
        dialog.exec_()


        
    def create_plane_and_value_selectors(self, plane_callback, value_callback):
        """Crée les sélecteurs pour le plan et la valeur."""
        selector_layout = QHBoxLayout()
        selector_layout.setContentsMargins(0, 0, 0, 0)

        # Sélecteur de plan
        plane_label = QLabel("Plan:")
        plane_selector = QComboBox()
        plane_selector.addItems(['x', 'y', 'z'])
        plane_selector.currentTextChanged.connect(plane_callback)

        # Sélecteur de valeur
        value_label = QLabel("Valeur:")
        value_selector = QComboBox()
        value_selector.currentTextChanged.connect(value_callback)

        selector_layout.addWidget(plane_label, alignment=Qt.AlignLeft)
        selector_layout.addWidget(plane_selector, alignment=Qt.AlignLeft)
        selector_layout.addWidget(value_label, alignment=Qt.AlignLeft)
        selector_layout.addWidget(value_selector, alignment=Qt.AlignLeft)

        return selector_layout, plane_selector, value_selector, plane_label

    def update_plane_selector_2d_values(self):
        """Met à jour les valeurs disponibles dans le sélecteur pour l'onglet 2D."""
        self.update_plane_selector_values(self.value_selector_2d, self.plane_selector_2d.currentText())

    def update_plane_selector_3d_values(self):
        """Met à jour les valeurs disponibles dans le sélecteur pour l'onglet 3D."""
        self.update_plane_selector_values(self.value_selector_3d, self.plane_selector_3d.currentText())

    def update_plane_selector_values(self, selector, plane):
        """Met à jour les valeurs disponibles dans le sélecteur de valeurs."""
        selector.blockSignals(True)
        selector.clear()
        unique_values = sorted(np.unique([p[plane] for p in self.simulation.points_haute_resolution])) if plane in ['x', 'y', 'z'] else []
        selector.addItems([f"{v:.5f}" for v in unique_values])
        selector.blockSignals(False)

        if unique_values:
            selector.setCurrentIndex(len(unique_values) // 2)

    def update_tab_2d_plane(self):
        """Met à jour les graphiques de l'onglet Champ dans le plan 2D."""
        plane = self.plane_selector_2d.currentText()
        if self.value_selector_2d.count() == 0:
            return

        try:
            value = float(self.value_selector_2d.currentText())
        except ValueError:
            return

        filtered_points = self.filter_points(plane, value)
        if not filtered_points:
            print(f"Aucun point trouvé pour le plan {plane}={value}")
            return

        coord1, coord2, H_total, H_component1_norm, H_component2_norm = self.prepare_plot_data(filtered_points, plane)

        coord1_grid, coord2_grid = np.meshgrid(np.unique(coord1), np.unique(coord2))

        H_total_grid = griddata((coord1, coord2), H_total, (coord1_grid, coord2_grid), method='linear', fill_value=0)
        H_component1_norm_grid = griddata((coord1, coord2), H_component1_norm, (coord1_grid, coord2_grid), method='linear', fill_value=0)
        H_component2_norm_grid = griddata((coord1, coord2), H_component2_norm, (coord1_grid, coord2_grid), method='linear', fill_value=0)

        self.figure_2d.clear()
        gs = self.figure_2d.add_gridspec(1, 3, width_ratios=[6, 0.4, 6])

        ax1 = self.figure_2d.add_subplot(gs[0, 0])
        contour = ax1.contourf(coord1_grid, coord2_grid, H_total_grid, levels=20, cmap='viridis')
        ax1.set_title(f"Norme du champ |H| ({plane})")
        ax1.set_xlabel(f'{plane} (m)')

        cbar_ax = self.figure_2d.add_subplot(gs[0, 1])
        self.figure_2d.colorbar(contour, cax=cbar_ax, label='|H| (A/m)')

        ax2 = self.figure_2d.add_subplot(gs[0, 2])
        quiver = ax2.quiver(coord1_grid, coord2_grid, H_component1_norm_grid, H_component2_norm_grid, color='red', scale=12)
        ax2.set_title(f"Direction du champ magnétique sur le plan {plane}")
        ax2.set_xlabel(f'{plane} (m)')
        ax2.set_aspect('equal')

        self.canvas_2d.draw()

    def update_tab_gaussian_and_radial(self):
        """Met à jour les graphiques de l'onglet Champ Amplitude et Vectoriel."""
        plane = self.plane_selector_3d.currentText()
        if self.value_selector_3d.count() == 0:
            return 

        try:
            value = float(self.value_selector_3d.currentText())
        except ValueError:
            return

        filtered_points = self.filter_points(plane, value)
        if not filtered_points:
            print(f"Aucun point trouvé pour le plan {plane}={value}")
            return

        coord1, coord2, H_total, H_component1_norm, H_component2_norm = self.prepare_plot_data(filtered_points, plane)

        coord1_grid, coord2_grid = np.meshgrid(np.unique(coord1), np.unique(coord2))

        H_total_grid = griddata((coord1, coord2), H_total, (coord1_grid, coord2_grid), method='linear', fill_value=0)
        H_component1_norm_grid = griddata((coord1, coord2), H_component1_norm, (coord1_grid, coord2_grid), method='linear', fill_value=0)
        H_component2_norm_grid = griddata((coord1, coord2), H_component2_norm, (coord1_grid, coord2_grid), method='linear', fill_value=0)

        self.figure_gaussian.clear()
        gs = self.figure_gaussian.add_gridspec(1, 2, width_ratios=[1, 1])

        ax1 = self.figure_gaussian.add_subplot(gs[0, 0], projection='3d')
        surf = ax1.plot_surface(coord1_grid, coord2_grid, H_total_grid, cmap='viridis', edgecolor='k', alpha=0.8)
        ax1.set_title(f'Amplitude du champ magnétique |H| ({plane})')
        ax1.set_xlabel(f'{plane} (m)')
        ax1.set_zlabel('Amplitude |H| (A/m)')
        self.figure_gaussian.colorbar(surf, ax=ax1, shrink=0.5, aspect=10)

        ax2 = self.figure_gaussian.add_subplot(gs[0, 1])
        quiver = ax2.quiver(coord1_grid, coord2_grid, H_component1_norm_grid, H_component2_norm_grid, scale=12)
        ax2.set_title(f"Champ vectoriel sur le plan {plane}")
        ax2.set_xlabel(f'{plane} (m)')
        ax2.set_aspect('equal')

        self.canvas_gaussian.draw()

    def filter_points(self, plane, value, epsilon=1e-5):
        """Filtre les points selon le plan et la valeur donnée."""
        return [p for p in self.simulation.points_haute_resolution if abs(p[plane] - value) < epsilon]

    def prepare_plot_data(self, filtered_points, plane):
        """Prépare les données de tracé pour les graphiques."""
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

        with np.errstate(divide='ignore', invalid='ignore'):
            H_component1_normalized = np.where(H_total != 0, H_component1 / H_total, 0)
            H_component2_normalized = np.where(H_total != 0, H_component2 / H_total, 0)

        return coord1, coord2, H_total, H_component1_normalized, H_component2_normalized
    
# Lancement de l'application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_window = MagneticFieldApp()
    main_window.show()
    sys.exit(app.exec_())