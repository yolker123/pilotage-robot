import sys
import numpy as np
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QTabWidget, QVBoxLayout, QWidget,
    QHBoxLayout, QLabel, QComboBox, QPushButton, QFileDialog, QRadioButton, QButtonGroup
)
from PyQt5.QtWidgets import QDialog, QLineEdit, QPushButton, QVBoxLayout, QFormLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtCore import Qt, QTimer

from scipy.interpolate import griddata

from MagneticFieldSimulation import MagneticFieldSimulation


class MagneticFieldApp(QWidget):
    def __init__(self):
        super().__init__()
        self.columns = None
        self.setWindowTitle("Simulation du Champ Magnétique")

        # Layout principal
        layout = QVBoxLayout(self)
        self.setLayout(layout)

        # Ajout des onglets
        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)

        self.simulation = MagneticFieldSimulation(resolution=1)
        print(self.simulation.resultats)
        print("ok")

        self.lines = None  # Store lines from file
        self.line_index = 0  # Line index
        self.timer = QTimer()
        self.initUI()
        if self.simulation.resultats:
            self.simulation.augmenter_resolution(self.simulation.resultats)  # Augmenter la résolution

    def select_file(self):
        # Ouvrir une boîte de dialogue pour sélectionner un fichier
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog  # Utiliser un dialogue natif selon votre système
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File", "", "All Files (*);;Text Files (*.txt)",
                                                   options=options)
        if file_path:
            self.simulation.file_path = file_path
            self.start_reading_file()  #

    def start_reading_file(self):
        # Initialisation du temporisateur pour lire les lignes à intervalles
        self.timer.timeout.connect(self.read_next_line)

        # Vérifier que le fichier peut être ouvert
        try:
            with open(self.simulation.file_path, 'r') as file:
                lines = file.readlines()
        except FileNotFoundError:
            print("Erreur : Le fichier ne peut pas être trouvé.")
            return
        except IOError:
            print("Erreur : Impossible d'ouvrir le fichier.")
            return

        # Vérifier le contenu du fichier
        if not lines:
            print("Erreur : Le fichier est vide.")
            return

        # Afficher les lignes lues (pour vérification)
        print("Lignes lues du fichier :")
        for line in lines:
            print(line.strip())

        # Extraire l'en-tête pour les colonnes
        self.columns = lines[0].strip().split(',')
        print(f"Noms des colonnes : {self.columns}")

        # Stocker le reste des lignes
        self.lines = lines[1:]  # Les données sans l'en-tête
        self.line_index = 0  # Démarrer à la première ligne de données

        # Démarrer le temporisateur
        self.timer.start(2000)

    def read_next_line(self):
        if self.line_index < len(self.lines):
            line = self.lines[self.line_index].strip()
            if line:
                self.simulation.read_file_and_calculate_point(line, self.columns)
                self.update_all_graphs()
            self.line_index += 1
        else:
            self.timer.stop()

    def initUI(self):
        self.add_tab_3d_vectors()
        self.add_tab_2d_plane()
        self.add_tab_gaussian_and_radial()

    def empty_graph(self):
        self.simulation.resultats = []
        self.simulation.points_haute_resolution = []
        self.update_all_graphs()

    def set_resolution(self, resolution_value, dialog, algorithm):
        try:
            resolution_value = int(resolution_value)
            self.simulation.resolution = resolution_value  # Modifier la résolution de la simulation
            self.simulation.augmenter_resolution(self.simulation.resultats)  # Augmenter la résolution
            self.update_all_graphs()  # Appeler une méthode pour rafraîchir les graphiques

            self.update_plane_selector_2d_values()
            self.update_plane_selector_3d_values()

            self.plane_selector_2d.currentTextChanged.emit(self.plane_selector_2d.currentText())
            self.plane_selector_3d.currentTextChanged.emit(self.plane_selector_3d.currentText())

            dialog.accept()
            print(f"Résolution modifiée à : {resolution_value}")
        except ValueError:
            print("Veuillez entrer un nombre entier valide.")

    def update_all_graphs(self):
        self.update_tab_2d_plane()  # Mettre à jour les graphiques 2D
        self.update_tab_gaussian_and_radial()  # Mettre à jour les graphiques 3D
        self.plot_3d_vectors()  # Mettre à jour les vecteurs 3D
        self.update_plane_selector_2d_values()
        self.update_plane_selector_3d_values()
        self.plane_selector_2d.currentTextChanged.emit(self.plane_selector_2d.currentText())
        self.plane_selector_3d.currentTextChanged.emit(self.plane_selector_3d.currentText())

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

        button_layout = QHBoxLayout()

        # Ajouter les boutons au layout horizontal
        file_button = QPushButton("Importer une mesure")
        file_button.clicked.connect(self.select_file)
        button_layout.addWidget(file_button)

        resolution_button = QPushButton("Modifier Résolution")
        resolution_button.clicked.connect(self.open_resolution_dialog)
        button_layout.addWidget(resolution_button)

        # Aligner le layout des boutons à droite
        button_layout.addStretch(1)  # Ajoute un espacement flexible à gauche des boutons

        # Ajouter le layout des boutons au layout principal (vertical)
        layout.addLayout(button_layout)
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
            self.ax_3d.quiver(x, y, z, Hx, Hy, Hz, color='b', length=0.1, normalize=True)
        # Tracer les vecteurs interpolés
        for point in points_interpolés:
            x, y, z = point['x'], point['y'], point['z']
            Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
            self.ax_3d.quiver(x, y, z, Hx, Hy, Hz, color='r', length=0.1, normalize=True)

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

    def add_tab_with_plane_selector(self, title, update_selector_function, update_plot_function, figure_attribute,
                                    canvas_attribute):
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
        file_button = QPushButton("Importer une mesure")
        file_button.clicked.connect(self.select_file)
        selector_layout.addWidget(file_button, alignment=Qt.AlignLeft)
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

        # Create radio buttons for linearity selection
        linear_button = QRadioButton("Linéaire")
        nonlinear_button = QRadioButton("Non Linéaire")
        form_layout.addRow(linear_button, nonlinear_button)

        # Group the radio buttons
        button_group = QButtonGroup()
        button_group.addButton(linear_button)
        button_group.addButton(nonlinear_button)

        # Add to the main layout
        layout.addLayout(form_layout)

        # Button to close the dialog
        ok_button = QPushButton("OK")
        ok_button.setEnabled(False)  # Initially disable the button
        layout.addWidget(ok_button)

        # Enable the OK button only when a radio button is selected
        def enable_ok_button():
            if linear_button.isChecked() or nonlinear_button.isChecked():
                ok_button.setEnabled(True)
            else:
                ok_button.setEnabled(False)

        linear_button.toggled.connect(enable_ok_button)
        nonlinear_button.toggled.connect(enable_ok_button)

        def on_ok_clicked():
            algorithm = "non-linear" if nonlinear_button.isChecked() else "linear"
            self.set_resolution(resolution_input.text(), dialog, algorithm)
            dialog.accept()

        ok_button.clicked.connect(on_ok_clicked)

        dialog.setLayout(layout)
        dialog.exec_()

    def create_plane_and_value_selectors(self, plane_callback, value_callback):
        """Crée les sélecteurs pour le plan et la valeur."""
        selector_layout = QHBoxLayout()
        selector_layout.setContentsMargins(1, 1, 1, 1)

        # Sélecteur de plan
        plane_label = QLabel("Plan:")
        plane_selector = QComboBox()
        plane_selector.addItems(['x', 'y', 'z'])
        plane_selector.currentTextChanged.connect(plane_callback)
        plane_selector.setMinimumWidth(75)  # Réglez la largeur minimale si besoin


        # Sélecteur de valeur
        value_label = QLabel("Valeur:")
        value_selector = QComboBox()
        value_selector.currentTextChanged.connect(value_callback)
        value_selector.setMinimumWidth(75)  # Ajustez la largeur minimale si nécessaire

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

        # rajouter les points resultats qui ne sont pas dans point haute resolution, fait des arrondis il peut y avoir un petit ecart
        # 2 points bugués
        for p in self.simulation.resultats:
            if p not in self.simulation.points_haute_resolution:
                self.simulation.points_haute_resolution.append(p)

        if plane in ['x', 'y', 'z']:
            # Extraire les valeurs uniques arrondies
            raw_values = [p[plane] for p in self.simulation.points_haute_resolution]
            unique_values = sorted(np.unique([round(v, 5) for v in raw_values]))
        else:
            unique_values = []
        print(unique_values)
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
        if not filtered_points or len(filtered_points) < 3:  # Ensure enough points for a plane
            print(f"Aucun point trouvé pour le plan {plane}={value}, ou pas assez de points.")
            return

        coord1, coord2, H_total, H_component1_norm, H_component2_norm = self.prepare_plot_data(filtered_points, plane)

        # If there are not enough points to interpolate smoothly, do a scatter plot
        if len(np.unique(coord1)) < 2 or len(np.unique(coord2)) < 2:
            print("Insufficient unique points for grid interpolation; using scatter plot.")
            self.figure_2d.clear()
            ax = self.figure_2d.add_subplot(111)
            scatter = ax.scatter(coord1, coord2, c=H_total, cmap='viridis', edgecolor='k')
            ax.quiver(coord1, coord2, H_component1_norm, H_component2_norm, color='red', scale=12)
            ax.set_title(f"Points sur le plan {plane}")
            if plane == "x":
                ax.set_xlabel(f'y (m)')
                ax.set_ylabel(f'z (m)')
            if plane == "y":
                ax.set_xlabel(f'x (m)')
                ax.set_ylabel(f'z (m)')
            if plane == "z":
                ax.set_xlabel(f'x (m)')
                ax.set_ylabel(f'y (m)')
            ax.set_ylabel('Other axis (m)')  # Change accordingly
            self.figure_2d.colorbar(scatter, ax=ax, label='|H| (A/m)')
            self.canvas_2d.draw()
            return

        # Proceed as usual if data is sufficient
        coord1_grid, coord2_grid = np.meshgrid(np.unique(coord1), np.unique(coord2))
        H_total_grid = griddata((coord1, coord2), H_total, (coord1_grid, coord2_grid), method='linear', fill_value=0)
        H_component1_norm_grid = griddata((coord1, coord2), H_component1_norm, (coord1_grid, coord2_grid),
                                          method='linear', fill_value=0)
        H_component2_norm_grid = griddata((coord1, coord2), H_component2_norm, (coord1_grid, coord2_grid),
                                          method='linear', fill_value=0)

        self.figure_2d.clear()
        gs = self.figure_2d.add_gridspec(1, 3, width_ratios=[6, 0.4, 6])
        ax1 = self.figure_2d.add_subplot(gs[0, 0])
        contour = ax1.contourf(coord1_grid, coord2_grid, H_total_grid, levels=20, cmap='viridis')
        ax1.set_title(f"Norme du champ |H| ({plane})")
        if plane == "x":
            ax1.set_xlabel(f'y (m)')
            ax1.set_ylabel(f'z (m)')
        if plane == "y":
            ax1.set_xlabel(f'x (m)')
            ax1.set_ylabel(f'z (m)')
        if plane == "z":
            ax1.set_xlabel(f'x (m)')
            ax1.set_ylabel(f'y (m)')
        cbar_ax = self.figure_2d.add_subplot(gs[0, 1])
        self.figure_2d.colorbar(contour, cax=cbar_ax, label='|H| (A/m)')
        ax2 = self.figure_2d.add_subplot(gs[0, 2])
        quiver = ax2.quiver(coord1_grid, coord2_grid, H_component1_norm_grid, H_component2_norm_grid, color='red',
                            scale=12)
        ax2.set_title(f"Direction du champ magnétique sur le plan {plane}")
        if plane == "x":
            ax2.set_xlabel(f'y (m)')
            ax2.set_ylabel(f'z (m)')
        if plane == "y":
            ax2.set_xlabel(f'x (m)')
            ax2.set_ylabel(f'z (m)')
        if plane == "z":
            ax2.set_xlabel(f'x (m)')
            ax2.set_ylabel(f'y (m)')
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
        if not filtered_points or len(filtered_points) < 3:
            print(f"Aucun point trouvé pour le plan {plane}={value}, ou pas assez de points.")
            return

        coord1, coord2, H_total, H_component1_norm, H_component2_norm = self.prepare_plot_data(filtered_points, plane)

        # Check if there are enough unique points for interpolation
        if len(np.unique(coord1)) < 2 or len(np.unique(coord2)) < 2:
            print("Insufficient unique points for grid interpolation; using scatter plot.")
            self.figure_gaussian.clear()
            ax1 = self.figure_gaussian.add_subplot(111, projection='3d')
            ax1.scatter(coord1, coord2, H_total, c=H_total, cmap='viridis', edgecolor='k', alpha=0.8)
            ax1.set_title(f"Points sur le plan {plane}")
            if plane == "x":
                ax1.set_xlabel(f'y (m)')
                ax1.set_ylabel(f'z (m)')
            if plane == "y":
                ax1.set_xlabel(f'x (m)')
                ax1.set_ylabel(f'z (m)')
            if plane == "z":
                ax1.set_xlabel(f'x (m)')
                ax1.set_ylabel(f'y (m)')
            ax1.set_zlabel('Amplitude |H| (A/m)')
            self.canvas_gaussian.draw()
            return

        # Proceed with grid interpolation if data is sufficient
        coord1_grid, coord2_grid = np.meshgrid(np.unique(coord1), np.unique(coord2))
        H_total_grid = griddata((coord1, coord2), H_total, (coord1_grid, coord2_grid), method='linear', fill_value=0)
        H_component1_norm_grid = griddata((coord1, coord2), H_component1_norm, (coord1_grid, coord2_grid),
                                          method='linear', fill_value=0)
        H_component2_norm_grid = griddata((coord1, coord2), H_component2_norm, (coord1_grid, coord2_grid),
                                          method='linear', fill_value=0)

        self.figure_gaussian.clear()
        gs = self.figure_gaussian.add_gridspec(1, 2, width_ratios=[1, 1])

        ax1 = self.figure_gaussian.add_subplot(gs[0, 0], projection='3d')
        surf = ax1.plot_surface(coord1_grid, coord2_grid, H_total_grid, cmap='viridis', edgecolor='k', alpha=0.8)
        ax1.set_title(f'Amplitude du champ magnétique |H| ({plane})')
        if plane == "x":
            ax1.set_xlabel(f'y (m)')
            ax1.set_ylabel(f'z (m)')
        if plane == "y":
            ax1.set_xlabel(f'x (m)')
            ax1.set_ylabel(f'z (m)')
        if plane == "z":
            ax1.set_xlabel(f'x (m)')
            ax1.set_ylabel(f'y (m)')
        ax1.set_zlabel('Amplitude |H| (A/m)')
        self.figure_gaussian.colorbar(surf, ax=ax1, shrink=0.5, aspect=10)

        ax2 = self.figure_gaussian.add_subplot(gs[0, 1])
        quiver = ax2.quiver(coord1_grid, coord2_grid, H_component1_norm_grid, H_component2_norm_grid, scale=12)
        ax2.set_title(f"Champ vectoriel sur le plan {plane}")
        if plane == "x":
            ax2.set_xlabel(f'y (m)')
            ax2.set_ylabel(f'z (m)')
        if plane == "y":
            ax2.set_xlabel(f'x (m)')
            ax2.set_ylabel(f'z (m)')
        if plane == "z":
            ax2.set_xlabel(f'x (m)')
            ax2.set_ylabel(f'y (m)')
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
