"""!
 * @file        GraphicsTab.py
 * @brief       Class corresponding to the graphics tab containing these display and update methods
 * @author      DEVAUX Baptiste | VOLPELLIERE Anthony | MULLOT Agathe
 * @version     0.1
 * @date        2025
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QTabWidget, QVBoxLayout, QWidget,
    QHBoxLayout, QLabel, QComboBox, QPushButton, QFileDialog, QRadioButton, QButtonGroup, QSlider, QCheckBox
)
from PyQt5.QtWidgets import QDialog, QLineEdit, QPushButton, QVBoxLayout, QFormLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtCore import Qt, QTimer

from scipy.interpolate import griddata

from MagneticFieldCalculation import MagneticFieldCalculation
import math


class GraphicsTab(QWidget):
    def __init__(self):
        """
        Initialise le widget d'onglet graphique avec des méthodes pour l'affichage et les mises à jour.
        """
        super().__init__()
        self.normalize_button_3D = None
        self.normalize_button_2D = None
        self.normalize_buttons_2D = []
        self.scale_input_3d = None
        self.scale_input_2d = None
        self.setWindowTitle("Simulation du Champ Magnétique")

        # Layout principal
        layout = QVBoxLayout(self)
        self.setLayout(layout)

        # Ajout des onglets
        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)

        # Variable
        self.vector_scale_2d = 12
        self.vector_length_3d = 1

        self.simulation = MagneticFieldCalculation(resolution=1)
        self.initUI()

    def initUI(self):
        """
        Initialise l'interface utilisateur en créant différents onglets.
        """
        self.create_tab_3d_vectors()
        self.add_tab_2d_plane()
        self.add_tab_gaussian_and_radial()

    def select_file(self):
        """
        Ouvre une boîte de dialogue pour sélectionner un fichier, puis commence à lire le fichier sélectionné.
        """
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select File",
            "",
            "CSV Files (*.csv);;All Files (*)",
            options=options
        )

        if file_path:
            hRobot = self.get_hrobot_from_filename(file_path)
            print(f"Valeur extraite : {hRobot}")  # Affiche la valeur extraite
            self.simulation.hRobot = hRobot
            self.start_reading_file(file_path)
            formated_file_path = self.format_file_path(file_path)
            self.label_file_path.setText(formated_file_path)
            self.label_file_path.setToolTip(file_path)
            # self.show_file_path(file_path)
            self.update_all_graphs()

    def format_file_path(self, path: str) -> str:
        # If the path length is less than or equal to 150 characters, return the original path.
        if len(path) <= 150:
            return path

        # Normalize and split the path to get individual parts
        parts = [p for p in os.path.normpath(path).split(os.path.sep) if p]

        # If there are no parts, return the original path
        if not parts:
            return path

        # Determine the first directory
        first_dir = parts[0]
        # Get the file name from the path
        file_name = os.path.basename(path)

        # Construct and return the shortened path string
        return f"{first_dir}/ ... /{file_name}"

    def get_hrobot_from_filename(self, file_path):
        """
        Extrait la valeur hRobot d'un nom de fichier donné.

        Parameters:
        file_path (str): Chemin du fichier CSV.

        Returns:
        str: Valeur extraite après le dernier underscore '_' dans le nom du fichier.
        """
        filename = os.path.basename(file_path)  # Récupérer le nom du fichier sans le chemin
        name_without_extension = filename.split('.csv')[0]  # Enlever l'extension .csv
        value = name_without_extension.rsplit('_', 1)[-1]  # Récupérer la partie après le dernier '_'
        return value

    def start_reading_file(self, file_path):
        """
        Lit un fichier et procède à l'analyse des mesures contenues.

        Parameters:
        file_path (str): Chemin du fichier CSV.
        """
        lines = None
        try:
            with open(file_path, 'r') as file:
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

        # Extraire l'en-tête pour les colonnes
        columns = lines[0].strip().split(',')
        print(f"Noms des colonnes : {columns}")

        lines = lines[1:]  # Enlever la ligne d'en-tête

        for line in lines:
            line = line.strip()
            if line:
                self.simulation.read_file_and_calculate_point(line, columns)




    def set_resolution(self, resolution_value, dialog, algorithm):
        """
        Définit la résolution de simulation, en utilisant l'algorithme spécifié.

        Parameters:
        resolution_value (str): Nouvelle valeur de résolution entrée par l'utilisateur.
        dialog (QDialog): Dialogue de saisie de résolution.
        algorithm (str): Type d'algorithme à utiliser ("linear" ou "non-linear").
        """
        try:
            resolution_value = int(resolution_value)
            self.simulation.resolution = resolution_value  # Modifier la résolution de la simulation
            if algorithm == "linear":
                self.simulation.augmenter_resolution(self.simulation.measuredPoints, algorithm)  # Augmenter la résolution
            else:
                points_proches = self.simulation.selectionner_points_proches()
                I_moyen = self.simulation.moyenne_I(points_proches)
                print(I_moyen)
                self.simulation.augmenter_resolution(self.simulation.measuredPoints, algorithm, I_moyen)
            print("calcul fini")
            self.update_all_graphs()  # Appeler une méthode pour rafraîchir les graphiques

            dialog.accept()
            print(f"Résolution modifiée à : {resolution_value}")
        except ValueError:
            print("Veuillez entrer un nombre entier valide.")

    def update_all_graphs(self):
        """
        Met à jour tous les graphiques de chaque onglet.
        """
        # Premièrement, on s'assure que les sélecteurs affichent les dernières valeurs disponibles
        self.update_plane_selector_2d_values()
        self.update_plane_selector_3d_values()

        # Emission manuelle du signal de changement pour forcer la mise à jour
        if self.plane_selector_2d.count() > 0:
            self.plane_selector_2d.currentTextChanged.emit(self.plane_selector_2d.currentText())
        if self.plane_selector_3d.count() > 0:
            self.plane_selector_3d.currentTextChanged.emit(self.plane_selector_3d.currentText())

        # Puis on met à jour les graphiques
        self.update_tab_2d_plane()
        self.update_tab_gaussian_and_radial()
        self.plot_3d_vectors()

    def clear_all_graphs(self):
        """
        Efface tous les points de mesure et réinitialise les graphiques.
        """
        self.simulation.points_haute_resolution = []
        self.simulation.measuredPoints = []
        self.label_file_path.setText("")
        self.update_all_graphs()

    def clear_interpolated_points(self):
        """
        Efface tous les points interpolés et réinitialise les graphiques.
        """
        self.simulation.points_haute_resolution = []
        self.label_file_path.setText("")
        self.update_all_graphs()

    def create_tab_3d_vectors(self):
        """
        Crée l'onglet d'affichage 3D des vecteurs.
        """
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

        self.label_file_path = QLabel("")
        button_layout.addWidget(self.label_file_path)

        # Aligner le layout des boutons à droite
        button_layout.addStretch(1)  # Ajoute un espacement flexible à gauche des boutons

        clear_button = QPushButton("Effacer")
        clear_button.clicked.connect(self.clear_all_graphs)
        button_layout.addWidget(clear_button)

        save_button = QPushButton("Sauvegarder")
        save_button.clicked.connect(self.saveToFile)
        button_layout.addWidget(save_button)

        # Ajouter le layout des boutons au layout principal (vertical)
        layout.addLayout(button_layout)

        scale_layout = self.create_scale_layout("3D")
        layout.addLayout(scale_layout)

        self.normalize_button_3D = QCheckBox("Normaliser les vecteurs")
        self.normalize_button_3D.setChecked(False)  # Valeur par défaut : False
        layout.addWidget(self.normalize_button_3D)
        self.normalize_button_3D.stateChanged.connect(self.update_all_graphs)

        layout.addWidget(self.canvas_3d)
        self.tabs.addTab(self.tab_3d, "Vecteurs 3D")
        self.add_vector_filter_layout(layout)

    def create_scale_layout(self, dimension):
        """
        Crée le layout en fonction de la dimension en reprenant la bonne fonction
        """
        scale_layout = QHBoxLayout()
        scale_label = QLabel("Vector scale :")
        scale_input = QLineEdit()
        scale_value_label = QLabel()
        scale_button = QPushButton("Valider")

        info_label = QLabel("3D : Higher is bigger / 2D Higher is smaller :")
        scale_label.setFixedSize(75, 25)
        scale_input.setFixedSize(50, 25)  # Définit une taille fixe pour éviter qu'il prenne trop de place
        scale_value_label.setFixedSize(50, 25)
        scale_button.setFixedSize(70, 25)
        info_label.setFixedSize(500, 25)
        if dimension == "2D":
            scale_input.setText(str(self.vector_scale_2d))
            scale_value_label.setText(str(self.vector_scale_2d))
            scale_button.clicked.connect(lambda: self.update_scale_2d(float(scale_input.text())))
            self.scale_input_2d = scale_input
        else:
            scale_input.setText(str(self.vector_length_3d))
            scale_value_label.setText(str(self.vector_length_3d))
            scale_button.clicked.connect(lambda: self.update_scale_3d(float(scale_input.text())))
            self.scale_input_3d = scale_input

        scale_layout.addWidget(scale_label)
        scale_layout.addWidget(scale_input)
        scale_layout.addWidget(scale_button)
        scale_layout.addWidget(info_label)

        scale_layout.setAlignment(Qt.AlignLeft)
        return scale_layout

    def update_scale_2d(self, value):
        """
        Met à jour l'échelle des graphiques 2D.
        """
        self.vector_scale_2d = value
        self.update_all_graphs()

    def update_scale_3d(self, value):
        """
        Met à jour l'échelle des graphiques 3D.
        """
        self.vector_length_3d = value
        self.update_all_graphs()

    def add_vector_filter_layout(self, layout):
        """
        Ajoute un ensemble de contrôles de filtre pour la visualisation des vecteurs 3D.
        """
        filter_layout = QHBoxLayout()

        # Filter type selector
        filter_label = QLabel("Filter by:")
        self.filter_type_selector = QComboBox()
        self.filter_type_selector.addItems(['None', '|H|', 'Hx', 'Hy', 'Hz'])
        filter_layout.addWidget(filter_label)
        filter_layout.addWidget(self.filter_type_selector)

        # Min value input
        min_label = QLabel("Min:")
        self.min_value_input = QLineEdit()
        self.min_value_input.setPlaceholderText("Minimum value")
        filter_layout.addWidget(min_label)
        filter_layout.addWidget(self.min_value_input)

        # Max value input
        max_label = QLabel("Max:")
        self.max_value_input = QLineEdit()
        self.max_value_input.setPlaceholderText("Maximum value")
        filter_layout.addWidget(max_label)
        filter_layout.addWidget(self.max_value_input)

        # Apply filter button
        apply_filter_button = QPushButton("Apply Filter")
        apply_filter_button.clicked.connect(self.apply_vector_filter)
        filter_layout.addWidget(apply_filter_button)

        reset_filter_button = QPushButton("Reset Filter")
        reset_filter_button.clicked.connect(self.reset_vector_filter)
        filter_layout.addWidget(reset_filter_button)

        layout.addLayout(filter_layout)

    def apply_vector_filter(self):
        """
        Applique un filtre au tracé vectoriel 3D basé sur les critères sélectionnés.
        """
        filter_type = self.filter_type_selector.currentText()

        # Get min and max values, defaulting to None if not provided
        try:
            min_val = float(self.min_value_input.text()) if self.min_value_input.text() else None
        except ValueError:
            min_val = None

        try:
            max_val = float(self.max_value_input.text()) if self.max_value_input.text() else None
        except ValueError:
            max_val = None

        # Filter points based on selected criteria
        if filter_type == 'None':
            filtered_points = self.simulation.points_haute_resolution
        else:
            filtered_points = self.filter_vector_points(
                self.simulation.points_haute_resolution,
                filter_type,
                min_val,
                max_val
            )
        # self.clear_all_graphs()
        # for p in filtered_points:
        #     print(p)
        self.simulation.points_haute_resolution = filtered_points

        self.update_all_graphs()

    def filter_vector_points(self, points, filter_type, min_val=None, max_val=None):
        """
        Filter points based on vector magnitude or component.

        Args:
            points (list): List of point dictionaries
            filter_type (str): Type of filter ('|H|', 'Hx', 'Hy', 'Hz')
            min_val (float, optional): Minimum value for filtering
            max_val (float, optional): Maximum value for filtering

        Returns:
            list: Filtered list of points
        """

        def point_meets_criteria(point):
            if filter_type == '|H|':
                value = np.sqrt(point['Hx'] ** 2 + point['Hy'] ** 2 + point['Hz'] ** 2)
            else:
                value = point[filter_type]
            # Check against min and max values
            if min_val is not None and value < min_val:
                return False
            if max_val is not None and value > max_val:
                return False
            print(value , "true")
            return True

        for point in points:
            if point_meets_criteria(point):
                point["display"] = True
            else:
                point["display"] = False

        return points

    def reset_vector_filter(self):
        for point in self.simulation.points_haute_resolution:
            point["display"] = True
        tmp = self.simulation.points_haute_resolution
        self.clear_interpolated_points()
        self.simulation.points_haute_resolution = tmp
        self.update_all_graphs()

    def plot_3d_vectors(self):
        """
        Trace les vecteurs 3D sur l'onglet correspondant.
        """
        self.ax_3d.clear()  # Effacer les anciens vecteurs

        # Données pour les vecteurs
        points = self.simulation.measuredPoints
        points_interpolés = self.simulation.points_haute_resolution

        # Tracer les vecteurs des points originaux
        for point in points:
            if point["display"]:
                x, y, z = point['x'], point['y'], point['z']
                Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
                self.ax_3d.quiver(x, y, z, Hx, Hy, Hz, color='b', length=self.vector_length_3d+1, normalize=self.normalize_button_3D.isChecked())
        # Tracer les vecteurs interpolés
        for point in points_interpolés:
            if point["display"]:
                x, y, z = point['x'], point['y'], point['z']
                Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
                self.ax_3d.quiver(x, y, z, Hx, Hy, Hz, color='r', length=self.vector_length_3d, normalize=self.normalize_button_3D.isChecked())

        # Configurer les axes
        self.ax_3d.set_title("Vecteurs 3D du champ magnétique")
        self.ax_3d.set_xlabel('X (mm)')
        self.ax_3d.set_ylabel('Y (mm)')
        self.ax_3d.set_zlabel('Z (mm)')

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
        scale_layout = self.create_scale_layout("2D")
        layout.addLayout(scale_layout)

        self.normalize_button_2D = QCheckBox("Normaliser les vecteurs")
        self.normalize_button_2D.setChecked(False)  # Valeur par défaut : False
        self.normalize_button_2D.stateChanged.connect(self.sync_checkboxes_2D)
        layout.addWidget(self.normalize_button_2D)

        self.normalize_buttons_2D.append(self.normalize_button_2D)
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

    def sync_checkboxes_2D(self, state):
        """Synchronise toutes les cases de normalisation 2D."""
        for checkbox in self.normalize_buttons_2D:
            checkbox.blockSignals(True)  # Évite une boucle infinie
            checkbox.setChecked(state)
            checkbox.blockSignals(False)

        self.update_all_graphs()

    def open_resolution_dialog(self):
        """
        Ouvre une boîte de dialogue pour changer la résolution de la simulation.
        """
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

        # Label for algorithm explanation
        explanation_label = QLabel("")
        explanation_label.setWordWrap(True)
        layout.addWidget(explanation_label)

        # Button to close the dialog
        ok_button = QPushButton("OK")
        ok_button.setEnabled(False)  # Initially disable the button
        layout.addWidget(ok_button)

        # Enable the OK button only when a radio button is selected
        # Function to update the explanation box with synthesized texts
        def update_explanation():
            if linear_button.isChecked():
                # Synthesized explanation for trilinear interpolation (linear)
                explanation_text = (
                    "Trilinear interpolation computes intermediate values by using the weighted average "
                    "of 8 adjacent points within a cube, assuming linear changes along each axis. "
                    "It efficiently increases grid resolution with simple linear approximations."
                )
                ok_button.setEnabled(True)
            elif nonlinear_button.isChecked():
                # Synthesized explanation for non-linear interpolation
                explanation_text = (
                    "Non-linear interpolation incorporates complex physical laws and corrections to capture "
                    "rapid variations in the field, resulting in a more precise estimation. "
                )
                ok_button.setEnabled(True)
            else:
                explanation_text = ""
                ok_button.setEnabled(False)
            explanation_label.setText(explanation_text)
            explanation_label.adjustSize()
            dialog.adjustSize()

        # Connect toggled signals to update the explanation
        linear_button.toggled.connect(update_explanation)
        nonlinear_button.toggled.connect(update_explanation)

        def on_ok_clicked():
            algorithm = "non-linear" if nonlinear_button.isChecked() else "linear"
            self.set_resolution(resolution_input.text(), dialog, algorithm)
            dialog.accept()

        ok_button.clicked.connect(on_ok_clicked)

        dialog.setLayout(layout)
        dialog.exec_()

    def saveToFile(self):
        """Save magnetic field data to CSV file."""
        # Open file dialog for saving
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save File",
            "",
            "CSV Files (*.csv);;All Files (*)",
            options=options
        )

        r = 0.01
        f = 13.56e6  # Fréquence en Hz
        omega = 2 * math.pi * f  # Pulsation angulaire en rad/s
        s = math.pi * r ** 2
        mu_0 = 4 * np.pi * 1e-7  # Perméabilité du vide (T·m/A)

        if file_path:
            # Add .csv extension if not present
            if not file_path.endswith('.csv'):
                file_path += '.csv'

            try:
                with open(file_path, 'w') as file:
                    # Write header
                    file.write('x,y,z,CH1_Max_Voltage,CH2_Max_Voltage,CH3_Max_Voltage\n')

                    # Write data points
                    for point in self.simulation.points_haute_resolution:
                        # Write line with dummy values for voltage channels
                        print(point)
                        ax = point['Hx'] * mu_0 * s * omega
                        ay = point['Hy'] * mu_0 * s * omega
                        az = point['Hz'] * mu_0 * s * omega
                        file.write(f"{point['x']},{point['y']},{point['z']},{ax},{ay},{az}\n")
            except Exception as e:
                print(f"Error saving file: {e}")

    def create_plane_and_value_selectors(self, plane_callback, value_callback):
        """Crée les sélecteurs pour le plan et la valeur."""
        selector_layout = QHBoxLayout()
        selector_layout.setContentsMargins(1, 1, 1, 1)

        # Sélecteur de plan
        plane_label = QLabel("Plan:")
        plane_selector = QComboBox()
        plane_selector.addItems(['YZ', 'XZ', 'XY'])
        plane_selector.currentTextChanged.connect(plane_callback)
        plane_selector.setMinimumWidth(75)  # Réglez la largeur minimale si besoin


        # Sélecteur de valeur
        value_label = QLabel("Coordonnée X:")
        value_selector = QComboBox()
        value_selector.currentTextChanged.connect(value_callback)
        value_selector.setMinimumWidth(75)  # Ajustez la largeur minimale si nécessaire

        # Update the value_label text based on selected plane
        def update_value_label(plane_text):
            if plane_text == 'XY':
                value_label.setText("Coordonnée Z (mm) :")
            elif plane_text == 'XZ':
                value_label.setText("Coordonnée Y (mm):")
            elif plane_text == 'YZ':
                value_label.setText("Coordonnée X (mm):")

        plane_selector.currentTextChanged.connect(update_value_label)

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
        for p in self.simulation.measuredPoints:
            if p not in self.simulation.points_haute_resolution:
                self.simulation.points_haute_resolution.append(p)

        if plane == 'XY':
            axe = 'z'
        if plane == 'XZ':
            axe = 'y'
        if plane == 'YZ':
            axe = 'x'
        if axe in ['x', 'y', 'z']:
            # Extraire les valeurs uniques arrondies
            raw_values = [p[axe] for p in self.simulation.points_haute_resolution]
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
        if plane == "XY":
            axe = "z"
        if plane == "XZ":
            axe = "y"
        if plane == "YZ":
            axe = "x"
        print(axe)
        if self.value_selector_2d.count() == 0:
            return

        try:
            value = float(self.value_selector_2d.currentText())
        except ValueError:
            return

        filtered_points = self.filter_points(axe, value)
        if not filtered_points or len(filtered_points) < 3:  # Ensure enough points for a plane
            print(f"Aucun point trouvé pour le plan {axe}={value}, ou pas assez de points.")
            return

        coord1, coord2, H_total, H_component1_norm, H_component2_norm, H_component1_base, H_component2_base = self.prepare_plot_data(filtered_points, axe)

        # If there are not enough points to interpolate smoothly, do a scatter plot
        if len(np.unique(coord1)) < 2 or len(np.unique(coord2)) < 2:
            print("Insufficient unique points for grid interpolation; using scatter plot.")
            self.figure_2d.clear()
            ax = self.figure_2d.add_subplot(111)
            scatter = ax.scatter(coord1, coord2, c=H_total, cmap='viridis', edgecolor='k')
            ax.quiver(coord1, coord2, H_component1_norm, H_component2_norm, color='red', scale=self.vector_scale_2d)
            ax.set_title(f"Points on plane {axe}")
            if axe == "x":
                ax.set_xlabel(f'Y (mm)')
                ax.set_ylabel(f'Z (mm)')
            if axe == "y":
                ax.set_xlabel(f'X (mm)')
                ax.set_ylabel(f'Z (mm)')
            if axe == "z":
                ax.set_xlabel(f'X (mm)')
                ax.set_ylabel(f'Y (mm)')
            ax.set_ylabel('Other axis (mm)')  # Change accordingly
            self.figure_2d.colorbar(scatter, ax=ax, label='|H| (A/m)')
            self.canvas_2d.draw()
            return

        # Proceed as usual if data is sufficient
        coord1_grid, coord2_grid = np.meshgrid(np.unique(coord1), np.unique(coord2))
        H_total_grid = griddata((coord1, coord2), H_total,
                                (coord1_grid, coord2_grid), method="linear", fill_value=0)
        H_component1_base_grid = griddata((coord1, coord2), H_component1_base,
                                          (coord1_grid, coord2_grid), method="linear", fill_value=0)
        H_component2_base_grid = griddata((coord1, coord2), H_component2_base,
                                          (coord1_grid, coord2_grid), method="linear", fill_value=0)
        H_component1_norm_grid = griddata((coord1, coord2), H_component1_norm,
                                          (coord1_grid, coord2_grid), method="linear", fill_value=0)
        H_component2_norm_grid = griddata((coord1, coord2), H_component2_norm,
                                          (coord1_grid, coord2_grid), method="linear", fill_value=0)

        self.figure_2d.clear()
        gs = self.figure_2d.add_gridspec(1, 3, width_ratios=[6, 0.4, 6])
        ax1 = self.figure_2d.add_subplot(gs[0, 0])
        contour = ax1.contourf(coord1_grid, coord2_grid, H_total_grid, levels=20, cmap='viridis')
        ax1.set_title(f"Norme du champ |H| ({axe})")
        if axe == "x":
            ax1.set_xlabel(f'Y (mm)')
            ax1.set_ylabel(f'Z (mm)')
        if axe == "y":
            ax1.set_xlabel(f'X (mm)')
            ax1.set_ylabel(f'Z (mm)')
        if axe == "z":
            ax1.set_xlabel(f'X (mm)')
            ax1.set_ylabel(f'Y (mm)')
        cbar_ax = self.figure_2d.add_subplot(gs[0, 1])
        self.figure_2d.colorbar(contour, cax=cbar_ax, label="|H| (A/m)")
        self.ax2 = self.figure_2d.add_subplot(gs[0, 2])

        # Draw only the quiver with vectors and no points.
        self.ax2.quiver(coord1_grid, coord2_grid, H_component1_norm_grid, H_component2_norm_grid,
                          color="red", scale=self.vector_scale_2d)
        self.ax2.set_title(f"Magnetic Field Direction on plane {axe}")
        if axe == "x":
            self.ax2.set_xlabel(f'Y (mm)')
            self.ax2.set_ylabel(f'Z (mm)')
        if axe == "y":
            self.ax2.set_xlabel(f'X (mm)')
            self.ax2.set_ylabel(f'Z (mm)')
        if axe == "z":
            self.ax2.set_xlabel(f'X (mm)')
            self.ax2.set_ylabel(f'Y (mm)')

        self.ax2.set_aspect("equal")
        pos = self.ax2.get_position()
        # Décaler ax2 vers la droite permet d'augmenter visuellement l'espace entre la légende et ce graphique,
        # tout en gardant un écart réduit entre le premier graphique et la légende.
        new_pos = [pos.x0 + 0.05, pos.y0, pos.width, pos.height]  # Valeur ajustable selon le besoin
        self.ax2.set_position(new_pos)

        # Store vector base positions and components for interactivity.

        self.vector_bases = np.array([coord1_grid.flatten(), coord2_grid.flatten()]).T
        self.vector_hcb1 = H_component1_base_grid.flatten()
        self.vector_hcb2 = H_component2_base_grid.flatten()
        self.vector_hcn1 = H_component1_norm_grid.flatten()
        self.vector_hcn2 = H_component2_norm_grid.flatten()
        self.vector_ht = H_total_grid.flatten()

        # Create the annotation once.
        self.annotation = self.ax2.annotate(
            text="",
            xy=(0, 0),
            xytext=(6, 15),
            textcoords="offset points",
            bbox={"boxstyle": "round", "fc": "w"},
            arrowprops={"arrowstyle": "->"}
        )
        self.annotation.set_visible(False)

        # Initialize interactivity for ax2
        self.init_interactivity(self.ax2, self.canvas_2d)


    def update_tab_gaussian_and_radial(self):
        """Met à jour les graphiques de l'onglet Champ Amplitude et Vectoriel."""
        plane = self.plane_selector_3d.currentText()
        if plane == "XY":
            axe = "z"
        if plane == "XZ":
            axe = "y"
        if plane == "YZ":
            axe = "x"
        if self.value_selector_3d.count() == 0:
            return

        try:
            value = float(self.value_selector_3d.currentText())
        except ValueError:
            return

        filtered_points = self.filter_points(axe, value)
        if not filtered_points or len(filtered_points) < 3:
            print(f"Aucun point trouvé pour le plan {axe}={value}, ou pas assez de points.")
            return

        coord1, coord2, H_total, H_component1_norm, H_component2_norm, H_component1_base, H_component2_base = self.prepare_plot_data(filtered_points, axe)

        # Check if there are enough unique points for interpolation
        if len(np.unique(coord1)) < 2 or len(np.unique(coord2)) < 2:
            print("Insufficient unique points for grid interpolation; using scatter plot.")
            self.figure_gaussian.clear()
            ax1 = self.figure_gaussian.add_subplot(111, projection='3d')
            ax1.scatter(coord1, coord2, H_total, c=H_total, cmap='viridis', edgecolor='k', alpha=0.8)
            ax1.set_title(f"Points sur le plan {axe}")
            if axe == "x":
                ax1.set_xlabel(f'Y (mm)')
                ax1.set_ylabel(f'Z (mm)')
            if axe == "y":
                ax1.set_xlabel(f'X (mm)')
                ax1.set_ylabel(f'Z (mm)')
            if axe == "z":
                ax1.set_xlabel(f'X (mm)')
                ax1.set_ylabel(f'Y (mm)')
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
        H_component1_base_grid = griddata((coord1, coord2), H_component1_base, (coord1_grid, coord2_grid),
                                          method='linear', fill_value=0)
        H_component2_base_grid = griddata((coord1, coord2), H_component2_base, (coord1_grid, coord2_grid),
                                          method='linear', fill_value=0)

        self.figure_gaussian.clear()
        gs = self.figure_gaussian.add_gridspec(1, 2, width_ratios=[1, 1])

        ax1 = self.figure_gaussian.add_subplot(gs[0, 0], projection='3d')
        surf = ax1.plot_surface(coord1_grid, coord2_grid, H_total_grid, cmap='viridis', edgecolor='k', alpha=0.8)
        ax1.set_title(f'Amplitude du champ magnétique |H| ({axe})')
        if axe == "x":
            ax1.set_xlabel(f'Y (mm)')
            ax1.set_ylabel(f'Z (mm)')
        if axe == "y":
            ax1.set_xlabel(f'X (mm)')
            ax1.set_ylabel(f'Z (mm)')
        if axe == "z":
            ax1.set_xlabel(f'X (mm)')
            ax1.set_ylabel(f'Y (mm)')
        ax1.set_zlabel('Amplitude |H| (A/m)')
        self.figure_gaussian.colorbar(surf, ax=ax1, shrink=0.5, aspect=10)

        self.ax_gaussian = self.figure_gaussian.add_subplot(gs[0, 1])
        self.ax_gaussian.quiver(coord1_grid, coord2_grid, H_component1_norm_grid, H_component2_norm_grid, scale=self.vector_scale_2d)
        self.ax_gaussian.set_title(f"Champ vectoriel sur le plan {axe}")
        if axe == "x":
            self.ax_gaussian.set_xlabel(f'Y (mm)')
            self.ax_gaussian.set_ylabel(f'Z (mm)')
        if axe == "y":
            self.ax_gaussian.set_xlabel(f'X (mm)')
            self.ax_gaussian.set_ylabel(f'Z (mm)')
        if axe == "z":
            self.ax_gaussian.set_xlabel(f'X (mm)')
            self.ax_gaussian.set_ylabel(f'Y (mm)')
        self.ax_gaussian.set_aspect('equal')
        # Store vector base positions and components for interactivity.

        self.gaussian_vector_bases = np.array([coord1_grid.flatten(), coord2_grid.flatten()]).T
        self.gaussian_vector_hcb1 = H_component1_base_grid.flatten()
        self.gaussian_vector_hcb2 = H_component2_base_grid.flatten()
        self.gaussian_vector_hcn1 = H_component1_norm_grid.flatten()
        self.gaussian_vector_hcn2 = H_component2_norm_grid.flatten()
        self.gaussian_vector_ht = H_total_grid.flatten()

        # Create the annotation once.
        self.gaussian_annotation = self.ax_gaussian.annotate(
            text="",
            xy=(0, 0),
            xytext=(6, 15),
            textcoords="offset points",
            bbox={"boxstyle": "round", "fc": "w"},
            arrowprops={"arrowstyle": "->"}
        )
        self.gaussian_annotation.set_visible(False)
        # Initialize interactivity for ax_gaussian
        self.init_interactivity(self.ax_gaussian, self.canvas_gaussian)

    def init_interactivity(self, ax, canvas):
        """Initialize interactivity for a given axis and canvas."""
        ax.set_autoscale_on(True)
        self.dragging = False
        self.previous_point = None

        canvas.mpl_connect("motion_notify_event", self.motion_hover)
        canvas.mpl_connect('scroll_event', self.on_scroll)
        canvas.mpl_connect('button_press_event', self.on_press)
        canvas.mpl_connect('button_release_event', self.on_release)
        canvas.mpl_connect('motion_notify_event', self.on_motion)

        canvas.draw()

    def on_scroll(self, event):
        """Handle scroll events for zooming."""
        if event.inaxes is None:
            return

        ax = event.inaxes
        canvas = ax.figure.canvas

        # Get the current x and y limits
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()

        # Get the current mouse position
        x, y = event.xdata, event.ydata

        # Zoom factor
        scale_factor = 1.1 if event.button == 'up' else 0.9

        # Set new limits
        new_width = (xlim[1] - xlim[0]) * scale_factor
        new_height = (ylim[1] - ylim[0]) * scale_factor

        ax.set_xlim([x - new_width * (x - xlim[0]) / (xlim[1] - xlim[0]),
                     x + new_width * (xlim[1] - x) / (xlim[1] - xlim[0])])
        ax.set_ylim([y - new_height * (y - ylim[0]) / (ylim[1] - ylim[0]),
                     y + new_height * (ylim[1] - y) / (ylim[1] - ylim[0])])

        canvas.draw_idle()

    def on_press(self, event):
        """Handle mouse button press for panning."""
        if event.inaxes is None:
            return

        if event.button == 1:  # Left mouse button
            self.dragging = True
            self.previous_point = (event.xdata, event.ydata)

    def on_release(self, event):
        """Handle mouse button release."""
        self.dragging = False
        self.previous_point = None

    def on_motion(self, event):
        """Handle mouse motion for panning and hover annotations."""
        ax = event.inaxes
        if ax is None:
            return

        canvas = ax.figure.canvas

        # First check if we're doing a panning operation
        if self.dragging and event.xdata and event.ydata:
            # Calculate the movement
            dx = self.previous_point[0] - event.xdata
            dy = self.previous_point[1] - event.ydata

            # Update the view limits
            ax.set_xlim(ax.get_xlim() + dx)
            ax.set_ylim(ax.get_ylim() + dy)

            # Update the previous point
            self.previous_point = (event.xdata, event.ydata)

            # Redraw the canvas
            canvas.draw_idle()
        # If not panning, handle hover functionality
        elif not self.dragging:
            # Your existing hover functionality
            self.motion_hover(event)

    def filter_points(self, plane, value, epsilon=1e-5):
        """Filtre les points selon le plan et la valeur donnée."""
        return [p for p in self.simulation.points_haute_resolution if abs(p[plane] - value) < epsilon]

    def prepare_plot_data(self, filtered_points, plane):
        """Prépare les données de tracé pour les graphiques 2D en fonction de la normalisation."""
        axes_vars = {'x': ['y', 'z'], 'y': ['x', 'z'], 'z': ['x', 'y']}
        axis1, axis2 = axes_vars[plane]

        coord1 = np.array([p[axis1] for p in filtered_points])
        coord2 = np.array([p[axis2] for p in filtered_points])

        H_total = np.sqrt(
            np.array([p['Hx'] for p in filtered_points]) ** 2 +
            np.array([p['Hy'] for p in filtered_points]) ** 2 +
            np.array([p['Hz'] for p in filtered_points]) ** 2
        )
        H_component1_base = np.array([p[f'H{axis1}'] for p in filtered_points])
        H_component2_base = np.array([p[f'H{axis2}'] for p in filtered_points])


        if self.normalize_button_2D.isChecked():
            with np.errstate(divide='ignore', invalid='ignore'):
                H_component1 = np.where(H_total != 0, H_component1_base / H_total, 0)
                H_component2 = np.where(H_total != 0, H_component2_base / H_total, 0)
        else:
            H_component1 = H_component1_base
            H_component2 = H_component2_base

        return coord1, coord2, H_total, H_component1, H_component2, H_component1_base, H_component2_base

    def motion_hover(self, event):
        if self.dragging:
            return
        if event.inaxes == self.ax2 and hasattr(self, "vector_bases"):
            self._handle_hover(event, "2d")
        elif event.inaxes == self.ax_gaussian and hasattr(self, "gaussian_vector_bases"):
            self._handle_hover(event, "3d")

    def _handle_hover(self, event, mode):
        """Handle hover logic for both 2D and 3D modes."""
        config = self._get_config(mode)
        P = np.array([event.xdata, event.ydata])
        best_index = self._find_closest_vector(P, config)

        if best_index is not None:
            self._update_annotation(best_index, config)
        else:
            self._hide_annotation(config)

    def _get_config(self, mode):
        """Get configuration based on mode."""
        if mode == "2d":
            return {
                'vectors': self.vector_bases,
                'hcn1': self.vector_hcn1,
                'hcn2': self.vector_hcn2,
                'hcb1': self.vector_hcb1,
                'hcb2': self.vector_hcb2,
                'ht': 'vector_ht',
                'annotation': self.annotation,
                'canvas': self.canvas_2d,
                'plane_selector': self.plane_selector_2d,
                'use_projection': True
            }
        return {
            'vectors': self.gaussian_vector_bases,
            'hcn1': self.gaussian_vector_hcn1,
            'hcn2': self.gaussian_vector_hcn2,
            'hcb1': self.gaussian_vector_hcb1,
            'hcb2': self.gaussian_vector_hcb2,
            'ht': 'gaussian_vector_ht',
            'annotation': self.gaussian_annotation,
            'canvas': self.canvas_gaussian,
            'plane_selector': self.plane_selector_3d,
            'use_projection': False
        }

    def _find_closest_vector(self, P, config):
        """Find the closest vector to point P."""
        tol = 1
        best_index = None
        best_distance = tol

        for i, base in enumerate(config['vectors']):
            A = np.array(base)
            hcn1, hcn2 = config['hcn1'][i], config['hcn2'][i]

            if config['use_projection']:
                distance = self._calculate_projection_distance(P, A, hcn1, hcn2)
            else:
                vector_tip = A + np.array([hcn1, hcn2])
                distance = np.linalg.norm(P - vector_tip)

            if distance < best_distance:
                best_distance = distance
                best_index = i
                if not config['use_projection']:
                    break

        return best_index

    def _calculate_projection_distance(self, P, A, hcn1, hcn2):
        """Calculate distance using vector projection."""
        B = A + np.array([hcn1, hcn2])
        AB = B - A
        AP = P - A

        if np.dot(AB, AB) == 0:
            return np.linalg.norm(AP)

        t = np.dot(AP, AB) / np.dot(AB, AB)
        if t < 0:
            closest = A
        elif t > 1:
            closest = B
        else:
            closest = A + t * AB
        return np.linalg.norm(P - closest)

    def _update_annotation(self, index, config):
        """Update annotation with vector information."""
        A = np.array(config['vectors'][index])
        hcn1, hcn2 = config['hcn1'][index], config['hcn2'][index]
        hcb1, hcb2 = config['hcb1'][index], config['hcb2'][index]
        tip = A + np.array([hcn1, hcn2])

        plane = config['plane_selector'].currentText()
        label_text = self._create_label_text(plane, hcb1, hcb2, config['ht'], index)

        annotation = config['annotation']
        annotation.xy = A
        annotation.set_text(label_text)
        annotation.get_bbox_patch().set_facecolor("yellow")
        annotation.set_alpha(0.8)
        annotation.set_visible(True)
        config['canvas'].draw_idle()

    def _create_label_text(self, plane, hcb1, hcb2, ht_attr, index):
        """Create label text for annotation."""
        base_text = f"H{plane[0]}={hcb1:.2f}\nH{plane[1]}={hcb2:.2f}"
        if hasattr(self, ht_attr):
            h = getattr(self, ht_attr)[index]
            return f"{base_text}\n|H|={h:.2f}"
        return base_text

    def _hide_annotation(self, config):
        """Hide annotation if visible."""
        if config['annotation'].get_visible():
            config['annotation'].set_visible(False)
            config['canvas'].draw_idle()