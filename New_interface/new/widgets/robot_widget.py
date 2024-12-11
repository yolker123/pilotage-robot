from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QPushButton, QSlider, QLabel, QComboBox, QHBoxLayout, QSpacerItem
)
from PyQt5.QtCore import Qt

class RobotWidget(QWidget):
    button_clicked = pyqtSignal(str)

    def __init__(self, visualisation_measure_widget):
        super().__init__()
        self.visualisation_measure_widget = visualisation_measure_widget
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)

        title_label = QLabel("Robot")
        layout.addWidget(title_label)

        layout.addWidget(self.create_port_dropdown())
        layout.addWidget(self.create_robot_dropdown())
        layout.addWidget(self.create_connection_button())
        layout.addWidget(self.create_home_position_button())

        layout.addWidget(QLabel("Vitesse"))
        layout.addLayout(self.create_speed_layout())

        layout.addWidget(QLabel("Contrôle du mouvement"))
        layout.addLayout(self.create_measure_method_choice_layout())
        layout.addLayout(self.create_step_layout())

        layout.addWidget(self.create_launch_button())

        layout.addWidget(self.create_emergency_button(), alignment=Qt.AlignCenter)

        self.setLayout(layout)

    # --- Widgets -----------------
    def create_port_dropdown(self) -> QComboBox:
        dropdown = QComboBox()
        dropdown.addItems(["COM1", "COM2", "COM3", "COM4"])
        # dropdown.currentIndexChanged.connect(on_dropdown_change)
        return dropdown

    def create_robot_dropdown(self) -> QComboBox:
        dropdown = QComboBox()
        dropdown.addItems(["Robot 5 axis", "Robot 6 axis"])
        # dropdown.currentIndexChanged.connect(on_dropdown_change)
        return dropdown

    def create_connection_button(self) -> QPushButton:
        button = QPushButton("Connection au robot")
        # button.click.connect(on_dropdown_change)
        return button

    def create_home_position_button(self) -> QPushButton:
        button = QPushButton("Position home")
        # button.click.connect(on_dropdown_change)
        return button

    def create_speed_layout(self) -> QVBoxLayout:
        layout = QVBoxLayout()

        self.slider_label = QLabel("Vitesse: 50")
        self.slider_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.slider_label)

        # Créer le slider
        slider = QSlider(Qt.Horizontal)
        slider.setRange(0, 100)
        slider.setValue(50)
        slider.valueChanged.connect(self.on_slider_change)

        layout.addWidget(slider)
        return layout

    def create_measure_method_choice_layout(self) -> QVBoxLayout:
        layout = QVBoxLayout()
        method_layout = QHBoxLayout()

        method_layout.addWidget(QLabel("Méthode de mesure : "))
        self.measure_method_dropdown = QComboBox()
        self.measure_method_dropdown.addItems(['Manual Measure', 'NFC', 'EMVCO', 'Custom Cube', 'One point', 'Custom Cylinder','Semi-sphere'])
        method_layout.addWidget(self.measure_method_dropdown)
        self.measure_method_dropdown.currentIndexChanged.connect(self.on_measure_method_change)

        layout.addLayout(method_layout)
        return layout

    def create_step_layout(self) -> QVBoxLayout:
        layout = QVBoxLayout()

        step_layout = QHBoxLayout()
        step_layout.addWidget(QLabel("Pas:"))
        self.step_combo = QComboBox()
        self.step_combo.addItems(["10 mm", "5 mm", "4 mm", "3 mm", "2 mm", "1 mm","0,5 mm"])
        step_layout.addWidget(self.step_combo)
        layout.addLayout(step_layout)

        line_layout = QHBoxLayout()
        # Contrôles X, Y, Z
        control_layout = QVBoxLayout()
        for axis in ['X', 'Y', 'Z']:
            axis_layout = QHBoxLayout()
            axis_layout.setAlignment(Qt.AlignCenter)
            axis_layout.addWidget(QLabel(f"{axis}:"))

            minus_btn = QPushButton("-")
            minus_btn.setFixedSize(40, 40)
            axis_layout.addWidget(minus_btn)

            coord_label = QLabel("0.0 mm")
            # coord_label.setAlignment(Qt.AlignCenter)
            axis_layout.addWidget(coord_label)

            plus_btn = QPushButton("+")
            plus_btn.setFixedSize(40, 40)
            axis_layout.addWidget(plus_btn)
            control_layout.addLayout(axis_layout)
        line_layout.addLayout(control_layout)
        line_layout.addStretch()

        # File section
        button_layout = QVBoxLayout()
        open_button = QPushButton("Open")
        open_button.setFixedSize(300, 40)
        save_button = QPushButton("Save")
        save_button.setFixedSize(300, 40)
        show_button = QPushButton("Show(x, y, z)")
        show_button.setFixedSize(300, 40)
        button_layout.addWidget(open_button)
        button_layout.addWidget(save_button)
        button_layout.addWidget(show_button)
        line_layout.addLayout(button_layout)

        layout.addLayout(line_layout)

        return layout

    def create_launch_button(self) -> QPushButton:
        button = QPushButton("Launch")

        return button

    def create_emergency_button(self) -> QPushButton:
        button = QPushButton("Emergency STOP")
        button.setFixedSize(500, 200)
        return button

    # --- Actions -----------------
    def on_slider_change(self, value):
        # Mettre à jour le texte de la vitesse en fonction de la valeur du slider
        vitesse = ""
        if value > 66:
            vitesse = "Rapide"
        elif value < 33:
            vitesse = "Lent"
        elif  value >= 33 and value <= 66:
            vitesse = "Moyen"
        else:
            "Erreur de valeur"
        self.slider_label.setText(f"{vitesse} : {value}")
        print(f"Slider changé à : {value}")

    # ------------------

    ### 1. Fonction pour créer le bouton ###
    def create_button(self):
        button_layout = QVBoxLayout()
        self.button = QPushButton("Lancer l'Action")
        self.button.clicked.connect(self.on_button_click)
        button_layout.addWidget(self.button)
        return button_layout


    ### 2. Fonction pour créer le slider ###
    ### Fonction pour créer le slider avec la mise à jour de la vitesse ###
    def create_slider(self):
        layout = QVBoxLayout()  # Layout principal pour le slider

        # Ajouter un label pour afficher la vitesse sélectionnée (initialisée à 50)
        self.slider_label = QLabel("Vitesse: 50")
        self.slider_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.slider_label)

        # Créer le slider
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setRange(0, 100)  # Définir la plage du slider de 0 à 100
        self.slider.setValue(50)  # Valeur initiale du slider
        self.slider.valueChanged.connect(self.on_slider_change)  # Connecter le signal du slider à la méthode

        # Diviser le slider en 3 sections (Lent, Moyen, Rapide) visuellement
        self.slider.setTickInterval(33)  # Diviser en 3 sections égales : [0-33], [34-66], [67-100]
        self.slider.setTickPosition(QSlider.TicksBelow)  # Positionner les ticks en dessous du slider

        layout.addWidget(self.slider)  # Ajouter le slider au layout

        return layout


    ### 3. Fonction pour créer le dropdown (ComboBox) ###
    def create_dropdown(self):
        dropdown_layout = QVBoxLayout()
        self.dropdown = QComboBox()
        self.dropdown.addItems(["Option 1", "Option 2", "Option 3"])
        self.dropdown.currentIndexChanged.connect(self.on_dropdown_change)

        dropdown_layout.addWidget(QLabel("Choisissez une option :"))
        dropdown_layout.addWidget(self.dropdown)
        return dropdown_layout


    ### Gestion des événements ###
    def on_button_click(self):
        print(f"Bouton cliqué! Valeur du slider : {self.slider.value()}")
        print(f"Option sélectionnée : {self.dropdown.currentText()}")


    ### Gestion de l'événement de changement de valeur du slider ###


    def on_measure_method_change(self, index):
        print(f"Dropdown changé à l'index {index} : {self.measure_method_dropdown.currentText()}")
        self.visualisation_measure_widget.change_image(self.measure_method_dropdown.currentText())
