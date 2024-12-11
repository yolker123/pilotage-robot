import os

from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel
from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt

class VisualisationMeasuresWidget(QWidget):
    def __init__(self):
        super().__init__()
        self.simulation = None
        self.image_label = None  # Label pour afficher l'image

        test = os.system('pwd')
        print(test)

        self.images = {
            'Manual Measure': 'assets/measure.png',
            'NFC': 'assets/nfc.jpg',
            'EMVCO': 'assets/emvco.jpg',
            'Custom Cube': 'assets/cube.png',
            'One point': 'assets/point.png',
            'Custom Cylinder': 'assets/point.png',
            'Semi-sphere': 'assets/point.png'
        }

        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)

        # Initialisation de la simulation (si nécessaire)
        # self.simulation = MagneticFieldSimulation_brouillon(self, width=8, height=6, dpi=100)

        # Créer le label pour afficher l'image
        self.image_label = QLabel(self)
        self.pixmap = QPixmap(self.images['Manual Measure'])
        self.image_label.setPixmap(self.pixmap)
        self.image_label.setAlignment(Qt.AlignCenter)  # Centrer l'image
        layout.addWidget(self.image_label)

    def change_image(self, image):
        if image in self.images:
            self.pixmap = QPixmap(self.images[image])
            self.image_label.setPixmap(self.pixmap)
        else:
            print(f"Image '{image}' non trouvée dans le dictionnaire des images.")
