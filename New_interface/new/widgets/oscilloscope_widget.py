from PyQt5.QtWidgets import QWidget, QVBoxLayout, QPushButton, QGridLayout, QLabel, QFrame
from PyQt5.QtCore import pyqtSignal

class OscilloscopeWidget(QWidget):

    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        frame = QFrame()
        layout = QVBoxLayout(frame)

        layout.addWidget(QLabel("Oscilloscope"))
        layout.addWidget(self.create_connection_button())


    def create_connection_button(self) -> QPushButton:
        button = QPushButton("Connecter l'oscilloscope")
        button.clicked.connect(self.action_connection())
        return button

    def action_connection(self):
        print("oscilloscope connection...")
