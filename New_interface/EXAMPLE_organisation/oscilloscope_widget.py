from PyQt5.QtWidgets import QWidget, QVBoxLayout, QPushButton, QGridLayout, QLabel
from PyQt5.QtCore import pyqtSignal

class OscilloscopeWidget(QWidget):
    button_clicked = pyqtSignal(str)

    def __init__(self):
        super().__init__()
        self.button_actions = {
                "Connection à l'oscilloscope": self.action_oscilloscope_connection,
            }
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        grid_layout = QGridLayout()
        self.setLayout(layout)

        title_label = QLabel("Oscilloscope")
        layout.addWidget(title_label)

        for index, (label, action) in enumerate(self.button_actions.items()):
            button = QPushButton(label)
            button.clicked.connect(lambda checked, lbl=label: self.handle_button_click(lbl))
            row = index // 3
            col = index % 3
            grid_layout.addWidget(button, row, col)

        layout.addLayout(grid_layout)

    def handle_button_click(self, label):
        self.button_clicked.emit(label)

        # Exécute l'action associée
        if label in self.button_actions:
            self.button_actions[label]()

    def action_oscilloscope_connection(self):
        print("oscilloscope connection...")
