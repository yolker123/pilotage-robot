from PyQt5.QtWidgets import QWidget, QVBoxLayout, QPushButton, QLabel
from logic.counter_logic import CounterLogic

class CounterTab(QWidget):
    def __init__(self):
        super().__init__()
        self.counter_logic = CounterLogic()
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        self.counter_label = QLabel(f"Compteur : {self.counter_logic.get_count()}")
        increment_button = QPushButton("Incrémenter")

        layout.addWidget(self.counter_label)
        layout.addWidget(increment_button)

        self.setLayout(layout)

        # Connexion du bouton à la logique
        increment_button.clicked.connect(self.increment_counter)

    def increment_counter(self):
        self.counter_logic.increment()
        self.counter_label.setText(f"Compteur : {self.counter_logic.get_count()}")
