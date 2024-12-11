from PyQt5.QtWidgets import QWidget, QVBoxLayout, QPushButton, QLineEdit, QLabel
from logic.calculator_logic import CalculatorLogic

class CalculatorTab(QWidget):
    def __init__(self):
        super().__init__()
        self.calculator = CalculatorLogic()
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        self.input_field = QLineEdit()
        self.result_label = QLabel("Résultat : ")
        calc_button = QPushButton("Calculer")

        layout.addWidget(self.input_field)
        layout.addWidget(calc_button)
        layout.addWidget(self.result_label)

        self.setLayout(layout)

        # Connexion du bouton à la logique
        calc_button.clicked.connect(self.calculate)

    def calculate(self):
        expression = self.input_field.text()
        result = self.calculator.evaluate_expression(expression)
        self.result_label.setText(f"Résultat : {result}")
