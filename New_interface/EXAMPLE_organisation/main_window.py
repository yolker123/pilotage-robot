from PyQt5.QtWidgets import QMainWindow, QTabWidget
from tabs.tab_calculator import CalculatorTab
from tabs.tab_tasks import TasksTab
from tabs.tab_counter import CounterTab

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Application Multi-Logique")
        self.setGeometry(100, 100, 800, 600)
        self.init_ui()

    def init_ui(self):
        tabs = QTabWidget()

        # Ajouter des onglets
        tabs.addTab(CalculatorTab(), "Calculatrice")
        tabs.addTab(TasksTab(), "Gestion des Tâches")
        tabs.addTab(CounterTab(), "Compteur de Clics")

        self.setCentralWidget(tabs)
