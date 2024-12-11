from PyQt5.QtWidgets import QMainWindow, QTabWidget
from tabs.measures_tab import MeasuresTab

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Robot Bureau d'Étude")
        # self.setGeometry(100, 100, 800, 600)
        self.showMaximized()
        self.init_ui()

    def init_ui(self):
        tabs = QTabWidget()

        # Ajouter des onglets
        tabs.addTab(MeasuresTab(), "Mesures")

        self.setCentralWidget(tabs)
