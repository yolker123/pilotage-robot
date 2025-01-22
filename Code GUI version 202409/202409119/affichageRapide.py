import sys
import numpy as np
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import (QApplication, QMainWindow, QVBoxLayout, QWidget,
                             QSpinBox, QLabel, QProgressBar)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from mpl_toolkits.mplot3d import Axes3D
from PyQt5.QtCore import QThread, pyqtSignal
import psutil
import time
from multiprocessing import Pool, cpu_count
from functools import partial


def generate_points_chunk(chunk_size, cpu_limit=80, ram_limit=80):
    """Génère un lot de points de manière optimisée"""
    while psutil.cpu_percent() > cpu_limit or psutil.virtual_memory().percent > ram_limit:
        time.sleep(0.05)

    # Générer tous les points d'un coup avec numpy
    points = np.random.rand(chunk_size, 3)
    return points


class PointsGenerator(QThread):
    points_ready = pyqtSignal(object)
    progress_updated = pyqtSignal(int)

    def __init__(self, n_points, cpu_limit, ram_limit):
        super().__init__()
        self.n_points = n_points
        self.cpu_limit = cpu_limit
        self.ram_limit = ram_limit
        self.running = True
        # Utiliser 75% des cœurs disponibles pour laisser de la marge pour l'interface
        self.n_processes = max(1, int(cpu_count() * 0.75))

    def run(self):
        # Calculer la taille optimale des chunks basée sur le nombre de processus
        chunk_size = max(2000, self.n_points // (self.n_processes * 4))
        n_chunks = (self.n_points + chunk_size - 1) // chunk_size

        # Créer un pool de processus
        with Pool(processes=self.n_processes) as pool:
            # Préparer la fonction partielle avec les limites
            gen_func = partial(generate_points_chunk,
                               cpu_limit=self.cpu_limit,
                               ram_limit=self.ram_limit)

            # Générer les chunks de points en parallèle
            points_array = np.empty((0, 3))
            for i, chunk_points in enumerate(pool.imap(gen_func, [chunk_size] * n_chunks)):
                if not self.running:
                    pool.terminate()
                    break

                # Ajouter les nouveaux points
                points_array = np.vstack((points_array, chunk_points))

                # Limiter au nombre de points demandé
                if len(points_array) > self.n_points:
                    points_array = points_array[:self.n_points]

                # Mettre à jour la progression
                progress = min((len(points_array) / self.n_points) * 100, 100)
                self.progress_updated.emit(int(progress))

                # Émettre les points pour mise à jour
                # Séparer les coordonnées x, y, z
                self.points_ready.emit((points_array[:, 0],
                                        points_array[:, 1],
                                        points_array[:, 2]))

        self.progress_updated.emit(100)

    def stop(self):
        self.running = False


class Simple3DVisualizer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Visualisation 3D Interactive Optimisée')
        self.setGeometry(100, 100, 800, 600)

        # Configuration de l'interface
        self._setup_ui()

        # Configuration du plot
        self._setup_plot()

        # Variables pour la gestion des ressources
        self.current_thread = None
        self.cpu_limit = 50
        self.ram_limit = 50

        # Démarrer la visualisation
        self.start_update()

    def _setup_ui(self):
        """Configure l'interface utilisateur"""
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout(central_widget)

        # Contrôle du nombre de points
        self.points_control = QSpinBox()
        self.points_control.setMinimum(1)
        self.points_control.setMaximum(2147483647)
        self.points_control.setValue(10000)
        self.points_control.valueChanged.connect(self.start_update)

        # Configuration du layout
        layout.addWidget(QLabel("Nombre de points:"))
        layout.addWidget(self.points_control)

        # Barre de progression
        self.progress_bar = QProgressBar()
        self.progress_bar.hide()
        layout.addWidget(self.progress_bar)

        # Canvas matplotlib
        self.fig = plt.figure()
        self.canvas = FigureCanvas(self.fig)
        layout.addWidget(self.canvas)

    def _setup_plot(self):
        """Configure le plot matplotlib"""
        self.ax = self.fig.add_subplot(111, projection='3d')

        # Optimiser les paramètres du plot
        self.fig.set_tight_layout(True)
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')

    def start_update(self):
        """Démarre la génération des points"""
        if self.current_thread and self.current_thread.isRunning():
            self.current_thread.stop()
            self.current_thread.wait()

        # Réinitialiser l'interface
        self.progress_bar.setValue(0)
        self.progress_bar.show()

        # Créer et démarrer le nouveau thread
        self.current_thread = PointsGenerator(
            self.points_control.value(),
            self.cpu_limit,
            self.ram_limit
        )
        self.current_thread.points_ready.connect(self.update_plot)
        self.current_thread.progress_updated.connect(self.progress_bar.setValue)
        self.current_thread.finished.connect(self.progress_bar.hide)
        self.current_thread.finished.connect(self.draw)
        self.current_thread.start()

    def update_plot(self, points_data):
        """Met à jour le plot avec les nouveaux points"""
        x, y, z = points_data

        # Effacer et tracer
        self.ax.clear()
        self.ax.scatter(x, y, z, alpha=0.6, s=1)

        # Mettre à jour le titre
        self.ax.set_title(f'Visualisation de {len(x)} points en 3D')

        # Rafraîchir le canvas
        # self.canvas.draw()
    def draw(self):
        self.canvas.draw()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = Simple3DVisualizer()
    window.show()
    sys.exit(app.exec_())