import sys
import numpy as np
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QTabWidget, QVBoxLayout, QWidget
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from scipy.constants import mu_0
from scipy.integrate import quad
import math


class MagneticFieldSimulation:
    def __init__(self, R=0.05, I=1.0):
        self.R = R
        self.I = I
        self.mu0 = 4 * math.pi * 1e-7
        self.resultats = []
        self.generate_simulated_points()

    def magnetic_field_helix(self, x, y, z, R=0.1, I=1.0):
        """Calcule le champ magnétique d'une spire au point (x, y, z)."""
        def dB_element(phi):
            x0 = R * np.cos(phi)
            y0 = R * np.sin(phi)
            z0 = 0  # Spire dans le plan xy
            r_vec = np.array([x - x0, y - y0, z - z0])
            r_mag = np.linalg.norm(r_vec)
            if r_mag == 0:  # Évite la singularité
                return np.array([0, 0, 0])
            dL = np.array([-R * np.sin(phi), R * np.cos(phi), 0])  # Direction du courant
            dB = (mu_0 * I / (4 * np.pi)) * np.cross(dL, r_vec) / (r_mag ** 3)
            return dB

        Bx, _ = quad(lambda phi: dB_element(phi)[0], 0, 2 * np.pi)
        By, _ = quad(lambda phi: dB_element(phi)[1], 0, 2 * np.pi)
        Bz, _ = quad(lambda phi: dB_element(phi)[2], 0, 2 * np.pi)
        return np.array([Bx, By, Bz])

    def generate_simulated_points(self):
        """Simule le champ magnétique dans une grille 3D."""
        x_min, x_max = -0.1, 0.1
        y_min, y_max = -0.1, 0.1
        z_min, z_max = -0.1, 0.1

        self.resultats = [{'x': x, 'y': y, 'z': z}
                          for x in np.linspace(x_min, x_max, 5)
                          for y in np.linspace(y_min, y_max, 20)
                          for z in np.linspace(z_min, z_max, 20)]

        for point in self.resultats:
            Bx, By, Bz = self.magnetic_field_helix(point['x'], point['y'], point['z'])
            Hx, Hy, Hz = Bx / self.mu0, By / self.mu0, Bz / self.mu0
            H_total = np.sqrt(Hx**2 + Hy**2 + Hz**2)
            point.update({'Hx': Hx, 'Hy': Hy, 'Hz': Hz, 'H_total': H_total})

        return self.resultats


class MagneticFieldApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Simulation du Champ Magnétique")
        self.setGeometry(100, 100, 1200, 600)

        # Initialise la simulation
        self.simulation = MagneticFieldSimulation()

        # Création de la disposition principale
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)

        # Ajouter les onglets
        self.add_tab_3d_vectors()
        self.add_tab_2d_plane()
        self.add_tab_gaussian_and_radial()

    def add_tab_3d_vectors(self):
        """Onglet 1 : Affichage 3D des vecteurs."""
        tab = QWidget()
        layout = QVBoxLayout()
        tab.setLayout(layout)

        figure = plt.figure()
        canvas = FigureCanvas(figure)
        ax = figure.add_subplot(111, projection='3d')

        # Données pour les vecteurs
        points = self.simulation.resultats
        for point in points:
            x, y, z = point['x'], point['y'], point['z']
            Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
            ax.quiver(x, y, z, Hx, Hy, Hz, color='b', length=0.005, normalize=True)

        ax.set_title("Vecteurs 3D du champ magnétique")
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        layout.addWidget(canvas)
        self.tabs.addTab(tab, "Vecteurs 3D")

    def add_tab_2d_plane(self):
        """Onglet 2 : Affichage 2D dans un plan."""
        tab = QWidget()
        layout = QVBoxLayout()
        tab.setLayout(layout)

        figure, axes = plt.subplots(1, 2, figsize=(12, 6))
        canvas = FigureCanvas(figure)

        # Filtrer les points pour le plan x=0
        points = [p for p in self.simulation.resultats if np.isclose(p['x'], 0)]
        Y = np.array([p['y'] for p in points])
        Z = np.array([p['z'] for p in points])
        H_total = np.array([p['H_total'] for p in points])
        Hy = np.array([p['Hy'] for p in points])
        Hz = np.array([p['Hz'] for p in points])

        # Reshape pour une grille 2D
        grid_size = int(np.sqrt(len(Y)))  # Supposer une grille carrée
        Y = Y.reshape((grid_size, grid_size))
        Z = Z.reshape((grid_size, grid_size))
        H_total = H_total.reshape((grid_size, grid_size))
        Hy = Hy.reshape((grid_size, grid_size))
        Hz = Hz.reshape((grid_size, grid_size))

        # Tracer |H| (lignes de niveau)
        contour = axes[0].contourf(Y, Z, H_total, levels=20, cmap='viridis')
        figure.colorbar(contour, ax=axes[0], label='|H| (A/m)')
        axes[0].set_title("Norme du champ |H| (y,z)")
        axes[0].set_xlabel('y (m)')
        axes[0].set_ylabel('z (m)')

        # Tracer le champ vectoriel
        quiver = axes[1].quiver(Y, Z, Hy, Hz, H_total, cmap='inferno', scale=100, scale_units='width', angles='xy')
        axes[1].set_title("Champ vectoriel (Hy, Hz)")
        axes[1].set_xlabel('y (m)')
        axes[1].set_ylabel('z (m)')
        axes[1].set_xlim(-0.1, 0.1)
        axes[1].set_ylim(-0.1, 0.1)
        axes[1].set_aspect('equal')  # Assure que les axes ont la même échelle

        # Ajouter une barre de couleur pour le quiver
        cbar = figure.colorbar(quiver, ax=axes[1], shrink=0.5, aspect=10, label='Magnitude |H| (A/m)')

        layout.addWidget(canvas)
        self.tabs.addTab(tab, "Champ dans le plan 2D")

    def add_tab_gaussian_and_radial(self):
        """Onglet 3 : Affichage du champ radial et de l'amplitude simulée."""
        tab = QWidget()
        layout = QVBoxLayout()
        tab.setLayout(layout)

        figure, axes = plt.subplots(1, 2, figsize=(12, 6))
        canvas = FigureCanvas(figure)

        # Filtrer les points pour le plan x=0
        points = [p for p in self.simulation.resultats if np.isclose(p['x'], 0)]
        H_total = np.array([p['H_total'] for p in points])
        Hy = np.array([p['Hy'] for p in points])
        Hz = np.array([p['Hz'] for p in points])
    
        Y = np.array([p['y'] for p in points])
        Z = np.array([p['z'] for p in points])
        
        # Reshape pour une grille 2D (y,z)
        grid_size = int(np.sqrt(len(points)))  # Supposer une grille carrée
        if grid_size ** 2 != len(points):
            raise ValueError("Le nombre de points ne forme pas une grille carrée.")
        H_total = H_total.reshape((grid_size, grid_size))
        Hy = Hy.reshape((grid_size, grid_size))
        Hz = Hz.reshape((grid_size, grid_size))

        # Créer une grille pour y et z basée sur les données simulées
        y_values = np.unique(Y := np.array([p['y'] for p in points]))
        z_values = np.unique(Z := np.array([p['z'] for p in points]))
        Y_grid, Z_grid = np.meshgrid(y_values, z_values)

        # Graphique de gauche : Amplitude du vecteur H en 3D
        ax1 = figure.add_subplot(121, projection='3d')
        surf = ax1.plot_surface(Y_grid, Z_grid, H_total, cmap='viridis', edgecolor='k', alpha=0.8)
        ax1.set_title('Amplitude du champ magnétique |H| (y,z)')
        ax1.set_xlabel('z (m)')
        ax1.set_ylabel('y (m)')
        ax1.set_zlabel('Amplitude')
        figure.colorbar(surf, ax=ax1, shrink=0.5, aspect=10, label='|H| (A/m)')

        # Graphique de droite : Champ vectoriel H (Hy, Hz)
        
        ax2 = axes[1]
        # Tracer le champ vectoriel avec coloration basée sur la magnitude
        quiver = ax2.quiver(Y, Z, Hy, Hz, H_total, cmap='inferno', scale=100, scale_units='width', angles='xy')
        ax2.set_title("Champ vectoriel (Hy, Hz)")
        ax2.set_xlabel('y (m)')
        ax2.set_ylabel('z (m)')
        ax2.set_xlim(-0.1, 0.1)
        ax2.set_ylim(-0.1, 0.1)
        ax2.set_aspect('equal')  # Assure que les axes ont la même échelle

        # Ajouter une barre de couleur pour le quiver
        cbar = figure.colorbar(quiver, ax=ax2, shrink=0.5, aspect=10, label='Magnitude |H| (A/m)')

        layout.addWidget(canvas)
        self.tabs.addTab(tab, "Champ Amplitude et Vectoriel")

# Lancement de l'application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_window = MagneticFieldApp()
    main_window.show()
    sys.exit(app.exec_())