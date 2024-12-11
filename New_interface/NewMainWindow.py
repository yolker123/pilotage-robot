import math
import sys

import numpy as np
import serial  # pip install pyserial
import serial.tools.list_ports
from PyQt5.QtCore import *
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import (QMainWindow, QTabWidget, QVBoxLayout, QWidget,
                             QHBoxLayout, QLabel, QComboBox, QPushButton, QSlider, QFrame)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from scipy import special
from scipy.constants import mu_0
from scipy.integrate import quad

from magneticFieldSimulation import MagneticFieldSimulation

"""
 * @brief Search the list of USB devices connected to the computer
"""


def find_USB_device():
    myports = [tuple(p) for p in list(serial.tools.list_ports.comports())]
    usb_port_list = [p[0:2] for p in
                     myports]  # Prendre p[0], p[1], no need p[2], because myports has 3 string parts.use print to see the value of myports.
    # for ex. myports[0]=('COM4', 'Périphérique série USB (COM4)', 'USB VID:PID=8087:0ACA SER=05022016'), only take first 2 parts.
    # p[0:2] take only p[0] ET P[1] . no p[2].   0<=i<2
    return usb_port_list


class MagneticFieldSimulation_brouillon(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax = fig.add_subplot(111, projection='3d')
        self.R = 0.1
        self.I = 1.0
        self.N = 10
        self.d = 0.01
        self.k = 10
        self.S = 1
        super().__init__(fig)

        self.resultats = []  # Points à afficher
        self.new_points_buffer = []  # Buffer pour les nouveaux points

    def magnetic_field_helix(self, x, y, z):
        """Calcul du champ magnétique pour un point (x, y, z)."""

        def dB_component(phi, z0, component):
            x0 = self.R * np.cos(phi)
            y0 = self.R * np.sin(phi)
            r_vec = np.array([x - x0, y - y0, z - z0])
            r_mag = np.linalg.norm(r_vec)
            if r_mag == 0:
                return 0
            dL = np.array([-self.R * np.sin(phi), self.R * np.cos(phi), 0])
            dB_vec = (mu_0 * self.I / (4 * np.pi)) * np.cross(dL, r_vec) / (r_mag ** 3)
            return dB_vec[component]

        def dB_component_Bx(phi, z0):
            return dB_component(phi, z0, 0)

        def dB_component_By(phi, z0):
            return dB_component(phi, z0, 1)

        def dB_component_Bz(phi, z0):
            return dB_component(phi, z0, 2)

        B_total = np.array([0.0, 0.0, 0.0])
        for i in range(self.N):
            z0 = i * self.d

            # Diviser l'intégrale en plusieurs sous-intervalles pour éviter les singularités
            Bx1, _ = quad(dB_component_Bx, 0, np.pi / 2, args=(z0,))
            Bx2, _ = quad(dB_component_Bx, np.pi / 2, np.pi, args=(z0,))
            Bx3, _ = quad(dB_component_Bx, np.pi, 3 * np.pi / 2, args=(z0,))
            Bx4, _ = quad(dB_component_Bx, 3 * np.pi / 2, 2 * np.pi, args=(z0,))

            By1, _ = quad(dB_component_By, 0, np.pi / 2, args=(z0,))
            By2, _ = quad(dB_component_By, np.pi / 2, np.pi, args=(z0,))
            By3, _ = quad(dB_component_By, np.pi, 3 * np.pi / 2, args=(z0,))
            By4, _ = quad(dB_component_By, 3 * np.pi / 2, 2 * np.pi, args=(z0,))

            Bz1, _ = quad(dB_component_Bz, 0, np.pi / 2, args=(z0,))
            Bz2, _ = quad(dB_component_Bz, np.pi / 2, np.pi, args=(z0,))
            Bz3, _ = quad(dB_component_Bz, np.pi, 3 * np.pi / 2, args=(z0,))
            Bz4, _ = quad(dB_component_Bz, 3 * np.pi / 2, 2 * np.pi, args=(z0,))

            # Additionner les résultats des sous-intervalles
            B_total += np.array([Bx1 + Bx2 + Bx3 + Bx4,
                                 By1 + By2 + By3 + By4,
                                 Bz1 + Bz2 + Bz3 + Bz4])

        return B_total

    def generate_simulated_points(self):
        """Simule les points pour le champ magnétique."""
        self.resultats = [{'x': x, 'y': y, 'z': z} for x in np.linspace(-0.1, 0.1, 5)
                          for y in np.linspace(-0.1, 0.1, 5)
                          for z in np.linspace(-0.1, 0.1, 5)]

    def afficher_vecteurs_3D(self):
        """Affiche les nouveaux points périodiquement (toutes les 250 ms)."""
        if not self.new_points_buffer:
            return

        # Ajouter les nouveaux points au graphique
        for point in self.new_points_buffer:
            x, y, z = point['x'], point['y'], point['z']
            Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
            self.ax.quiver(x, y, z, Hx, Hy, Hz, color='b', length=0.01, normalize=True)

        self.new_points_buffer.clear()
        self.draw()

    def calcul_r_teta(self, x, y, z):
        r = math.sqrt(x ** 2 + y ** 2 + z ** 2)
        teta = math.degrees(math.acos(z / r)) if r != 0 else 0.0
        return r, teta

    def selectionner_points_proches(self):
        resultats_sans_origine = [r for r in self.resultats if not (r['x'] == 0 and r['y'] == 0 and r['z'] == 0)]
        points_selectionnes = []

        for axis in ['x', 'y', 'z']:
            axis_points = sorted(resultats_sans_origine,
                                 key=lambda r: abs(r[axis]) if all(r[a] == 0 for a in 'xyz' if a != axis) else float(
                                     'inf'))[:2]
            if len(axis_points) > 1 and (axis_points[0][axis] * axis_points[1][axis] > 0):
                axis_points = [axis_points[0]]
            points_selectionnes.extend(axis_points)
        return points_selectionnes

    def moyenne_I(self, points_proches):
        I_valeurs = []
        for point in points_proches:
            r, teta = self.calcul_r_teta(point['x'], point['y'], point['z'])
            H_r, H_theta = point['Hx'], point['Hy']  # Simplified for now
            I_r = self.calculer_I(H_r, r, teta)
            I_theta = self.calculer_I_Htetha(H_theta, r, teta)
            I_valeurs.append((I_r + I_theta) / 2 if I_r and I_theta else (I_r or I_theta))
        return sum(I_valeurs) / len(I_valeurs) if I_valeurs else None

    def calculer_I(self, H_r, r, theta):
        facteur1 = 1j / (self.k ** 2 * r ** 2)
        facteur2 = 1 / (self.k ** 3 * r ** 3)
        denominateur = self.S * (self.k ** 3) * (facteur1 + facteur2) * math.cos(math.radians(theta))
        return abs(2 * math.pi * H_r / denominateur) if denominateur else None

    def calculer_I_Htetha(self, H_theta, r, theta):
        sin_theta = np.sin(np.radians(theta))
        facteur = -1 / (self.k * r) + 1j / (self.k ** 2 * r ** 2) + 1 / (self.k ** 3 * r ** 3)
        return abs(4 * np.pi * H_theta / (self.S * self.k ** 3 * facteur * sin_theta)) if sin_theta else None

    def calculer_H_total(self, x, y, z, I):
        Hr, Htheta = self.calculer_Hr(x, y, z, I), self.calculer_Htheta(x, y, z, I)
        return math.sqrt(Hr ** 2 + Htheta ** 2)

    def calculer_Hr(self, x, y, z, I):
        r, theta = self.calcul_r_teta(x, y, z)
        facteur = (1j / (self.k ** 2 * r ** 2)) + (1 / (self.k ** 3 * r ** 3))
        Hr = (I * self.S * self.k ** 3 / (2 * math.pi)) * facteur * math.cos(math.radians(theta)) * np.exp(
            -1j * self.k * r)
        return abs(Hr)

    def calculer_Htheta(self, x, y, z, I):
        r, theta = self.calcul_r_teta(x, y, z)
        sin_theta = np.sin(np.radians(theta))
        facteur = (-1 / (self.k * r) + 1j / (self.k ** 2 * r ** 2) + 1 / (self.k ** 3 * r ** 3))
        Htheta = (I * self.S * self.k ** 3 / (4 * math.pi)) * facteur * sin_theta * np.exp(-1j * self.k * r)
        return abs(Htheta)

    def augmenter_resolution(self, points, I_moyen):
        interpolated_points = []
        for i in range(len(points) - 1):
            p1, p2 = points[i], points[i + 1]
            interpolated_points.extend(self.interpoler_points(p1, p2, I_moyen, self.resolution))
        self.points_haute_resolution = interpolated_points

    def interpoler_points(self, p1, p2, I, resolution):
        interpolated = []
        for i in range(1, resolution):
            fraction = i / resolution
            x = p1['x'] + (p2['x'] - p1['x']) * fraction
            y = p1['y'] + (p2['y'] - p1['y']) * fraction
            z = p1['z'] + (p2['z'] - p1['z']) * fraction
            Hr = self.calculer_Hr(x, y, z, I)
            Htheta = self.calculer_Htheta(x, y, z, I)
            # put in cartesian coordinates
            Hx = Hr * np.cos(np.radians(y))
            Hy = Hr * np.sin(np.radians(y))
            Hz = Htheta
            interpolated.append({'x': x, 'y': y, 'z': z, 'Hx': Hx, 'Hy': Hy, 'Hz': Hz})
        return interpolated


class Worker(QObject):
    point_generated = pyqtSignal(dict)

    def __init__(self, simulation):
        super().__init__()
        self.simulation = simulation
        self.running = True

    def run(self):
        """Générer des points en continu, en arrière-plan."""
        for x in np.linspace(-0.1, 0.1, 5):
            for y in np.linspace(-0.1, 0.1, 5):
                for z in np.linspace(-0.1, 0.1, 5):
                    if not self.running:
                        return
                    Bx, By, Bz = self.simulation.magnetic_field_helix(x, y, z)
                    Hx, Hy, Hz = Bx / mu_0, By / mu_0, Bz / mu_0
                    point = {'x': x, 'y': y, 'z': z, 'Hx': Hx, 'Hy': Hy, 'Hz': Hz}
                    self.simulation.new_points_buffer.append(point)
                    self.point_generated.emit(point)
                    # Calcul de chaque point et ajout dans le buffer

    def stop(self):
        """Arrêter le calcul du worker."""
        self.running = False


class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.ax = fig.add_subplot(111, projection='3d')
        super().__init__(fig)
        self.points = []  # Liste pour stocker les points aléatoires

    def plot_initial_surface(self):
        """Trace une surface initiale avec SciPy."""
        x = np.linspace(-5, 5, 100)
        y = np.linspace(-5, 5, 100)
        x, y = np.meshgrid(x, y)
        z = special.j0(np.sqrt(x ** 2 + y ** 2))

        self.ax.plot_surface(x, y, z, cmap='viridis', alpha=0.6)
        self.ax.set_title("Visualisation 3D avec SciPy")
        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        self.ax.set_zlabel("Z")

    def add_random_point(self):
        """Ajoute un point aléatoire au graphique."""
        x = np.random.uniform(-5, 5)
        y = np.random.uniform(-5, 5)
        z = special.j0(np.sqrt(x ** 2 + y ** 2))
        self.points.append((x, y, z))

        # Effacer les anciens points avant de redessiner
        self.ax.cla()
        self.plot_initial_surface()

        # Ajouter les nouveaux points
        for px, py, pz in self.points:
            self.ax.scatter(px, py, pz, color='red', s=50)

        # Rafraîchir le canvas pour afficher les nouveaux points
        self.draw()


class MeasureInterface(QMainWindow):

    def __init__(self):
        super().__init__()
        self.interface = None
        self.buttons = {"oscilloscope" : {}, "robot" : {}, "demo" : None}
        self.style_h1 = "font-size: 14pt; font-weight: bold; color: black;"
        self.style_h2 = "font-size: 10pt; font-weight: bold; color: Red;"
        self.speed_value = 1000  # Valeur initiale de la vitesse
        self.createInterface()

    def getInterface(self) -> QWidget:
        return self.interface

    def createInterface(self):
        widget = QWidget()
        layout = QHBoxLayout()
        widget.setLayout(layout)

        left_ui_widget = self.create_left_ui_widget()
        right_ui_widget = self.create_right_ui_widget()

        layout.addWidget(left_ui_widget, stretch=4)
        layout.addWidget(right_ui_widget, stretch=1)

        self.init_ui_status_bar()
        self.connect_signals()
        self.interface = widget

    def create_left_ui_widget(self) -> QFrame:
        left_panel = QFrame()
        left_layout = QHBoxLayout(left_panel)

        # Création du canvas pour le graphique 3D
        self.simulation = MagneticFieldSimulation_brouillon(self, width=8, height=6, dpi=100)
        left_layout.addWidget(self.simulation)

        # self.measure_tab_main_layout.addWidget(left_panel, stretch=4)
        return left_panel

    def create_right_ui_widget(self) -> QFrame:
        # === PANNEAU DROIT (Contrôles) ===
        right_panel = QFrame()
        right_layout = QVBoxLayout(right_panel)

        oscilloscope_layout = self.create_oscilloscope_widget()
        robot_section = self.create_robot_widget()

        right_layout.addWidget(oscilloscope_layout)
        right_layout.addWidget(robot_section)

        self.buttons["demo"] = QPushButton("Démo")
        self.buttons["demo"].setStyleSheet("""
                    QPushButton {
                        background-color: orange;
                        color: white;
                        font-weight: bold;
                    }
                    QPushButton:hover {
                        background-color: black;
                    }
                """)
        right_layout.addWidget(self.buttons["demo"])

        # Ajouter un espace extensible
        right_layout.addStretch()

        # self.measure_tab_main_layout.addWidget(right_panel, stretch=1)
        return right_panel

    def create_oscilloscope_widget(self) -> QWidget:
        widget = QWidget()
        oscilloscope_layout = QVBoxLayout()
        widget.setLayout(oscilloscope_layout)

        conn_title = QLabel("Oscilloscope")
        conn_title.setStyleSheet(self.style_h1)
        oscilloscope_layout.addWidget(conn_title)

        self.buttons["oscilloscope"]["connection"] = QPushButton("Connecter Oscilloscope")
        self.buttons["oscilloscope"]["connection"].setStyleSheet("""
                    QPushButton {
                        background-color: #3b82f6;
                        color: white;
                    }
                    QPushButton:hover {
                        background-color: #2563eb;
                    }
                """)
        oscilloscope_layout.addWidget(self.buttons["oscilloscope"]["connection"])

        # self.right_layout.addLayout(oscilloscope_layout)
        return widget

    def create_robot_widget(self) -> QWidget:
        robot_title = QLabel("Robot")
        robot_title.setStyleSheet(self.style_h1)

        # Sélection du port
        self.port_combo = QComboBox()
        self.port_combo.addItems(["COM1", "COM2", "COM3", "COM4"])
        # self.portlist = find_USB_device()

        # self.port_combo.addItems([p[1] for p in self.portlist])  # take only the p[1] (the second value) of each p in portlist to the list items[].
        # self.port_combo.append("Choice COM Port KEOLABS Robx")  # add end list one item

        robot_list = QComboBox()
        robot_list.addItems(["Robot 5 axis", "Robot 6 axis"])

        self.buttons["robot"]["connection"] = QPushButton("Robot connection")
        self.buttons["robot"]["connection"].setStyleSheet("""
                    QPushButton {
                        background-color: #3b82f6;
                        color: white;
                    }
                    QPushButton:hover {
                        background-color: #2563eb;
                    }
                """)

        self.buttons["robot"]["home"] = QPushButton("Home position")
        self.buttons["robot"]["home"].setEnabled(False)

        # --- Section Vitesse ---
        speed_frame = QFrame()
        speed_layout = QVBoxLayout(speed_frame)

        speed_title = QLabel("Vitesse")
        speed_title.setStyleSheet(self.style_h2)
        speed_layout.addWidget(speed_title)

        # Création d'un layout horizontal pour le slider et la valeur
        speed_control_layout = QHBoxLayout()

        self.speed_slider = QSlider(Qt.Horizontal)
        self.speed_slider.setMinimum(100)
        self.speed_slider.setMaximum(2000)
        self.speed_slider.setValue(self.speed_value)
        speed_control_layout.addWidget(self.speed_slider, stretch=4)

        # Amélioration du label de vitesse
        self.speed_value_label = QLabel(f"{self.speed_value} mm/s")
        self.speed_value_label.setStyleSheet("""
                           QLabel {
                               background-color: #e5e7eb;
                               border-radius: 4px;
                               min-width: 100px;
                               text-align: center;
                               font-weight: bold;
                           }
                       """)
        self.speed_value_label.setAlignment(Qt.AlignCenter)
        # speed_control_layout.addWidget(self.speed_value_label, stretch=1)

        speed_layout.addWidget(self.speed_value_label, stretch=4)
        speed_layout.addLayout(speed_control_layout)

        # Ajout des labels de vitesse avec un style amélioré
        speed_labels = QHBoxLayout()
        speed_labels_text = ["Very slow", "Slow", "Medium", "Fast", "Very fast"]
        for text in speed_labels_text:
            label = QLabel(text)
            label.setStyleSheet("""
                               QLabel {
                                   color: #4b5563;
                                   font-size: 12px;
                               }
                           """)
            label.setAlignment(Qt.AlignCenter)
            speed_labels.addWidget(label)
        speed_layout.addLayout(speed_labels)

        # --- Section Contrôle du Mouvement ---
        movement_frame = QFrame()
        movement_layout = QVBoxLayout(movement_frame)

        move_title = QLabel("Contrôle du Mouvement")
        move_title.setStyleSheet(self.style_h2)
        movement_layout.addWidget(move_title)

        # Sélection du pas
        step_layout = QHBoxLayout()
        step_layout.addWidget(QLabel("Pas:"))
        self.step_combo = QComboBox()
        self.step_combo.addItems(["1 mm", "5 mm", "10 mm"])
        step_layout.addWidget(self.step_combo)
        movement_layout.addLayout(step_layout)

        # Contrôles X, Y, Z
        for axis in ['X', 'Y', 'Z']:
            axis_layout = QHBoxLayout()
            axis_layout.addWidget(QLabel(f"{axis}:"))

            minus_btn = QPushButton("-")
            minus_btn.setFixedSize(40, 40)
            axis_layout.addWidget(minus_btn)

            coord_label = QLabel("0.0 mm")
            coord_label.setAlignment(Qt.AlignCenter)
            axis_layout.addWidget(coord_label)

            plus_btn = QPushButton("+")
            plus_btn.setFixedSize(40, 40)
            axis_layout.addWidget(plus_btn)

            setattr(self, f"{axis.lower()}_minus_btn", minus_btn)
            setattr(self, f"{axis.lower()}_plus_btn", plus_btn)
            setattr(self, f"{axis.lower()}_label", coord_label)

            movement_layout.addLayout(axis_layout)

        # Boutons Sauvegarder/Charger
        save_load_layout = QHBoxLayout()
        self.save_btn = QPushButton("Sauvegarder")
        self.load_btn = QPushButton("Charger")
        save_load_layout.addWidget(self.save_btn)
        save_load_layout.addWidget(self.load_btn)
        movement_layout.addLayout(save_load_layout)

        # --- Bouton d'arrêt d'urgence ---
        self.emergency_btn = QPushButton("ARRÊT D'URGENCE")
        self.emergency_btn.setStyleSheet("""
                    QPushButton {
                        background-color: #dc2626;
                        color: white;
                        font-weight: bold;
                    }
                    QPushButton:hover {
                        background-color: #b91c1c;
                    }
                """)
        widget = QWidget()
        robot_layout = QVBoxLayout()
        widget.setLayout(robot_layout)

        robot_layout.addWidget(robot_title)
        robot_layout.addWidget(self.port_combo)
        robot_layout.addWidget(robot_list)
        robot_layout.addWidget(self.buttons["robot"]["connection"])
        robot_layout.addWidget(self.buttons["robot"]["home"])
        robot_layout.addWidget(speed_frame)
        robot_layout.addWidget(movement_frame)
        robot_layout.addWidget(self.emergency_btn)

        # self.right_layout.addLayout(robot_layout)
        return widget

    def init_ui_status_bar(self):
        # --- Barre de statut ---
        self.status_bar = self.statusBar()
        self.status_bar.setStyleSheet("""
                    QStatusBar {
                        background-color: #1f2937;
                        color: white;
                    }
                """)
        self.status_bar.showMessage("État: Déconnecté")

    def connect_signals(self):
        # Connexion des boutons aux slots
        self.buttons["oscilloscope"]["connection"].clicked.connect(self.toggleOscilloscope)
        self.buttons["robot"]["connection"].clicked.connect(self.toggleRobot)
        self.buttons["robot"]["home"].clicked.connect(self.goHome)
        # self.emergency_btn.clicked.connect(self.emergencyStop)
        # self.save_btn.clicked.connect(self.savePosition)
        self.load_btn.clicked.connect(self.loadPosition)
        self.step_combo.currentIndexChanged.connect(self.updateStep)
        self.speed_slider.valueChanged.connect(self.updateSpeed)
        self.buttons["demo"].clicked.connect(self.start_demo)

    @pyqtSlot()
    def start_demo(self):
        # Worker et thread pour le calcul en parallèle
        self.worker = Worker(self.simulation)
        self.thread = QThread()
        self.worker.moveToThread(self.thread)

        self.worker.point_generated.connect(self.on_point_generated)

        # Démarrer le calcul immédiatement dans un thread séparé
        self.start_simulation()

        # Mise à jour du graphique toutes les 250 ms (pour l'affichage fluide)
        self.update_timer = QTimer(self)
        self.update_timer.timeout.connect(self.simulation.afficher_vecteurs_3D)
        self.update_timer.start(500)  # Mise à jour tous les 250 ms (4 fois par seconde)

    def start_simulation(self):
        """Démarrer le thread du worker pour le calcul en parallèle."""
        self.thread.started.connect(self.worker.run)
        self.thread.start()

    def on_point_generated(self, point):
        """Fonction appelée lorsqu'un point est généré."""
        self.simulation.new_points_buffer.append(point)

    def closeEvent(self, event):
        """Arrêter proprement le thread lors de la fermeture."""
        # self.worker.stop()
        # self.worker.thread.quit()
        # self.worker.thread.wait()
        super().closeEvent(event)

    @pyqtSlot(int)
    def updateSpeed(self, value):
        """
        Met à jour la vitesse du robot et l'affichage
        """
        self.speed_value = value
        # Mise à jour du label avec formatage amélioré
        self.speed_value_label.setText(f"{value} mm/s")

        # Mise à jour de la couleur en fonction de la vitesse
        if value < 500:
            color = "#22c55e"  # vert pour vitesse lente
        elif value < 1000:
            color = "#eab308"  # jaune pour vitesse moyenne
        else:
            color = "#ef4444"  # rouge pour vitesse rapide

        self.speed_value_label.setStyleSheet(f"""
                  QLabel {{
                      padding: 5px 10px;
                      background-color: {color};
                      color: white;
                      border-radius: 4px;
                      min-width: 100px;
                      text-align: center;
                      font-weight: bold;
                  }}
              """)

        # Mise à jour du statut
        self.status_bar.showMessage(f"Vitesse mise à jour : {value} mm/s", 2000)

    @pyqtSlot()
    def toggleOscilloscope(self):
        connected = self.oscilloscope_btn.text() == "Connecter Oscilloscope"
        self.oscilloscope_btn.setText("Déconnecter Oscilloscope" if connected else "Connecter Oscilloscope")
        self.oscilloscope_btn.setStyleSheet("""
              QPushButton {
                  background-color: #22c55e;
                  color: white;
              }
          """ if connected else """
              QPushButton {
                  background-color: #3b82f6;
                  color: white;
              }
          """)

    @pyqtSlot()
    def toggleRobot(self):
        connected = self.robot_connection_btn.text() == "Connecter Robot"
        self.robot_connection_btn.setText("Déconnecter Robot" if connected else "Connecter Robot")
        self.robot_connection_btn.setStyleSheet("""
              QPushButton {
                  background-color: #22c55e;
                  color: white;
              }
          """ if connected else """
              QPushButton {
                  background-color: #3b82f6;
                  color: white;
              }
          """)
        self.robot_home_position.setEnabled(connected)
        self.status_bar.showMessage(f"État: {'Connecté' if connected else 'Déconnecté'}")

    @pyqtSlot()
    def goHome(self):
        # Implémenter le retour à la position initiale
        pass

    @pyqtSlot()
    def emergencyStop(self):
        # Implémenter l'arrêt d'urgence
        self.status_bar.showMessage("ARRÊT D'URGENCE ACTIVÉ", 5000)

    @pyqtSlot()
    def savePosition(self):
        # Implémenter la sauvegarde de position
        pass

    @pyqtSlot()
    def loadPosition(self):
        # Implémenter le chargement de position
        pass

    @pyqtSlot(int)
    def updateStep(self, index):
        # Mettre à jour la valeur du pas
        step_values = [1, 5, 10]
        self.step_value = step_values[index]


class RobotInterface(QMainWindow):
    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        # Configuration de base de la fenêtre
        self.setWindowTitle("Robot Bureau d'Étude")

        self.showMaximized()
        self.setStyleSheet("""
            QPushButton {
                padding: 8px;
                border-radius: 4px;
                background-color: #f3f4f6;
            }
            QPushButton:hover {
                background-color: #e5e7eb;
            }
            QLabel {
                font-size: 14px;
            }
            QFrame {
                border-radius: 8px;
                padding: 10px;
            }
        """)

        # Widget central
        central_widget = QTabWidget()
        self.setCentralWidget(central_widget)

        measure_interface = MeasureInterface()
        # calculs_widget = self.createCalculsWidget()

        central_widget.addTab(measure_interface.getInterface(), "Measurement")
        # central_widget.addTab(calculs_widget, "Calculs")



if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = RobotInterface()
    window.show()
    sys.exit(app.exec_())
