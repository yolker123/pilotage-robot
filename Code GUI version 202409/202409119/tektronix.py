"""!
 * @file        tektronix.py
 * @brief       Contains methods for connecting, configuring and measuring with the Tektronix oscilloscope.
 * @author      DEVAUX Baptiste | VOLPELLIERE Anthony | MULLOT Agathe
 * @version     0.1
 * @date        2025
"""

import atexit
from tm_devices.helpers import (
    DMConfigOptions,
    SYSTEM_DEFAULT_VISA_BACKEND,
)
import datetime
from datetime import datetime as dt
import pyvisa
from tm_devices import DeviceManager
from tm_devices.drivers import MSO6B
from tm_devices.helpers import PYVISA_PY_BACKEND
import os
import time
from bddSetupOscilloscope import *

# Configuration initiale
OSCILLOSCOPE_IP = "172.16.115.218"
visa_address = f"TCPIP::{OSCILLOSCOPE_IP}::INSTR"

pwd = os.getcwd() # Répertoire de travail actuel


global wf_img_config # Configuration globale pour les captures

class Tektronix:
    def __init__(self, main_window):
        """
        Initialise les variables nécessaires pour interagir avec l'oscilloscope Tektronix.
        """
        self.header_line = None
        self.scope = None
        self.measurement_number = None
        self.channel_measurements = None
        self.id_map = []
        self.filename = ""
        self.current_time_str = None
        self.main_window = main_window
        self.device_manager = None

    def init_connection(self):
        """
        Initialise la connexion avec l'oscilloscope, réessaye si nécessaire en cas d'échec.
        """
        self.measurement_number = 0
        while self.scope is None:
            try:
                CONFIG_OPTIONS = DMConfigOptions(
                    setup_cleanup=True,  # update the value for this option, all other options will remain untouched
                    teardown_cleanup=True,
                )
                self.device_manager = DeviceManager(verbose=False, config_options=CONFIG_OPTIONS)
                atexit.register(self.device_manager.close)

                self.device_manager.visa_library = SYSTEM_DEFAULT_VISA_BACKEND
                # Note: USB and GPIB connections are not supported with PyVISA-py backend
                scope: MSO6B = self.device_manager.add_scope("USB0::0x0699::0x0530::C071483::INSTR")

                # print("Connected to:", scope.idn_string)
                self.scope = scope
                self.main_window.validationText.append("Connected to the oscilloscope Tektronix")
            except Exception as e:
                self.main_window.validationText.append("Oscilloscope connection attempt failed, retrying...")
                self.main_window.fail_box("Oscilloscope connection attempt failed, retrying...")
                print("Oscilloscope connection attempt failed. Retrying...")
                print(f"Erreur lors de la connexion à l'oscilloscope : {e}")

    def remove_device(self):
        if self.device_manager is not None:
            self.device_manager.remove_all_devices()
            self.scope = None


    # Fonction pour mesurer les tensions maximales
    def measure_channels(self):
        """
        Mesure les mesures moyenne sur les canaux configurés et retourne les résultats.
        """
        result = []
        for idx in range(self.measurement_number):
            max_voltage = self.scope.commands.measurement.meas[idx + 1].results.currentacq.mean.query()
            channel, meas_type = self.id_map[idx]
            result.append({"channel": channel, "meas_type": meas_type, "value": max_voltage})
        return result

    def set_parameters(self, config):
        """
        Définit les paramètres d'acquisition selon la configuration fournie.
        """
        # self.scope.commands.acquire.state.write("OFF")

        self.scope.write("FPANEL:PRESS DEFaultsetup")
        self.scope.commands.acquire.mode.write("Sample")
        # Dictionnaires pour stocker les mesures par canal
        channel_map = {
            "C1": "CH1",
            "C2": "CH2",
            "C3": "CH3",
            "C4": "CH4"
        }

        expected_types = ["PK2PK", "MEAN", "MAXIMUM", "MINIMUM", "RMS", "PERIOD", "AMPLITUDE", "FREQUENCY"]

        # Dictionnaire pour stocker les mesures spécifiques à chaque canal
        self.channel_measurements = {}

        for item in config:
            ch = item.get('ch')
            info = item.get('info')

            if ch in channel_map and info in expected_types:
                mapped_channel = channel_map[ch]
                if mapped_channel not in self.channel_measurements:
                    self.channel_measurements[mapped_channel] = []

                if info not in self.channel_measurements[mapped_channel]:
                    self.channel_measurements[mapped_channel].append(info)

        # Configurer l'oscilloscope et préparer l'enregistrement des mesures
        unique_id = 1
        self.id_map = []

        self.scope.turn_channel_off("CH1")
        self.scope.turn_channel_off("CH2")
        self.scope.turn_channel_off("CH3")
        self.scope.turn_channel_off("CH4")
        for channel, measures in self.channel_measurements.items():
            print(channel)
            self.scope.turn_channel_on(channel)
            for meas_type in measures:
                print(f"Configuration de la mesure {meas_type} sur le canal {channel}")
                print(f"ID unique : {unique_id}", meas_type, channel)

                self.scope.add_new_measurement(f"MEAS{unique_id}", meas_type, channel)
                self.scope.commands.measurement.meas[unique_id].source.write(channel)
                self.id_map.append((channel, meas_type))
                unique_id += 1




        # Met à jour le nombre total de mesures configurées
        self.measurement_number = len(self.id_map)

        self.main_window.validationText.append(f"Configuration completed : {self.channel_measurements}")
        self.main_window.validationText.append(f"Total number of measures : {self.measurement_number}")
        print("Configuration completed :", self.channel_measurements)
        print(f"Total number of measures : {self.measurement_number}")

        # Préparer les en-têtes des fichiers de mesure
        headers = ["Timestamp","X","Y","Z"]  # Les en-têtes fixes
        headers += [f"{ch}_{mt}" for ch, mt in self.id_map]
        self.header_line = ",".join(headers) + "\n"

        self.main_window.validationText.append("Oscilloscope Configured")


    def get_measures(self, point, delay=1):
        """
        Lance une mesure sur les canaux configurés et enregistre les résultats dans un fichier.
        """
        point_ = point.strip("()")
        values = point_.split()
        x, y, z = map(float, values)
        # print(f"Run {run + 1}/{num_runs}")
        self.scope.commands.acquire.state.write("ON")

        time.sleep(delay)
        values = self.measure_channels()

        print(values)
        # Obtenir l'horodatage actuel
        timestamp = dt.now().strftime("%Y-%m-%d %H:%M:%S")

        # Préparer la ligne à écrire
        data_line_str = f"{timestamp}, {x}, {y}, {z}"
        for value in values:
            data_line_str += f", {value['value']}"

        # Écrire les valeurs mesurées dans le fichier
        with open(self.filename, 'a') as f:
            f.write(data_line_str + "\n")

        print(f"Les valeurs maximales ont été sauvegardées dans {self.filename}")

        # Gestion des fichiers pour les captures d'écran et formes d'ondes
        self.scope.write('FILESystem:MOUNT:DRIVE "L:;192.168.30.31;')

        if 'wf_img_config' in globals():
            for kt, it in wf_img_config.items():
                if it.get("wf"):
                    waveform_directory = "C:/Users/Public/Tektronix/TekScope/WaveForm/"
                    waveform_directory_measures = waveform_directory + f"logWaveform_tektronix_{self.current_time_str}/"
                    cwd_command_waveform = 'FILESystem:CWD "C:/Users/Public/Tektronix/TekScope/WaveForm"'
                    create_directory_command_waveform = 'FILESystem:MKDir "' + waveform_directory_measures + '"'
                    cwd_command_waveform_subfolder = 'FILESystem:CWD "' + waveform_directory_measures + '"'
                    create_directory_command_waveform2 = 'FILESystem:MKDir "' + kt + '"'
                    cwd_command_waveform_subfolder2 = 'FILESystem:CWD "' + waveform_directory_measures + kt + '"'
                    self.scope.write(cwd_command_waveform)
                    self.scope.write(create_directory_command_waveform)
                    self.scope.write(cwd_command_waveform_subfolder)
                    self.scope.write(create_directory_command_waveform2)
                    self.scope.write(cwd_command_waveform_subfolder2)
                    print("waveform activé pour le canal", kt)
                    self.captureWF_tektronix(kt, f"{x},{y},{z}")

                if it.get("img"):
                    # self.scope.commands.acquire.state.write("OFF")
                    print("image activé pour le canal", kt)
                    screenshots_directory = "C:/Users/Public/Tektronix/TekScope/Screenshots/"
                    screenshots_directory_measures = screenshots_directory + f"logScreenshot_tektronix_{self.current_time_str}/"
                    screenshot_directory_measures_channel = screenshots_directory_measures + kt + "/"
                    cwd_command_screenshot = 'FILESystem:CWD "C:/Users/Public/Tektronix/TekScope/Screenshots"'
                    create_directory_command_screenshot = 'FILESystem:MKDir "' + screenshots_directory_measures + '"'
                    cwd_command_screenshot_subfolder = 'FILESystem:CWD "' + screenshots_directory_measures + '"'
                    create_directory_command_screenshot2 = 'FILESystem:MKDir "' + kt + '"'
                    cwd_command_screenshot_subfolder2 = 'FILESystem:CWD "' + screenshot_directory_measures_channel + '"'
                    self.scope.write(cwd_command_screenshot)
                    self.scope.write(create_directory_command_screenshot)
                    self.scope.write(cwd_command_screenshot_subfolder)
                    self.scope.write(create_directory_command_screenshot2)
                    self.scope.write(cwd_command_screenshot_subfolder2)
                    self.captureScreen_tektronix(screenshot_directory_measures_channel, kt, f"{x},{y},{z}")

        else:
            print("La configuration 'wf_img_config' n'est pas définie")

        return values
    def create_file(self, hauteur, form):
        """
        Crée un fichier CSV pour sauvegarder les mesures et écrit l'en-tête des colonnes.

        Parameters:
        hauteur (str): Utilisé pour faire fonctionner l'algorithme non linéaire qui a besoin de la hauteur du robot au moment de la mesure.
        """
        # Ouvrir le fichier pour écrire l'en-tête
        hauteur = float(f"{hauteur:.2f}")
        self.current_time_str = dt.now().strftime("%Y%m%d_%H%M%S")
        if not os.path.exists(os.path.join(pwd, "Measures_tektronix")):
            os.makedirs(os.path.join(pwd, "Measures_tektronix"))
        self.filename = os.path.join(pwd, "Measures_tektronix", f"{form}_{self.current_time_str}_{hauteur}.csv")
        with open(self.filename, 'w') as f:
            f.write(self.header_line)

    def captureWF_tektronix(self, chan, nomPoint):
        """
        Capture la forme d'onde pour un canal donné et sauvegarde le fichier.

        Parameters:
        chan (str): Identifiant du canal (e.g., "C1", "C2").
        nomPoint (str): Point de mesure actuel (inclut les coordonnées x, y, z).
        """
        channel_map = {
            "C1": "CH1",
            "C2": "CH2",
            "C3": "CH3",
            "C4": "CH4"
        }

        mapped_channel = channel_map.get(chan, chan)
        print(chan)
        print(mapped_channel)
        current_time = datetime.datetime.now()
        instant = current_time.strftime('%Y-%m-%d_%H-%M-%S')
        waveform_filename = f"logWaveform_{mapped_channel}_{nomPoint}_{instant}.csv"

        self.scope.commands.save.waveform.write(f'{mapped_channel}, "{waveform_filename}"')


    def captureScreen_tektronix(self, acquisition_directory, chan, nomPoint):
            """
            Capture une capture d'écran pour un canal donné et sauve le fichier image.

            Parameters:
            acquisition_directory (str): Répertoire où stocker le fichier.
            chan (str): Identifiant du canal (e.g., "C1", "C2").
            nomPoint (str): Point de mesure actuel (inclut les coordonnées x, y, z).
            """
            current_time = datetime.datetime.now()
            instant = current_time.strftime('%Y-%m-%d_%H-%M-%S')
            screenshot_filename = f"logScreenshot_tektronix_{chan}_{nomPoint}_{instant}.png"
            # We can put the %S for the seconds to avoid rewrite the screenshot if we take two measures in the same minute
            # Specify the full path to the screenshot folder
            # Check if the screenshot folder exists, create it if it doesn't
            # Create the full path of the screenshot file
            screenshot_filepath = acquisition_directory + screenshot_filename

            write_screenshot = 'SAVE:IMAGe \"' + screenshot_filepath + '"'
            self.scope.write(write_screenshot)

            print(f"Screenshot captured and saved to {screenshot_filepath}")

