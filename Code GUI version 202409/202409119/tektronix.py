import atexit

import datetime
from datetime import datetime as dt
import pyvisa
from tm_devices import DeviceManager
from tm_devices.drivers import MSO6B
from tm_devices.helpers import PYVISA_PY_BACKEND
import os
import time

# Initialisation des paramètres
OSCILLOSCOPE_IP = "172.16.115.218"
visa_address = f"TCPIP::{OSCILLOSCOPE_IP}::INSTR"
num_runs = 1
current_time_str = dt.now().strftime("%Y%m%d_%H%M%S")
# measurement_types = ["PK2PK", "MEAN", "MAXIMUM", "MINIMUM", "RMS", "PERIOD", "AMPLITUDE", "FREQUENCY"]
# channels = ['CH1', 'CH2', 'CH3']

pwd = os.getcwd()
filename = os.path.join(pwd, "Measure", f"MaxAmp_{current_time_str}.txt")
print(f"FILENAME: {filename}")

global wf_img_config

class Tektronix:
    def __init__(self):
        self.measurement_number = None
        self.channel_measurements = None
        self.id_map = []

    def init_connection(self):
        self.scope = None
        self.measurement_number = 0
        while self.scope is None:
            try:
                device_manager = DeviceManager(verbose=False)
                atexit.register(device_manager.close)
                device_manager.visa_library = PYVISA_PY_BACKEND
                device_manager.setup_cleanup_enabled = False
                device_manager.teardown_cleanup_enabled = False

                scope: MSO6B = device_manager.add_scope(OSCILLOSCOPE_IP)
                print("Connected to:", scope.idn_string)
                self.scope = scope

            except Exception as e:
                print("Connection attempt failed. Retrying...")
                print(f"Erreur lors de la connexion à l'oscilloscope : {e}")


    # Fonction pour mesurer les tensions maximales
    def measure_channels(self):
        result = []
        for idx in range(self.measurement_number):
            max_voltage = self.scope.commands.measurement.meas[idx + 1].results.allacqs.mean.query()
            channel, meas_type = self.id_map[idx]
            result.append({"channel": channel, "meas_type": meas_type, "value": max_voltage})
        return result

    def set_parameters(self, config):
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

        for channel, measures in self.channel_measurements.items():
            for meas_type in measures:
                print(f"Configuration de la mesure {meas_type} sur le canal {channel}")
                print(f"ID unique : {unique_id}", meas_type, channel)

                self.scope.add_new_measurement(f"MEAS{unique_id}", meas_type, channel)
                self.scope.commands.measurement.meas[unique_id].source.write(channel)
                self.id_map.append((channel, meas_type))
                unique_id += 1

        # Ouvrir le fichier pour écrire l'en-tête
        with open(filename, 'w') as f:
            headers = ["Timestamp","x","y","z"]  # Les en-têtes fixes
            headers += [f"{ch}_{mt}" for ch, mt in self.id_map]
            f.write(",".join(headers) + "\n")

        # Met à jour le nombre total de mesures configurées
        self.measurement_number = len(self.id_map)

        print("Configuration complétée :", self.channel_measurements)
        print(f"Nombre total de mesures : {self.measurement_number}")

    def get_measures(self, x, y, z):
        # print(f"Run {run + 1}/{num_runs}")
        self.scope.commands.acquire.state.write("OFF")
        time.sleep(1)
        self.scope.commands.acquire.mode.write("Sample")
        self.scope.commands.acquire.state.write("ON")
        time.sleep(1)
        self.scope.commands.acquire.state.write("OFF")

        values = self.measure_channels()
        print(values)
        # Obtenir l'horodatage actuel
        timestamp = dt.now().strftime("%Y-%m-%d %H:%M:%S")

        # Préparer la ligne à écrire
        data_line_str = f"{timestamp}, {x}, {y}, {z}"
        for value in values:
            data_line_str += f", {value['value']}"

#        data_line = [timestamp, x, y, z, values[0]['value'], values[1]['value'], values[2]['value']]
 #       data_line_str = ", ".join(map(str, data_line))

        with open(filename, 'a') as f:
            f.write(data_line_str + "\n")

        print(f"Les valeurs maximales ont été sauvegardées dans {filename}")

        for kt, it in wf_img_config.items():
            if it.get("wf"):
                print("waveform activé pour le canal", kt)
                self.captureWF_tektronix("Measure", kt, f"{x},{y},{z}")
            if it.get("img"):
                print("image activé pour le canal", kt)
                self.captureScreen_tektronix("Measure", kt, f"{x},{y},{z}")

        return values

    def captureWF_tektronix(self, acquisition_directory, chan, nomPoint):

        current_time = datetime.datetime.now()
        instant = current_time.strftime('%Y-%m-%d_%H-%M-%S')
        waveform_filename = f"logWaveform_{chan}_{nomPoint}_{instant}"
        waveform_directory = acquisition_directory + "/logWaveform"
        if not os.path.exists(waveform_directory):
            os.makedirs(waveform_directory)
        waveform_filepath = os.path.join(waveform_directory, waveform_filename + ".csv")

        # save the waveform on CH1 as example.wfm
        self.scope.commands.save.waveform.write(chan, waveform_filepath)
        # # Enregistrement et tranfert de la Waveform
        # wf_dir = "D:\\Waveforms"
        # wf_name = waveform_filename
        # wf_path = wf_dir + "\\" + wf_name
        # lecroy.saveWaveform(wf_dir, chan, wf_name)
        # wf_reelpath = wf_dir + "\\" + chan + wf_name + "00000.csv"
        # lecroy.transferFile(wf_reelpath, waveform_filepath)

    def captureScreen_tektronix(self, acquisition_directory, chan, nomPoint):
            current_time = datetime.datetime.now()
            instant = current_time.strftime('%Y-%m-%d_%H-%M-%S')
            screenshot_filename = f"logScreenshot_tektronix_{chan}_{nomPoint}_{instant}.png"
            # We can put the %S for the seconds to avoid rewrite the screenshot if we take two measures in the same minute
            # Specify the full path to the screenshot folder
            screenshot_directory = acquisition_directory + "/logScreenshot"
            # Check if the screenshot folder exists, create it if it doesn't
            if not os.path.exists(screenshot_directory):
                os.makedirs(screenshot_directory)
            # Create the full path of the screenshot file
            screenshot_filepath = os.path.join(screenshot_directory, screenshot_filename)
            # Capture and save the screenshot
            self.scope.save_screenshot(
                screenshot_filepath,
                colors="INVERTED",
                local_folder=screenshot_directory,
                device_folder="./test",
                keep_device_file=True,
            )
            print(f"Screenshot captured and saved to {screenshot_filepath}")

