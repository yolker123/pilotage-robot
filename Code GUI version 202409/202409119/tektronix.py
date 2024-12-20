import atexit
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

def tektronix_connection():
    try:
        device_manager = DeviceManager(verbose=False)
        atexit.register(device_manager.close)
        device_manager.visa_library = PYVISA_PY_BACKEND
        device_manager.setup_cleanup_enabled = False
        device_manager.teardown_cleanup_enabled = False

        scope: MSO6B = device_manager.add_scope(OSCILLOSCOPE_IP)
        print("Connected to:", scope.idn_string)
        return scope

    except Exception as e:
        print(f"Erreur lors de la connexion à l'oscilloscope : {e}")
        return None


# Fonction pour mesurer les tensions maximales
def measure_channels(scope, nb_measurements):
    max_values = []
    for idx in range(nb_measurements):
        max_voltage = scope.commands.measurement.meas[idx + 1].results.allacqs.mean.query()
        max_values.append(max_voltage)
    return max_values


def tektronix_set_parameters(scope, config):
    measureTypes = []
    channels = []
    for item in config:
        ch = item.get('ch')  # Extract the channel
        info = item.get('info')  # Extract the measurement info

        channel_map = {
            "C1": "CH1",
            "C2": "CH2",
            "C3": "CH3",
            "C4": "CH4"
        }
        if ch is not None and info is not None:
            # Add channel to channels list if it's not already there
            if ch not in channels:
                if ch in channel_map:
                    channels.append(channel_map[ch])

            # Add measurement type to measureTypes if it conforms to expected values
            expected_types = ["PK2PK", "MEAN", "MAXIMUM", "MINIMUM", "RMS", "PERIOD", "AMPLITUDE",
                              "FREQUENCY"]
            if info in expected_types and info not in measureTypes:
                measureTypes.append(info)

    # Configuration de l'oscilloscope
    id_map = []
    for meas_type_index, meas_type in enumerate(measureTypes):
        for channel_index, channel in enumerate(channels):
            print(f"Configuration de la mesure {meas_type} sur le canal {channel}")
            unique_id = meas_type_index * len(channels) + channel_index + 1
            print(f"ID unique : {unique_id}", meas_type, channel)
            print(f"scope:{scope}")
            scope.add_new_measurement(f"MEAS{unique_id}", meas_type, channel)
            scope.commands.measurement.meas[unique_id].source.write(channel)
            id_map.append((channel, meas_type))

    # Ouvrir le fichier pour écrire l'en-tête
    with open(filename, 'w') as f:
        headers = ["Timestamp", "x", "y", "z"]
        headers += [f"{ch}_{mt}" for ch, mt in id_map]
        f.write(", ".join(headers) + "\n")

    return len(channels) * len(measureTypes)


def tektronix_get_measures(scope, measurementNumber, x, y, z):
    # Boucle de mesure

    # print(f"Run {run + 1}/{num_runs}")
    scope.commands.acquire.state.write("OFF")
    time.sleep(0.2)
    scope.commands.acquire.mode.write("Sample")
    scope.commands.acquire.state.write("ON")
    time.sleep(0.2)
    scope.commands.acquire.state.write("OFF")

    max_values = measure_channels(scope, measurementNumber)

    # Obtenir l'horodatage actuel
    timestamp = dt.now().strftime("%Y-%m-%d %H:%M:%S")

    # Préparer la ligne à écrire
    data_line = [timestamp, x, y, z] + max_values
    data_line_str = ", ".join(map(str, data_line))

    with open(filename, 'a') as f:
        f.write(data_line_str + "\n")

    print(f"Les valeurs maximales ont été sauvegardées dans {filename}")
