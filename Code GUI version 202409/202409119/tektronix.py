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
channels = ['CH1', 'CH2', 'CH3']

filename = os.path.join(
    "C:", os.sep, "Users", "antho", "OneDrive - yncréa",
    "Documents", "BE", "Code GUI version 202409",
    "202409119", "Measure", f"MaxAmp_{current_time_str}.txt"
)

def tektronix_connection():
    rm = pyvisa.ResourceManager()
    try:
        oscilloscope = rm.open_resource(visa_address)
        print(oscilloscope.query('*IDN?'))
    except Exception as e:
        print(f"Erreur lors de la connexion à l'oscilloscope : {e}")
    finally:
        # Toujours fermer le gestionnaire de ressources PyVISA
        rm.close()

    with DeviceManager(verbose=True) as device_manager:
        # Activer la bibliothèque PyVISA-Py
        device_manager.visa_library = PYVISA_PY_BACKEND
        device_manager.setup_cleanup_enabled = True
        device_manager.teardown_cleanup_enabled = True

        scope: MSO6B = device_manager.add_scope(OSCILLOSCOPE_IP)
        print("Connected to:", scope.idn_string)
        return scope


# Fonction pour mesurer les tensions maximales
def measure_channels(scope, nb_measurements):
    max_values = []
    for idx in range(nb_measurements):
        max_voltage = scope.commands.measurement.meas[idx + 1].results.allacqs.mean.query()
        max_values.append(max_voltage)
    return max_values

def tektronix_set_parameters(scope, measurement_types):
    # Connexion à l'oscilloscope et configuration
    print("Connected to:", scope.idn_string)

    # Configuration de l'oscilloscope
    id_map = []
    for meas_type_index, meas_type in enumerate(measurement_types):
        for channel_index, channel in enumerate(channels):
            unique_id = meas_type_index * len(channels) + channel_index + 1
            scope.add_new_measurement(f"MEAS{unique_id}", meas_type, channel)
            scope.commands.measurement.meas[unique_id].source.write(channel)
            id_map.append((channel, meas_type))
    return id_map

def tektronix_get_measures(scope, id_map, measurement_types):
    # Ouvrir le fichier pour écrire l'en-tête
    with open(filename, 'w') as f:
        headers = ["Timestamp", "Run"]
        headers += [f"{ch}_{mt}" for ch, mt in id_map]
        f.write(", ".join(headers) + "\n")

    # Boucle de mesure
    for run in range(num_runs):
        print(f"Run {run + 1}/{num_runs}")
        scope.commands.acquire.state.write("OFF")
        time.sleep(0.2)
        scope.commands.acquire.mode.write("Sample")
        scope.commands.acquire.state.write("ON")
        time.sleep(0.2)
        scope.commands.acquire.state.write("OFF")

        max_values = measure_channels(scope, len(channels) * len(measurement_types))

        # Obtenir l'horodatage actuel
        timestamp = dt.now().strftime("%Y-%m-%d %H:%M:%S")

        # Préparer la ligne à écrire
        data_line = [timestamp, str(run + 1)] + max_values
        data_line_str = ", ".join(map(str, data_line))

        with open(filename, 'a') as f:
            f.write(data_line_str + "\n")

        print(f"Les valeurs maximales ont été sauvegardées dans {filename}")