"""!
 * @file        GraphicsTab.py
 * @brief       Methods for converting measurements made with lecroy to the same format as tektronix
 * @author      DEVAUX Baptiste | VOLPELLIERE Anthony
 * @version     0.1
 * @date        2025
"""

import os
import re
import sys
from typing import Dict, List, Tuple

def convert_formats(input_dir: str, output_file: str):
    """Fonction principale de conversion."""
    try:
        # Traitement de tous les fichiers
        all_data = process_directory(input_dir)
        # Écriture dans le nouveau format
        write_new_format(output_file, all_data)
    except Exception as e:
        print(f"Error while convert_formats(): {e}")
        return
    print("Conversion finished.")

def parse_coordinates(line: str) -> Tuple[float, float, float]:
    """Extrait les coordonnées x, y, z de la ligne."""
    coords = re.findall(r'(-?\d+\.?\d*)', line)
    return tuple(float(coord) for coord in coords[:3])


def parse_channel_data(lines: List[str]) -> Dict[str, str]:
    """Parse les données des channels (Mean, RMS, etc.)."""
    channel_data = {}
    for line in lines:
        if '_' in line:
            channel, value = line.strip().split(':')
            if value == "No Data Available":
                value = None
            channel_data[channel] = value
    return channel_data


def process_old_format_file(filepath: str) -> Tuple[Tuple[float, float, float], Dict[str, str]]:
    """Traite un fichier au ancien format et retourne les coordonnées et données."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    # Première ligne contient les coordonnées
    coordinates = parse_coordinates(lines[0])

    # Reste des lignes contient les données des channels
    channel_data = parse_channel_data(lines[1:])

    return coordinates, channel_data


def process_directory(input_dir: str) -> List[Tuple[Tuple[float, float, float], Dict[str, str]]]:
    """Traite tous les fichiers d'un répertoire."""
    all_data = []
    for filename in os.listdir(input_dir):
        if filename.endswith('.txt'):  # Ajustez l'extension selon vos fichiers
            filepath = os.path.join(input_dir, filename)
            measurement_data = process_old_format_file(filepath)
            all_data.append(measurement_data)
    return all_data


def write_new_format(output_file: str, data: List[Tuple[Tuple[float, float, float], Dict[str, str]]]):
    """Écrit les données dans le nouveau format CSV."""

    # Fonction pour effectuer le remplacement et convertir en majuscules
    def replace_and_uppercase(strings):
        return [re.sub(r'\bC(\d+)', r'CH\1', s).upper() for s in strings]

    # Création de l'en-tête
    channels = set()
    for _, channel_data in data:
        channels.update(channel_data.keys())

    header = ['x', 'y', 'z'] + sorted(list(channels))
    # print(f"Avant : {header}")
    header = replace_and_uppercase(header)
    # print(f"Après : {header}")

    with open(output_file, 'w') as f:
        # Écriture de l'en-tête
        f.write(','.join(header) + '\n')

        # Écriture des données
        for coordinates, channel_data in data:
            line_data = list(coordinates)
            for channel in header[3:]:  # Skip x, y, z
                line_data.append(str(channel_data.get(channel, '')))
            f.write(', '.join(map(str, line_data)) + '\n')
