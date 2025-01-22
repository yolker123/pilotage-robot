import os
import re
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constantes
mu0 = 4 * math.pi * 1e-7  # Perméabilité magnétique du vide en T·m/A
F = 13.56e6  # Fréquence en Hz
omega = 2 * math.pi * F  # Pulsation angulaire en rad/s
S = 1e-4  # Surface de l'antenne en m² (1 cm²)


# Fonction pour extraire les amplitudes à partir du contenu d'un fichier
def extraire_amplitude_et_position(file_content, filename):
    matches = re.findall(r"C(\d+)_pkpk:(\d+\.\d+)", file_content)
    amplitudes = {f"C{axis}": float(amp) for axis, amp in matches}
    pos_match = re.search(r"logMeasureOscillo_\((-?\d+\.?\d*) (-?\d+\.?\d*) (-?\d+\.?\d*)\)", filename)

    if pos_match:
        x, y, z = map(float, pos_match.groups())
    else:
        x = y = z = 0.0  # Valeurs par défaut si la position est absente

    return amplitudes['C1'], 3, 4, x, y, z  # Bx, By, Bz, position x, y, z


# Fonction pour calculer les valeurs B, H et les propriétés du vecteur H
def calcul_champ_magnetique(amp_x, amp_y, amp_z):
    Bx = amp_x / (S * omega)
    By = amp_y / (S * omega)
    Bz = amp_z / (S * omega)

    # Calcul des composantes H
    Hx = Bx / mu0
    Hy = By / mu0
    Hz = Bz / mu0

    # Norme de H
    H_magnitude = math.sqrt(Hx ** 2 + Hy ** 2 + Hz ** 2)

    theta = math.degrees(math.acos(Hz / H_magnitude)) if H_magnitude != 0 else 0.0

    phi = math.degrees(math.atan2(Hy, Hx))

    return (Bx, By, Bz), (Hx, Hy, Hz), H_magnitude, theta, phi


# Fonction principale pour traiter les fichiers et calculer les champs
def traiter_fichiers(dossier):
    resultats = []
    for filename in os.listdir(dossier):
        if filename.startswith("logMeasureOscillo"):
            filepath = os.path.join(dossier, filename)
            with open(filepath, 'r') as file:
                file_content = file.read()

            # Extraction des amplitudes et de la position
            amp_x, amp_y, amp_z, x, y, z = extraire_amplitude_et_position(file_content, filename)

            # Calcul des champs et vecteurs
            (Bx, By, Bz), (Hx, Hy, Hz), H_magnitude, theta, phi = calcul_champ_magnetique(amp_x, amp_y, amp_z)

            # Sauvegarde des résultats pour chaque point
            resultats.append({
                'fichier': filename,
                'x': x, 'y': y, 'z': z,
                'Bx': Bx, 'By': By, 'Bz': Bz,
                'Hx': Hx, 'Hy': Hy, 'Hz': Hz,
                '|H|': H_magnitude,
                'theta': theta,
                'phi': phi
            })

    return resultats


def afficher_vecteurs_3D(resultats):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for resultat in resultats:
        x, y, z = resultat['x'], resultat['y'], resultat['z']
        Hx, Hy, Hz = resultat['Hx'], resultat['Hy'], resultat['Hz']
        H_magnitude = resultat['|H|']
        # Origine (x, y, z) et vecteur H
        ax.quiver(x, y, z, Hx, Hy, Hz, length=H_magnitude / 100, normalize=True)

        ax.scatter(x, y, z, color='red')  # Indiquer l'origine de chaque vecteur

    # Paramètres du graphique
    ax.set_xlabel('X Position (mm)')
    ax.set_ylabel('Y Position (mm)')
    ax.set_zlabel('Z Position (mm)')
    plt.title("Vecteurs H dans un espace 3D")
    plt.show()

# Fonction pour sélectionner les 6 points les plus proches du centre (excluant le point (0, 0, 0))
def selectionner_points_proches(resultats):
    # Filtrer les points qui ne sont pas à l'origine (0, 0, 0)
    resultats_sans_origine = [r for r in resultats if not (r['x'] == 0 and r['y'] == 0 and r['z'] == 0)]

    # Sélectionner les 2 points les plus proches en x
    points_x = sorted(resultats_sans_origine,
                      key=lambda r: abs(r['x']) if r['y'] == 0 and r['z'] == 0 else float('inf'))[:2]
    # Filtrer pour ne conserver qu'un seul point si les signes sont les mêmes
    if len(points_x) > 1 and (points_x[0]['x'] * points_x[1]['x'] > 0):
        points_x = [points_x[0]]

    # Sélectionner les 2 points les plus proches en y
    points_y = sorted(resultats_sans_origine,
                      key=lambda r: abs(r['y']) if r['x'] == 0 and r['z'] == 0 else float('inf'))[:2]
    # Filtrer pour ne conserver qu'un seul point si les signes sont les mêmes
    if len(points_y) > 1 and (points_y[0]['y'] * points_y[1]['y'] > 0):
        points_y = [points_y[0]]

    # Sélectionner les 2 points les plus proches en z
    points_z = sorted(resultats_sans_origine,
                      key=lambda r: abs(r['z']) if r['x'] == 0 and r['y'] == 0 else float('inf'))[:2]
    # Filtrer pour ne conserver qu'un seul point si les signes sont les mêmes
    if len(points_z) > 1 and (points_z[0]['z'] * points_z[1]['z'] > 0):
        points_z = [points_z[0]]

    # Regrouper les points sélectionnés
    points_selectionnes = points_x + points_y + points_z
    return points_selectionnes


# Fonction pour convertir les coordonnées cartésiennes en coordonnées sphériques
def convertir_cartesiennes_en_spherique(points):
    resultats_spherique = []

    for point in points:
        x = point['x']
        y = point['y']
        z = point['z']

        # Calcul de r
        r = math.sqrt(x ** 2 + y ** 2 + z ** 2)

        # Calcul de theta
        theta = math.atan2(y, x)  # Retourne l'angle en radians

        # Ajout des résultats dans la liste
        resultats_spherique.append({
            'x': x,
            'y': y,
            'z': z,
            'r': r,
            'theta': math.degrees(theta)  # Convertir en degrés pour la lisibilité
        })

    return resultats_spherique


# Exemple d'utilisation
dossier = "logMeasure"  # Dossier contenant les fichiers de mesure
resultats = traiter_fichiers(dossier)

# Affichage des résultats pour chaque point
for resultat in resultats:
    # Conversion des coordonnées cartésiennes en coordonnées sphériques
    r = math.sqrt(resultat['x'] ** 2 + resultat['y'] ** 2 + resultat['z'] ** 2)
    teta = math.degrees(math.atan2(resultat['y'], resultat['x']))  # θ en coordonnées sphériques

    print(f"Fichier: {resultat['fichier']}")
    print(f"Bx = {resultat['Bx']:.4e} T, By = {resultat['By']:.4e} T, Bz = {resultat['Bz']:.4e} T")
    print(f"Hx = {resultat['Hx']:.4e} A/m, Hy = {resultat['Hy']:.4e} A/m, Hz = {resultat['Hz']:.4e} A/m")
    print(f"|H| = {resultat['|H|']:.4e} A/m")
    print(f"θ (angle avec l'axe z) = {resultat['theta']:.2f}°")
    print(f"φ (angle dans le plan xy) = {resultat['phi']:.2f}°")

    # Affichage des coordonnées sphériques
    print(f"r = {r:.4f}")
    print(f"θ (coordonnée sphérique) = {teta:.2f}°")
    print("-------------------------------------------------")

# Sélection des points les plus proches (si nécessaire)
points_proches = selectionner_points_proches(resultats)

# Affichage des vecteurs pour les points sélectionnés
afficher_vecteurs_3D(points_proches)
