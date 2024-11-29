import os
import re
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import cmath

import numpy as np
# Constantes
c = 3e8  # Vitesse de la lumière en m/s
mu0 = 4 * math.pi * 1e-7  # Perméabilité magnétique du vide en T·m/A
F = 13.56e6  # Fréquence en Hz
omega = 2 * math.pi * F  # Pulsation angulaire en rad/s
S = 1e-4  # Surface de l'antenne en m² (1 cm²)
k = omega / c  # Nombre d'onde en rad/m

# Fonction pour extraire les amplitudes à partir du contenu d'un fichier
def extraire_amplitude_et_position(file_content, filename):
    matches = re.findall(r"C(\d+)_pkpk:(\d+\.\d+)", file_content)
    amplitudes = {f"C{axis}": float(amp) for axis, amp in matches}
    pos_match = re.search(r"logMeasureOscillo_\((-?\d+\.?\d*) (-?\d+\.?\d*) (-?\d+\.?\d*)\)", filename)
    
    if pos_match:
        x, y, z = map(float, pos_match.groups())
    else:
        x = y = z = 0.0  # Valeurs par défaut si la position est absente

    return amplitudes['C1'], 0, 0, x, y, z  # Bx, By, Bz, position x, y, z


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
    H_magnitude = math.sqrt(Hx**2 + Hy**2 + Hz**2)

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
            (Bx, By, Bz), (Hx, Hy, Hz), H_magnitude, Htheta, Hphi = calcul_champ_magnetique(amp_x, amp_y, amp_z)

            # Sauvegarde des résultats pour chaque point
            resultats.append({
                'fichier': filename,
                'x': x, 'y': y, 'z': z,
                'Bx': Bx, 'By': By, 'Bz': Bz,
                'Hx': Hx, 'Hy': Hy, 'Hz': Hz,
                '|H|': H_magnitude,
                'Htheta': Htheta,
                'Hphi': Hphi
            })

    return resultats

def afficher_vecteurs_3D(resultats):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for resultat in resultats:
        x, y, z = resultat['x'], resultat['y'], resultat['z']
        print(resultat)
        Hx, Hy, Hz = resultat['Hx'], resultat['Hy'], resultat['Hz']
        H_magnitude = resultat['|H|']
        # Origine (x, y, z) et vecteur H
        ax.quiver(x, y, z, Hx, Hy, Hz, length=H_magnitude/100, normalize=True)
        
        ax.scatter(x, y, z, color='red')  # Indiquer l'origine de chaque vecteur

    # Paramètres du graphique
    ax.set_xlabel('X Position (mm)')
    ax.set_ylabel('Y Position (mm)')
    ax.set_zlabel('Z Position (mm)')
    plt.title("Vecteurs H dans un espace 3D")
    plt.show()
    
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

def calcul_r_teta(x, y, z):
    r = math.sqrt(x ** 2 + y ** 2 + z ** 2)
    teta = math.degrees(math.atan2(y, x))
    return r, teta

# Calculer I basé sur H_r
def calculer_I(H_r, r, theta):
    numerateur = 2 * math.pi * H_r
    k = omega / 3e8  # Nombre d'onde
    facteur1 = 1j / (k**2 * r**2)
    facteur2 = 1 / (k**3 * r**3)
    denominateur = S * (k**3) * (facteur1 + facteur2) * math.cos(math.radians(theta))
    
    if denominateur == 0:
        return None
    I = numerateur / denominateur
    return abs(I)

# Calculer I basé sur H_theta
def calculer_I_Htetha(H_theta, r, theta):
    # Vérifier si theta est proche de 0 pour éviter division par zéro dans sin(theta)
    if np.isclose(theta, 0):
        return None  # Ignorer ce calcul si theta est proche de 0

    k = omega / 3e8  # Nombre d'onde
    sin_theta = np.sin(np.radians(theta))

    # Calcul de la partie complexe de l'équation
    facteur = -1 / (k * r) + 1j / (k**2 * r**2) + 1 / (k**3 * r**3)

    # Calcul de I
    I = (4 * np.pi * H_theta) / (S * k**3 * facteur * sin_theta)
    return abs(I)

# Fonction pour calculer la moyenne des valeurs de I
def moyenne_I(points_proches):
    I_valeurs = []
    for point in points_proches:
        H_r = point['|H|'] 
        H_theta = point['Htheta']
        r, theta = point['r'], point['tetha']
        
        # Calculer I pour H_r et H_theta
        I_r = calculer_I(H_r, r, theta)
        I_theta = calculer_I_Htetha(H_theta, r, theta)
        print("I_r:", I_r)
        print("I_theta:", I_theta)
        # Si les deux valeurs existent, on prend la moyenne pour ce point
        if I_r is not None and I_theta is not None:
            I_moyenne_point = (I_r + I_theta) / 2
            I_valeurs.append(I_moyenne_point)
        elif I_r is not None:
            I_valeurs.append(I_r)
        elif I_theta is not None:
            I_valeurs.append(I_theta)

    # Calcul de la moyenne globale
    if I_valeurs:
        print("I_valeurs:", I_valeurs)
        return sum(I_valeurs) / len(I_valeurs)
    else:
        return None
    
    
def calculer_Hr(x, y, z, I_moyenne):
    r = math.sqrt(x**2 + y**2 + z**2)
    
    # Éviter la division par zéro
    if r == 0:
        return 0  # Retourne 0 ou une autre valeur par défaut si r est nul
    
    theta = math.degrees(math.acos(z / r)) if r != 0 else 0

    facteur = (1j / (k**2 * r**2)) + (1 / (k**3 * r**3))
    Hr = (I_moyenne * S * k**3 / (2 * math.pi)) * facteur * math.cos(math.radians(theta)) * np.exp(-1j * k * r)

    return abs(Hr)  # Prendre le module pour obtenir une valeur réelle

# Fonction pour calculer Hθ pour un point donné (x, y, z)
def calculer_Htheta(x, y, z, I_moyenne):
    r = math.sqrt(x**2 + y**2 + z**2)
    theta = math.degrees(math.acos(z / r)) if r != 0 else 0

    sin_theta = np.sin(np.radians(theta))
    if np.isclose(sin_theta, 0):
        return 0  # Éviter une division par zéro si theta est proche de 0

    facteur = (-1 / (k * r) + 1j / (k**2 * r**2) + 1 / (k**3 * r**3))
    Htheta = (I_moyenne * S * k**3 / (4 * math.pi)) * facteur * sin_theta * np.exp(-1j * k * r)
    
    return abs(Htheta)   

# Fonction pour calculer H_total à partir de Hr et Hθ
def calculer_H_total(x, y, z, I):
    Hr = calculer_Hr(x, y, z, I)
    Htheta = calculer_Htheta(x, y, z, I)
    H_total = math.sqrt(Hr**2 + Htheta**2)
    return H_total
    

# Exemple d'utilisation
dossier = "logMeasure"  # Dossier contenant les fichiers de mesure
resultats = traiter_fichiers(dossier)           
# Sélection des points les plus proches
points_proches = selectionner_points_proches(resultats)
for point in points_proches:
    r, tetha = calcul_r_teta(point['x'], point['y'], point['z'])
    point['r'] = r
    point['tetha'] = tetha
    
    
I_moyenne = moyenne_I(points_proches)
print("I_moyenne:", I_moyenne)
print("Points les plus proches:", points_proches)
# Affichage des vecteurs pour les points sélectionnés

# Exemple d'utilisation
x, y, z = 10.0, 0.0, 0.0  # Position dans l'espace en mm

Htotal_val = calculer_H_total(x, y, z, I_moyenne)

print(f"H_total = {Htotal_val:.4e} A/m")

afficher_vecteurs_3D(points_proches)
# Afficher les résultats
for point in resultats:
    print(f"Fichier: {point['fichier']}")
    print(f"Bx = {point['Bx']:.4e} T, By = {point['By']:.4e} T, Bz = {point['Bz']:.4e} T")
    print(f"Hx = {point['Hx']:.4e} A/m, Hy = {point['Hy']:.4e} A/m, Hz = {point['Hz']:.4e} A/m")
    print(f"|H| = {point['|H|']:.4e} A/m")
    print(f"Hθ (angle avec l'axe z) = {point['Htheta']:.2f}°")
    print(f"Hφ (angle dans le plan xy) = {point['Hphi']:.2f}°")
    #Affichage des coordonnées sphériques
    print(f"r = {r:.4f}")
    print(f"θ (coordonnée sphérique) = {tetha:.2f}°")
    print("-------------------------------------------------")
    

# retourner q'un point si meme signe dans point proche
# fonction prend un point en entrée et retourne r et teta (pour r, passer x y et z, formule rac(x²+y²+z²)), pour tetat, atan2(y,x)

# avec Hr = H total -> isoler I pour chaque point et faire la moyenne/ OK mais incohérent
# Hteta = arctan(Hy/Hx) -> 12 valeurs de I -> moyenne
# Faire une fonction pour calculer Hr et Hteta avec notre I
# Calculer Htotal = sqrt(Hr² + Hteta²)

# Visualiser les points générer
# creer une variable resolution, qui augmente la densité de points

# Fonction pour interpoler les points entre deux points existants
def interpoler_points(p1, p2, resolution, I_moyenne):
    points_interpolés = []
    for i in range(1, resolution):
        fraction = i / resolution
        x = p1['x'] + (p2['x'] - p1['x']) * fraction
        y = p1['y'] + (p2['y'] - p1['y']) * fraction
        z = p1['z'] + (p2['z'] - p1['z']) * fraction

        # Calculer |H|, r et tetha pour les points interpolés
        r = math.sqrt(x**2 + y**2 + z**2)
        tetha = math.degrees(math.atan2(y, x))
        H_magnitude = calculer_Hr(x, y, z, I_moyenne)

        # Créer le point interpolé avec tous les attributs nécessaires
        point_interpolé = {
            'x': x,
            'y': y,
            'z': z,
            '|H|': H_magnitude,
            'r': r,
            'tetha': tetha,
            'Hx': p1['Hx'],  # ou interpoler en fonction des besoins
            'Hy': p1['Hy'],
            'Hz': p1['Hz'],
            'Bx': p1['Bx'],
            'By': p1['By'],
            'Bz': p1['Bz']
        }
        points_interpolés.append(point_interpolé)
    
    return points_interpolés


# Fonction pour augmenter la résolution du champ magnétique
def augmenter_resolution(resultats, resolution, I_moyenne):
    points_haute_resolution = []

    for i in range(len(resultats) - 1):
        p1 = resultats[i]
        p2 = resultats[i + 1]

        # Ajouter le point initial
        points_haute_resolution.append(p1)
        print("p1:", p1)
        # Ajouter les points interpolés
        points_interpolés = interpoler_points(p1, p2, resolution, I_moyenne)
        points_haute_resolution.extend(points_interpolés)

    # Ajouter le dernier point
    points_haute_resolution.append(resultats[-1])

    return points_haute_resolution

# Exemple d'utilisation avec une résolution choisie
resolution = 5
points_haute_resolution = augmenter_resolution(points_proches, resolution, I_moyenne)

# Affichage des vecteurs pour les points interpolés avec la haute résolution
afficher_vecteurs_3D(points_haute_resolution)

# Afficher les résultats pour vérifier
for point in points_haute_resolution:
    print(f"x = {point['x']:.2f}, y = {point['y']:.2f}, z = {point['z']:.2f}")
    print(f"|H| = {point['|H|']:.4e} A/m")
    print(f"r = {point['r']:.4f}, tetha = {point['tetha']:.2f}°")
    print("-------------------------------------------------")
