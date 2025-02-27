#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import io
import matplotlib.pyplot as plt
import matplotlib as mpl
from PyQt5.QtGui import QPixmap, QImage
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QLabel, QHBoxLayout


def render_latex_to_pixmap(latex_str, dpi=150):
    """
    Génère une image à partir d'une formule LaTeX et renvoie un QPixmap.

    Parameters:
        latex_str (str): La chaîne contenant la formule LaTeX, par exemple r"$\vec{H}_x$"
        dpi (int): La résolution de l'image générée.
    Returns:
        QPixmap: L'image contenant le rendu de la formule.
    """
    # Désactiver l'utilisation de LaTeX externe pour éviter l'erreur
    mpl.rc('text', usetex=False)
    mpl.rc('font', family='serif')

    fig = plt.figure(figsize=(0.01, 0.01))
    # Placer le texte sur la figure (nous n'avons pas besoin d'axes)
    fig.text(0, 0, latex_str, fontsize=12)

    # Supprimer les axes
    plt.axis('off')

    # Sauvegarder la figure dans un tampon en mémoire
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close(fig)
    buf.seek(0)

    # Charger le tampon dans une QImage
    qimg = QImage()
    qimg.loadFromData(buf.read(), "PNG")

    # Convertir la QImage en QPixmap
    pixmap = QPixmap.fromImage(qimg)
    return pixmap


def create_vector_display(window, tab_type='2d'):
    """
    Crée un widget pour afficher les vecteurs du champ magnétique

    Parameters:
        window: La fenêtre principale
        tab_type (str): Le type d'onglet ('2d' ou '3d')
    """
    vector_widget = QWidget()
    vector_layout = QHBoxLayout(vector_widget)
    vector_layout.setContentsMargins(5, 0, 5, 0)

    # Label d'explication
    vector_label = QLabel("Vecteurs :")
    vector_layout.addWidget(vector_label)

    # Conteneur pour les formules LaTeX - sauvegardé avec un nom spécifique au type d'onglet
    formula_label = QLabel()
    vector_layout.addWidget(formula_label)

    # Stocker la référence dans l'objet window avec un nom unique basé sur le type d'onglet
    if tab_type == '2d':
        window.vector_formula_label_2d = formula_label
    else:  # '3d'
        window.vector_formula_label_3d = formula_label

    # Initialiser avec une formule par défaut
    update_vector_display(window, 'YZ', tab_type)

    return vector_widget


def update_vector_display(window, plane, tab_type='2d'):
    """
    Met à jour l'affichage des vecteurs en fonction du plan sélectionné

    Parameters:
        window: La fenêtre principale
        plane (str): Le plan sélectionné ('XY', 'XZ', ou 'YZ')
        tab_type (str): Le type d'onglet ('2d' ou '3d')
    """
    if plane == 'XY':
        latex_formula = r"$\vec{H}_x + \vec{H}_y$"
    elif plane == 'XZ':
        latex_formula = r"$\vec{H}_x + \vec{H}_z$"
    elif plane == 'YZ':
        latex_formula = r"$\vec{H}_y + \vec{H}_z$"
    else:
        latex_formula = r""

    pixmap = render_latex_to_pixmap(latex_formula)

    # Utiliser la bonne référence de label selon le type d'onglet
    if tab_type == '2d':
        if hasattr(window, 'vector_formula_label_2d'):
            window.vector_formula_label_2d.setPixmap(pixmap)
    else:  # '3d'
        if hasattr(window, 'vector_formula_label_3d'):
            window.vector_formula_label_3d.setPixmap(pixmap)