#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import io
import matplotlib.pyplot as plt
import matplotlib as mpl
from PyQt5.QtGui import QPixmap, QImage
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLabel


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
    fig.text(0, 0, latex_str, fontsize=20)

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


class LatexLabelDemo(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Affichage de formule LaTeX dans QLabel")
        layout = QVBoxLayout()

        # La formule LaTeX à afficher : H avec une flèche au-dessus et un indice x.
        latex_formula = r"$\vec{H}_x$"
        pixmap = render_latex_to_pixmap(latex_formula)

        # Créer un QLabel et y placer le pixmap
        label = QLabel()
        label.setPixmap(pixmap)
        layout.addWidget(label)

        self.setLayout(layout)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    demo = LatexLabelDemo()
    demo.show()
    sys.exit(app.exec_())