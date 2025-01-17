import sys
from PyQt5.Qt import *
from PyQt5.QtCore import *
from bddSetupOscilloscope import *  # Assurez-vous que ces modules sont correctement configurés


# Supposons que measure_config et wf_img_config sont des variables globales définies ailleurs
global measure_config
global wf_img_config

# Constantes
NB_MAX_P = 8  # Nombre maximal de P

# Dictionnaire des options pour les mesures
OPTIONS = {
    "freq": {"lecroy": "Frequency", "tektronix": "FREQUENCY"},
    "max": {"lecroy": "Maximum", "tektronix": "MAXIMUM"},
    "min": {"lecroy": "Minimum", "tektronix": "MINIMUM"},
    "ampl": {"lecroy": "Amplitude", "tektronix": "AMPLITUDE"},
    "RMS": {"lecroy": "RMS", "tektronix": "RMS"},
    "pkpk": {"lecroy": "Pic To Pic", "tektronix": "PK2PK"},
    "period": {"lecroy": "Period", "tektronix": "PERIOD"},
    "mean": {"lecroy": "Mean", "tektronix": "MEAN"},
    "duty cycle": {"lecroy": "Duty Cycle"},
    "slew": {"lecroy": "Slew"},
    "NBPW": {"lecroy": "NBPW"},
    "rise28": {"lecroy": "rise28"},
    "PKS": {"lecroy": "Number of Peaks"},
    "PNTS": {"lecroy": "PKS"}
}


class MeasureSetupPopup(QWidget):
    onClose = pyqtSignal()
    onValidate = pyqtSignal(object)

    def __init__(self, oscilloName):
        super().__init__()
        self.oscilloName = oscilloName
        self.scope = None

        # Initialisation des layouts
        self.ch_grid = QGridLayout()
        self.opt_grid = QGridLayout()
        self.btm_grid = QGridLayout()

        # Lignes séparatrices
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        line2 = QFrame()
        line2.setFrameShape(QFrame.HLine)
        line2.setFrameShadow(QFrame.Sunken)

        # Initialisation des boutons de canal
        self.ch = 0
        self.buttonCh = []
        for i in range(1, 5):
            btn = QPushButton(f"CH{i}")
            self.buttonCh.append(btn)

        # Connexion des boutons de canal à changeChannel avec des connexions uniques
        for i in range(1, 5):
            # Utilisation de lambda avec un argument par défaut pour capturer la valeur actuelle de i
            self.buttonCh[i - 1].clicked.connect(lambda checked, j=i: self.changeChannel(j - 1))

        # Bouton de validation
        self.validateButton = QPushButton("Validate")
        self.validateButton.clicked.connect(self.validate)
        self.validateButton.setEnabled(True)  # Initialement activé

        # Initialisation des options de mesure
        self.options = OPTIONS  # Référence à OPTIONS
        self.opt_checkboxes = {}  # Stocke les checkboxes pour les options
        self.buildRows()  # Initialise self.opt_layout avec les checkboxes

        # Checkboxes Waveform/Image
        self.cbImage = QCheckBox("Screenshot")
        self.cbWaveform = QCheckBox("Waveform")
        self.cbImage.setVisible(False)  # Initialement cachés
        self.cbWaveform.setVisible(False)  # Initialement cachés
        self.cbImage.clicked.connect(self.onCheckImg)
        self.cbWaveform.clicked.connect(self.onCheckWf)

        # Layout inférieur contenant cbImage, cbWaveform et le bouton de validation
        btm_layout = [
            [self.cbImage, self.cbWaveform, None, self.validateButton]
        ]

        # Layout principal
        self.layout = QVBoxLayout()

        # Layout des canaux
        ch_layout = [
            [(QLabel("Choose the channel to measure"), 4)],
            [self.buttonCh[0], self.buttonCh[1], self.buttonCh[2], self.buttonCh[3]]
        ]
        self.ajustLayout(ch_layout)
        self.buildGrid(self.ch_grid, ch_layout)

        # Layout des options
        # Note : buildRows a déjà configuré self.opt_layout
        self.ajustLayout(self.opt_layout)
        self.buildGrid(self.opt_grid, self.opt_layout)

        # Layout inférieur
        self.ajustLayout(btm_layout)
        self.buildGrid(self.btm_grid, btm_layout)

        # Ajout de tous les layouts au layout principal avec les lignes séparatrices
        self.layout.addLayout(self.ch_grid)
        self.layout.addWidget(line)
        self.layout.addLayout(self.opt_grid)
        self.layout.addWidget(line2)
        self.layout.addLayout(self.btm_grid)

        self.setLayout(self.layout)

    def closeEvent(self, event):
        self.onClose.emit()
        super().closeEvent(event)

    def buildRows(self):
        """
        Créer une matrice de checkboxes à partir du dictionnaire OPTIONS, pour le layout des options.
        """
        self.opt_layout = []
        self.opt_checkboxes = {}
        line = 0
        col = 0
        for key, value in self.options.items():
            # Vérifie si oscilloName est présent dans la valeur
            if self.oscilloName in value:
                if col >= 4:  # Supposons 4 colonnes max
                    col = 0
                    line += 1
                    # Initialise la ligne si nécessaire
                if len(self.opt_layout) <= line:
                    self.opt_layout.append([])
                # Crée la checkbox
                checkbox_label = value.get("lecroy", key)
                checkbox = QCheckBox(checkbox_label)
                checkbox.clicked.connect(self.onChangeValue)
                self.opt_layout[line].append(checkbox)
                self.opt_checkboxes[key] = checkbox  # Stocke la checkbox avec la clé de l'option
                col += 1
        return self.opt_layout

    def ajustLayout(self, layout):
        """
        Ajoute des None dans un tableau contenant des widgets devant prendre plusieurs cases d'une grid.
        Note : La méthode originale n'est pas très claire ; elle est conservée pour compatibilité.
        """
        x = 0
        y = 0
        while y < len(layout):
            while x < len(layout[y]):
                item = layout[y][x]
                if item is not None:
                    if isinstance(item, (list, tuple)):
                        if len(item) > 1:
                            # Supposons que le deuxième élément du tuple est l'étendue de colonne
                            colspan = item[1]
                            if isinstance(colspan, int) and colspan > 1:
                                for _ in range(colspan - 1):
                                    layout[y].insert(x + 1, None)
                        x += 1
                    else:
                        x += 1
                else:
                    x += 1
            y += 1

    def buildGrid(self, grid, layout):
        """
        Construire une grid à partir d'une matrice de widgets.
        Args:
            grid (QGridLayout): La grille où ajouter les widgets.
            layout (list of list): La matrice contenant les widgets.
        """
        for y in range(len(layout)):
            for x in range(len(layout[y])):
                item = layout[y][x]
                if item is not None:
                    if isinstance(item, (tuple, list)) and isinstance(item[0], QWidget):
                        # Supposons le format : (widget, colspan, rowspan, alignment)
                        widget = item[0]
                        if len(item) == 4:
                            colspan, rowspan, alignment = item[1], item[2], item[3]
                            grid.addWidget(widget, y, x, rowspan, colspan, alignment)
                        elif len(item) == 3:
                            colspan, rowspan = item[1], item[2]
                            grid.addWidget(widget, y, x, rowspan, colspan)
                        elif len(item) == 2:
                            colspan = item[1]
                            grid.addWidget(widget, y, x, 1, colspan)
                        else:
                            grid.addWidget(widget, y, x)
                    elif isinstance(item, QWidget):
                        grid.addWidget(item, y, x)
                else:
                    grid.addWidget(QFrame(), y, x)  # Frame vide pour espacement

    def changeChannel(self, channel):
        """
        Modifier les checkboxes pour correspondre au canal choisi.
        Args:
            channel (int): Index du canal, 0-based.
        """
        self.ch = channel + 1  # Canaux 1-based : C1 à C4
        self.cbImage.show()
        self.cbWaveform.show()

        # Bloquer les signaux avant de modifier les checkboxes
        self.cbImage.blockSignals(True)
        self.cbWaveform.blockSignals(True)

        # Mettre à jour les checkboxes 'Screenshot' et 'Waveform' en fonction de wf_img_config
        channel_key = f"C{self.ch}"
        img_checked = wf_img_config.get(channel_key, {}).get("img", False)
        wf_checked = wf_img_config.get(channel_key, {}).get("wf", False)
        self.cbImage.setChecked(img_checked)
        self.cbWaveform.setChecked(wf_checked)

        # Débloquer les signaux après modification
        self.cbImage.blockSignals(False)
        self.cbWaveform.blockSignals(False)

        # Mettre à jour l'état des checkboxes d'options en fonction de measure_config
        for key, checkbox in self.opt_checkboxes.items():
            # Bloquer les signaux pour éviter l'appel à onChangeValue
            checkbox.blockSignals(True)
            # Vérifier si cette option est présente dans measure_config pour le canal actuel
            matched = any(
                msr.get("ch") == channel_key and msr.get("info") == self.options[key].get(self.oscilloName, key)
                for msr in measure_config
            )
            checkbox.setChecked(matched)
            checkbox.show()
            checkbox.blockSignals(False)

        # Mettre à jour les boutons de canal : désactiver le canal actuel, activer les autres
        for i, btn in enumerate(self.buttonCh):
            if i == channel:
                btn.setEnabled(False)
            else:
                btn.setEnabled(True)

    def onCheckImg(self):
        """
        Lorsque l'image est cochée/décochée.
        """
        if self.ch < 1 or self.ch > len(self.buttonCh):
            return
        channel_key = f"C{self.ch}"
        wf_img_config.setdefault(channel_key, {})
        wf_img_config[channel_key]["img"] = self.cbImage.isChecked()

    def onCheckWf(self):
        """
        Lorsque la waveform est cochée/décochée.
        """
        if self.ch < 1 or self.ch > len(self.buttonCh):
            return
        channel_key = f"C{self.ch}"
        wf_img_config.setdefault(channel_key, {})
        wf_img_config[channel_key]["wf"] = self.cbWaveform.isChecked()

    def onChangeValue(self):
        """
        Lorsque une checkbox d'option est cochée/décochée, mettre à jour measure_config.
        """
        if self.ch < 1 or self.ch > len(self.buttonCh):
            return
        channel_key = f"C{self.ch}"

        # Récupérer l'objet checkbox émetteur
        sender = self.sender()
        if not isinstance(sender, QCheckBox):
            return

        # Identifier quelle option correspond à cette checkbox
        option_key = None
        for key, checkbox in self.opt_checkboxes.items():
            if checkbox is sender:
                option_key = key
                break
        if option_key is None:
            return
        print(f"{option_key} : {self.options[option_key].get(self.oscilloName)} / {sender.isChecked()}")
        # Récupérer le nom de l'information correspondant
        info_name = self.options[option_key].get(self.oscilloName, option_key)

        if sender.isChecked():
            # Vérifier si l'ajout dépasse le maximum autorisé
            current_count = len(measure_config)
            if current_count >= NB_MAX_P:
                # Revenir en arrière et avertir l'utilisateur
                sender.blockSignals(True)
                sender.setChecked(False)
                sender.blockSignals(False)
                QMessageBox.warning(self, "Maximum Atteint", f"Impossible d'ajouter plus de {NB_MAX_P} mesures.")
                return
            # Ajouter à measure_config si non déjà présent
            if not any(msr.get("ch") == channel_key and msr.get("info") == info_name for msr in measure_config):
                measure_config.append({"ch": channel_key, "info": info_name})
        else:
            # Supprimer de measure_config toutes les entrées correspondant à ce canal et cette info
            measure_config[:] = [msr for msr in measure_config if
                                 not (msr.get("ch") == channel_key and msr.get("info") == info_name)]

        # Mettre à jour l'état du bouton de validation en fonction de measure_config
        if len(measure_config) > NB_MAX_P:
            self.validateButton.setEnabled(False)
        else:
            self.validateButton.setEnabled(True)

    def validate(self):
        """
        Lorsque le bouton de validation est cliqué.
        """
        self.onValidate.emit(measure_config)
        self.close()