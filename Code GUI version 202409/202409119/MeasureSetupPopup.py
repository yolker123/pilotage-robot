import sys
from PyQt5.Qt import *
from PyQt5.QtCore import *

from bddSetupOscilloscope import *

global measure_config
global wf_img_config #dictionnaire pour les acquisitions de waveform et d'image
NB_MAX_P=8 #nombre maximal de P

class MeasureSetupPopup(QWidget):
    onClose = pyqtSignal()
    onValidate = pyqtSignal(object)
    

    def __init__(self):
        QWidget.__init__(self)
        self.ch_grid = QGridLayout()
        self.opt_grid = QGridLayout()
        self.btm_grid = QGridLayout()
        line = QFrame()
        line2 = QFrame()
        line.setFrameShape(QFrame.HLine)
        line2.setFrameShape(QFrame.HLine)
        self.ch = 0
        self.buttonCh = []
        for i in range(1, 5):
            btn = QPushButton(f"CH{i}")
            self.buttonCh.append(btn)
            def getLambda(j):
                k = int(j)
                l= lambda: self.changeChannel(k-1)
                return l
            self.buttonCh[i-1].clicked.connect(getLambda(i))
            
        
        self.validateButton = QPushButton("Validate")
        self.validateButton.clicked.connect(self.validate)
        #'pkpk', 'ampl', 'delay', 'freq',  'period', 'max',   'min, 'RMS', 'mean', 'duty cycle', 'rise28',  , 'slew', 'NBPW' 'PKS', 'PNTS'
        self.options = {
            "freq" : {
                "row" : 0,
                "name" : "Frequency"
            },
            "max" : {
                "row" : 0,
                "name" : "Maximum"
            },
            "min" : {
                "row" : 0,
                "name" : "Minimum"
            },
            "ampl":{
                "row" : 0,
                "name" : "Amplitude"
            },
        
            "RMS" : {
                "row" : 1,
                "name" : "RMS"
            },
            "pkpk" : {
                "row" : 1,
                "name" : "Pic To Pic"
            },
            "period" : {
                "row" : 1,
                "name" : "Period"
            },
            "mean" : {
                "row" : 1,
                "name" : "Mean"
            },
            
            "duty cycle" : {
                "row" : 2,
                "name" : "Duty Cycle"
            },
            "slew" : {
                "row" : 2,
                "name" : "Slew"
            },
            "NBPW" : {
                "row" : 2,
                "name" : "NBPW"
            },
            "rise28" : {
                "row" : 2,
                "name" : "rise28"
            },
        
            "PKS" : {
                "row" : 3,
                "name" : "Number of Peaks"
            },
            "PNTS" : {
                "row" : 3,
                "name" : "PKS"
            }
        
        }
        
        #channel layout
        self.layout = QVBoxLayout()
        
        ch_layout =  [
            [(QLabel("Choose the channel to measure"), 4)],
            [self.buttonCh[0], self.buttonCh[1], self.buttonCh[2], self.buttonCh[3]]
            ]
        
        #options layout
        opt_layout = self.buildRows(self.options)
        
        for _, it in self.options.items():
            sizePolicy = it["checkBox"].sizePolicy()
            sizePolicy.setRetainSizeWhenHidden(True)
            it["checkBox"].setSizePolicy(sizePolicy)
            it["checkBox"].hide()
        
        #waveform/image/validation layout
        self.cbImage = QCheckBox("Screenshot")
        self.cbWaveform = QCheckBox("Waveform")
        for it in [self.cbImage, self.cbWaveform]:
            sizePolicy = it.sizePolicy()
            sizePolicy.setRetainSizeWhenHidden(True)
            it.setSizePolicy(sizePolicy)
            it.hide()
        btm_layout = [
            [self.cbImage, self.cbWaveform, None, self.validateButton]
        ]
        
        self.ajustLayout(ch_layout)
        self.buildGrid(self.ch_grid, ch_layout)
        
        self.ajustLayout(opt_layout)
        self.buildGrid(self.opt_grid, opt_layout)
        
        self.ajustLayout(btm_layout)
        self.buildGrid(self.btm_grid, btm_layout)
        
        #popup layout
        self.layout.addLayout(self.ch_grid)
        self.layout.addWidget(line)
        self.layout.addLayout(self.opt_grid)
        self.layout.addWidget(line2)
        self.layout.addLayout(self.btm_grid)
        
        self.setLayout(self.layout)
        
        self.cbImage.clicked.connect(self.onCheckImg)
        self.cbWaveform.clicked.connect(self.onCheckWf)
        
        
        
        
        
        
    #lorsque l'image est cliquée
    def onCheckImg(self):
        wf_img_config.get(f"C{self.ch}")["img"] = self.cbImage.isChecked()
    
    #lorsque la waveforme est cliquée
    def onCheckWf(self):
        wf_img_config.get(f"C{self.ch}")["wf"] = self.cbWaveform.isChecked()
        
        
    #lorsqu'on ferme la feunetre sans valider
    def closeEvent(self, event):
        self.onClose.emit()
        super().closeEvent(event)
    
    """
     * @brief renvoie une matrice contenant les checkBox créée à partir d'un dictionnaire passé
    """
    def buildRows(self, dictionnary):
        rows = []
        for key, value in dictionnary.items():
            if value.get("row") is None:
                value["row"] = 0
            while len(rows) <= value["row"]:
                rows.append([])
            value["checkBox"] = QCheckBox(value["name"])
            value["checkBox"].clicked.connect(self.onChangeValue)
            rows[value["row"]].append(value["checkBox"])
            
        return rows
    
    """
     * ajoute des None dans un tableau contenant des widgets devant prendre plusieurs cases d'une grid
    """
    def ajustLayout(self, layout):
        x = 0
        y = 0
        while y < len(layout):
            while x < len(layout[y]):
                item = layout[y][x]
                if item is not None:
                    if isinstance(item, list) or isinstance(item, tuple):
                        if len(item) > 1:
                            nb = abs(item[1])-1
                            for _ in range(nb):
                                layout[y].insert(x+1, None)
                x += 1
            y += 1
    
    """
     * construit une grid à partir d'une matrice de widgets
    """
    def buildGrid(self, grid, layout):
        for y in range(len(layout)):
            for x in range(len(layout[y])):
                item = layout[y][x]
                if item is not None:
                    if isinstance(item, tuple) or isinstance(item, list):
                        if len(item) > 3:
                            grid.addWidget(item[0],y,x, item[2], item[1], item[3])
                        if len(item) > 2:
                            grid.addWidget(item[0],y,x, item[2], item[1])
                        elif len(item) > 1:
                            grid.addWidget(item[0],y,x, 1, item[1])
                        else:
                            grid.addWidget(item[0],y,x)
                    else:
                        grid.addWidget(item,y,x)
                else:
                    grid.addWidget(QFrame(),y,x)
    
    """
     * @brief permet de modifier les checkBox pour correspondre au channel choisi
    """
    def changeChannel(self, channel):
        self.ch = channel+1
        self.cbImage.show()
        self.cbWaveform.show()
        self.cbImage.setChecked(wf_img_config.get(f"C{self.ch}").get("img"))
        self.cbWaveform.setChecked(wf_img_config.get(f"C{self.ch}").get("wf"))

        for k_opt, opt in self.options.items():
            opt["checkBox"].setChecked(False)
            opt["checkBox"].show()
                
        for it in self.buttonCh:
            if self.buttonCh[channel] is it:
                it.setEnabled(False)
            else:
                it.setEnabled(True)
                for msr in measure_config:
                    if msr.get("ch") == f"C{channel+1}":
                        for k_opt, opt in self.options.items():
                            if k_opt == msr.get("info"):
                                opt["checkBox"].setChecked(True)
    
    """
     * @brief lorsqu'on coche/decoche une options
    """
    def onChangeValue(self):
        if self.ch < 1 or self.ch > len(self.buttonCh): return
        channel = self.ch-1
        
        for index, msr in reversed(list(enumerate(measure_config))):
            if msr.get("ch") == f"C{channel+1}":
                measure_config.pop(index)
        
        for k_opt, opt in self.options.items():
            if opt["checkBox"].isChecked():
                measure_config.append({"ch":f"C{channel+1}", "info": k_opt})
        self.changeChannel(channel)
        if len(measure_config) > NB_MAX_P:
            self.validateButton.setEnabled(False)
        else:
            self.validateButton.setEnabled(True)
    
    #lorsqu'on clique sur le bouton de validation
    def validate(self):
        self.onValidate.emit(measure_config)
        self.close()