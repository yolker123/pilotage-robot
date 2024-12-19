"""!
 * @file        MainWindows.py
 * @brief       Principal page to launch the QT application
 * @author      Clément Rouvier | Lisa Duterte
 * @version     0.1
 * @date        2023
"""
# pour ouvrir RobXManager, ouvrir un invite de commande dans le dossier KEOLABS/RobX/Bin 
#et entrer la commande "Rob6xManager P S192.168.30.34 DENSO VS-AV6"

import sys
import serial,serial.tools.list_ports  #pip install pyserial
# from Robot import *
# from aquisitionLoops import *
# import pythoncom
# import win32com.client
import threading
from concurrent.futures import Future
import time
from MeasureSetupPopup import *
from NewWindowWithData import MagneticFieldApp
# from oscilloscopeAcquisition import *
from PyQt5 import QtCore, QtGui, QtWidgets

speed = 2000
COM = "COM22"     # Port du cable#
initOscilloscope_flag = False   
initRobot = False
idOscilloscope = 'USB0::0x05FF::0x1023::3561N16324::INSTR'  #Indentifant de l'osiloscope:.. /support/about/remote
global workers #contient un dictionnaire pour les threads permettant de garantir leur durée de vie et de stocker l'objet COM pour KEOLABS
global emergencyStop

"""
 * @brief Search the list of USB devices connected to the computer
"""
def find_USB_device():
        myports = [tuple(p) for p in list(serial.tools.list_ports.comports())] 
        usb_port_list = [p[0:2] for p in myports]     #Prendre p[0], p[1], no need p[2], because myports has 3 string parts.use print to see the value of myports.
                                                   #for ex. myports[0]=('COM4', 'Périphérique série USB (COM4)', 'USB VID:PID=8087:0ACA SER=05022016'), only take first 2 parts.
                                                   # p[0:2] take only p[0] ET P[1] . no p[2].   0<=i<2
        return usb_port_list

"""
 * @brief permet de tester si l'objet passé est convertible en un nombre flottant positif/ to verify if the parameter is positive, example the size of the cube.
"""  
def is_positive_float(element: any) -> bool:
    #If you expect None to be passed:
    if element is None: 
        return False
    try:
        f = float(element)
        return f > 0
    except Exception:
        return False
""" 
 * to verify the parameter if it can be covert to float type
"""     
def is_float(element: any) -> bool:
    if element is None: 
        return False
    try:
        f = float(element)
        return True
    except Exception:
        return False


"""
 * @brief annotation pour qu'une méthode soit lancée en parallèle.
 * L'interface graphique ne doit pas être modifiée hors du thread principal. 
 * Il est possible de lever le signal executeFunction de MainWindow pour executer une fonction dans le thread principal.
"""
def asynchrone(fn):
    def wrapper(self, *args, **kwargs):
        future = Future()
        def Worker(future, fn, this, *args, **kwargs):
            workers[threading.current_thread().ident] = [None, threading.current_thread()]
            result = None
            pythoncom.CoInitialize()
            try:
                result =fn(this, *args, **kwargs)
            except Exception as e:
                print("error : ", e)
                future.set_exception(e)
            
            pythoncom.CoUninitialize()
            future.set_result(result)
            try:
                workers.pop(threading.current_thread().ident)
            except Exception as e:
                print("error del : ", e)

        worker = threading.Thread(target=Worker, args=(future, fn, self, *args),kwargs=kwargs)
        """multithread creation"""
        
        worker.start()
        
        return future
    return wrapper

"""
 * @brief Class managing the QT application 
"""
class MainWindow(QMainWindow):
    
    executeFunction = pyqtSignal(object) #signal permettant de faire executer une fonction dans le thread de la feunetre principal
    
    def handleExecuteFunction(self, func):
        func()

    def testButton(self):
        print("testButton")

    def __init__(self):
        super().__init__()
        self.executeFunction.connect(self.handleExecuteFunction)
        self.distanceMove=5
        #change the title of the window
        self.setWindowTitle("Robot bureau d'étude")
        self.RobID = "DENSO"
        # self.robot = createRobot(self.RobID)
        # ------------------- Definition of central point for acquisition ---
        self.x_ptr = 0   # point de reference  --prt. begin point of robort arm.
        self.y_ptr = 0     # unit mm
        self.z_ptr = 0   #localisation initial xyz
        
        #point de reference fixé pour la figure One Point
        self.fixed_x = 0 
        self.fixed_y = 0
        self.fixed_z = 0
        
        
        current_time = datetime.datetime.now()
        #emplacement des acquisitions pour la figure One Point
        self.point_log_dir = f"logAcquisition/points_{current_time.strftime('%Y-%m-%d_%H-%M-%S')}"
        
        #prédicat pour savoir si le robot peut bouger
        self.canMove = False
        
        #prédicat pour savoir si l'oscilloscope peut mesurer
        self.canMeasure = False

        self.showMaximized()

        self.portlist=find_USB_device()
        self.items=[p[1] for p in self.portlist]    #take only the p[1] (the second value) of each p in portlist to the list items[].
        self.items.append("Choice COM Port KEOLABS Robx")    #    add end list one item
        
        centralArea = QWidget()
        self.layout_central_V = QVBoxLayout(centralArea)
        # centralArea.setStyleSheet("background: '#ECECEC'")     #  #ECECEC color gray

        tabs = QTabWidget()
        tabs.addTab(centralArea, "Mesures")
        mfa = MagneticFieldApp()
        tabs.addTab(mfa, "Algo")
        self.setCentralWidget(tabs)  # The defaut value of centralArea is set in self.resize(1200, 800)

        self.layout_central_top = QHBoxLayout()
        self.layout_central_V.addLayout(self.layout_central_top, stretch=2)
        self.left_widget = QWidget()
        self.layout_central_top.addWidget(self.left_widget, stretch=2)

        self.right_widget = QWidget()
        self.layout_right_V = QVBoxLayout(self.right_widget)
        self.layout_central_top.addWidget(self.right_widget, stretch=1)
        # right_widget.setGeometry(800, 0, 400, 1000)

        bottom_widget = QWidget()
        self.layout_central_V.addWidget(bottom_widget, stretch=1)
        # bottom_widget.setGeometry(10, 800, 790, 200)
        
        #------------------------COM Ports----------------------
        port_widget = QWidget()
        self.layout_right_V.addWidget(port_widget)
        # port_widget.setGeometry(0, 0, 400, 80)
        port_widget.selectlbl = QLabel("Select port:")  #creat a lable for select tool bar
        #label
        port_widget.typeBox=QComboBox()  # creat toolbar box
        port_widget.typeBox.addItems(self.items) # add items to toolbar
        port_widget.typeBox.setCurrentIndex(port_widget.typeBox.count()-1) 
        # set current index to the last item in the typebox (count()-1 -->last one)
       
        fields=QGridLayout()   #Dispose les widgets dans une grille.
        fields.addWidget(port_widget.selectlbl,0,0)  #  line 0 colone 0
        fields.addWidget(port_widget.typeBox,1,0)   # line 1 colone 0
        port_widget.setLayout(fields)  # set is to apply this fields variable
        #signal connection 
        COM = port_widget.typeBox.currentIndexChanged.connect(self.valuePort)
        ##when release mouse, go to function valuePort in the class to get the current value
        
        
        #------------------------Speed----------------------        
        speed_widget = QWidget()
        self.layout_right_V.addWidget(speed_widget)
        # speed_widget.setGeometry(0, 80, 400, 90)
                
        #Slider speed
        title_speed = QLabel("Select speed:")
        slider = QSlider(Qt.Horizontal, speed_widget) # creat a slider in speed_widget
        labelVerySlow = QLabel("Very Slow")
        labelSlow = QLabel("  Slow")
        labelMedium = QLabel("  Medium")
        labelFast = QLabel("   Fast")
        labelVeryFast = QLabel("     Very Fast")
        grid = QGridLayout()                    #Dispose les widgets dans une grille.
        grid.addWidget(title_speed, 0, 0, 1, 5)  #start location from line 0,colone0 occupy the 1st and 5th colones (from 1 to 5)
        grid.addWidget(slider, 1, 0, 1, 5)       # from line 1,colone0 occupy the 1st and 5th colones (from 1 to 5)
        grid.addWidget(labelVerySlow, 2, 0)
        grid.addWidget(labelSlow, 2, 1)  # line 2 colone 1
        grid.addWidget(labelMedium, 2, 2)
        grid.addWidget(labelFast, 2, 3)
        grid.addWidget(labelVeryFast, 2, 4)  # totally 5 colones from 0 to 4
        slider.setMinimum(100)  # Set slider value from 1000 TO 2000 for speed
        slider.setMaximum(2000)
        slider.setValue(speed)  #initial value speed 2000
        speed_widget.setLayout(grid)  # set is to apply this grid variable
        #connection du signal
        slider.sliderReleased.connect(self.valueChanged)     # save value to valuechanged function when release the slider with mouse. 
        #when release mouse, go to function valueChanged in the class

        # fonction autosetup odcilloscope


        #------------------------Buttons----------------------

        # ----- MOCHEEEEE ------
        btn_widgetAxes2 = QWidget()
        self.layout_right_V.addWidget(btn_widgetAxes2)
        grid_btnAxes2 = QGridLayout() # creat gridlayout
        self.buttons = None
        self.buttonLecroy = QPushButton('Lecroy')
        self.buttonLecroy.clicked.connect(self.lecroy)
        self.buttonTektronix = QPushButton('Tektronix')
        self.buttonTektronix.clicked.connect(self.tektronix)
        grid_btnAxes2.addWidget(self.buttonLecroy , 0, 0)  # add widgets to ths gridlayout
        grid_btnAxes2.addWidget(self.buttonTektronix, 0, 1)
        btn_widgetAxes2.setLayout(grid_btnAxes2)
        # ----------------------

        btn_widget = QWidget()
        self.layout_right_V.addWidget(btn_widget)
        # btn_widget.setGeometry(0, 220, 400, 350)  #begin from (0,250) in right_widget, length 400 and height 250
        
        #button to connect oscilloscope
        self.buttonConnectOscilloscope = QPushButton(' Oscilloscope COM Initialization')
        self.buttonConnectOscilloscope.clicked.connect(self.initOscilloscope)     #when click buttonConnectOscilloscope, go to function initOscilloscope in the class, so self.initO.... , if function outside this class, no 'self.'
        self.buttonConnectOscilloscope.setEnabled(False)

        #button to open popup to setup measures
        self.buttonSetupOscilloscope = QPushButton('Setup Oscilloscope Parameters')
        self.buttonSetupOscilloscope.clicked.connect(self.setupMeasure)
        self.buttonSetupOscilloscope.setEnabled(False)

        # button to select 5 axes robot
        self.buttonRobot5Axes = QPushButton('Robot 5 axes')
        self.buttonRobot5Axes.clicked.connect(self.robotAxis5)

        # button to select 6 axes robot
        self.buttonRobot6Axes = QPushButton('Robot 6 axes')
        self.buttonRobot6Axes.clicked.connect(self.robotAxis6)

        btn_widgetAxes = QWidget()
        self.layout_right_V.addWidget(btn_widgetAxes)
        # btn_widgetAxes.setGeometry(0, 170, 400, 50)
        grid_btnAxes = QGridLayout() # creat gridlayout
        grid_btnAxes.addWidget(self.buttonRobot5Axes, 0, 0)  # add widgets to ths gridlayout
        grid_btnAxes.addWidget(self.buttonRobot6Axes, 0, 1)
        btn_widgetAxes.setLayout(grid_btnAxes)       # set this gridlayout(grid_btnAxes) to btn_widgetAxes widget which is defined in line 143
        # setLayout is to display the buttons, if not the button you set will not display on the screen
        
        #button to connect robot
        self.buttonConnectRobot = QPushButton('Robot Initialization')
        self.buttonConnectRobot.clicked.connect(self.initialization)
        self.buttonConnectRobot.setEnabled(False)
        
        #button to move robot to home position (defined by the robot)
        self.buttonHomePosition = QPushButton('Home Position ')
        self.buttonHomePosition.clicked.connect(lambda: self.moveHome())  # Go to function S_Gotohome, it's outside this class, so no self. S_GoToHome() make robort arm to initial position
        self.buttonHomePosition.setEnabled(False)
        
        #button to disconnect robot
        self.buttonDisconnectRobot = QPushButton('Close Robot communication')
        self.buttonDisconnectRobot.clicked.connect(self.disconnect)
        self.buttonDisconnectRobot.setEnabled(False)


        grid_btn = QGridLayout()  # creat gridlayout named grid_btn

        grid_btn.addWidget(self.buttonConnectOscilloscope, 0, 0)
        grid_btn.addWidget(self.buttonSetupOscilloscope, 1, 0)
        grid_btn.addWidget(self.buttonConnectRobot, 2, 0)
        grid_btn.addWidget(self.buttonHomePosition, 3, 0)
        grid_btn.addWidget(self.buttonDisconnectRobot, 4, 0)

        btn_widget.setLayout(grid_btn)               # display the buttons by putting grid_btn gridlayout to btn_widget widget
        
        #------------------------Control widget----------------------
        self.ctrl_widget = QWidget(bottom_widget)
        # self.ctrl_widget.setGeometry(0, 0, 750, 150)  #begin from (0,250) in right_widget, length 400 and height 250.
        
        self.ctrl_title = QLabel("Robot movement control:")

        # AXE X
        self.ctrl_X = QLabel("X")
        
        #button + on x
        self.xPlus = QPushButton('+')
        self.xPlus.clicked.connect(lambda: self.moveXplus())
        self.xPlus.setEnabled(False)
        
        #button - on x
        self.xMinus = QPushButton('-')
        self.xMinus.clicked.connect(lambda: self.moveXminus())
        self.xMinus.setEnabled(False)

        
        # AXE Y
        self.ctrl_Y = QLabel("Y")
        
        #button + on Y
        self.yPlus = QPushButton('+')
        self.yPlus.clicked.connect(lambda: self.moveYplus())
        self.yPlus.setEnabled(False)
        
        #button - on Y
        self.yMinus = QPushButton('-')
        self.yMinus.clicked.connect(lambda: self.moveYminus())
        self.yMinus.setEnabled(False)


        # AXES Z
        self.ctrl_Z = QLabel("Z")
        
        #button + on Z
        self.zPlus = QPushButton('+')
        self.zPlus.clicked.connect(lambda: self.moveZplus())
        self.zPlus.setEnabled(False)
        
        #button - on z
        self.zMinus = QPushButton('-')
        self.zMinus.clicked.connect(lambda: self.moveZminus())
        self.zMinus.setEnabled(False)
        
        #button to show actual coordinates of the robot
        self.buttonShowCoords = QPushButton("Show (x,y,z)")
        self.buttonShowCoords.clicked.connect(self.showCoords)
        self.buttonShowCoords.setEnabled(False)
        
        #button to fix actual coordinates of the robot as reference point for One Point figure
        self.buttonFixPoint = QPushButton("Fix Point")
        self.buttonFixPoint.clicked.connect(self.fixPoint)
        self.buttonFixPoint.setEnabled(False)
        
        #button to save actual coordinates of the robot in file
        self.buttonSavePoint = QPushButton("Save")
        self.buttonSavePoint.clicked.connect(self.savePoint)
        self.buttonSavePoint.setEnabled(False)
        
        #button to load coordinates of the robot from file. It moves the robot
        self.buttonOpenPoint = QPushButton("Open")
        self.buttonOpenPoint.clicked.connect(self.openPoint)
        self.buttonOpenPoint.setEnabled(False)
        
        #hide fix button but keep space
        for w in [self.buttonFixPoint]:
            sizePolicy = w.sizePolicy()
            sizePolicy.setRetainSizeWhenHidden(True)
            w.setSizePolicy(sizePolicy)
            w.hide()
        

        ### Define step indicator's Label
        self.stepLabel = QLabel("Step:")
    
        # Define step selector's Combo Box
        self.stepSelector = QComboBox()
        self.stepSelector.currentIndexChanged.connect(self.changeStep)
        self.stepSelector.addItems(self.getListStep()[0])

        grid_ctrl = QGridLayout()  # creat gridlayout named grid_ctrl
        grid_ctrl.addWidget(self.ctrl_title, 0, 0, 1, 3)

        # Line 2 laying out steps widgets
        grid_ctrl.addWidget(self.stepLabel, 1, 0, 1, 1) 
        grid_ctrl.addWidget(self.stepSelector, 1, 2)

        # Line 3 laying out x coordinate widgets
        grid_ctrl.addWidget(self.xMinus, 2, 0)
        grid_ctrl.addWidget(self.ctrl_X, 2, 1)
        grid_ctrl.addWidget(self.xPlus, 2, 2)

        # Line 4 laying out y coordinate widgets
        grid_ctrl.addWidget(self.yMinus, 3, 0)
        grid_ctrl.addWidget(self.ctrl_Y, 3, 1)
        grid_ctrl.addWidget(self.yPlus, 3, 2)

        # Line 5 laying out z coordinate widgets
        grid_ctrl.addWidget(self.zMinus, 4, 0)
        grid_ctrl.addWidget(self.ctrl_Z, 4, 1)
        grid_ctrl.addWidget(self.zPlus, 4, 2)

        # Other widgets
        grid_ctrl.addWidget(self.buttonShowCoords, 4, 3)
        grid_ctrl.addWidget(self.buttonSavePoint, 3, 3)
        grid_ctrl.addWidget(self.buttonOpenPoint, 2, 3)

        # Define widgets for figure grid
        self.labelFigureSelector =  QLabel("Choose Measurement Method :")
        self.figureSelector = QComboBox()
        self.figureSelector.addItems(self.getFigureList()[0])
        self.figureSelector.currentIndexChanged.connect(lambda:self.onFigureSelected())
        self.bottom_hbox = QHBoxLayout()
        self.bottom_hbox.addLayout(grid_ctrl)
        # Creating an empyty label for add a new line to the layout (for align figureGrid and grid_ctrl's layout)
        self.empty1Label = QLabel("")


        # Laying out the figureGrid
        self.figureGrid = QGridLayout()
        self.figureGrid.addWidget(self.labelFigureSelector, 0, 0)
        self.figureGrid.addWidget(self.figureSelector, 1, 0)
        self.figureGrid.addWidget(self.buttonFixPoint, 2, 0)
        
        #button to launch a figure
        self.buttonLaunch = QPushButton("Launch")
        self.buttonLaunch.clicked.connect(self.figureLaunch)
        
        self.figureGrid.addWidget(self.buttonLaunch, 3, 0)
        # Laying out an empty line
        self.figureGrid.addWidget(self.empty1Label, 4, 0)

        self.label_x_cube_size = QLabel("x size :")
        self.label_y_cube_size = QLabel("y size :")
        self.label_z_cube_size = QLabel("z size :")
        
        #fields for custom cube figure
        self.x_cube_size = QLineEdit()
        self.x_cube_size.textChanged.connect(self.udpdateEnable)
        self.y_cube_size = QLineEdit()
        self.y_cube_size.textChanged.connect(self.udpdateEnable)
        self.z_cube_size = QLineEdit()
        self.z_cube_size.textChanged.connect(self.udpdateEnable)
        
        self.label_x_cube_points = QLabel("x points :")
        self.label_y_cube_points = QLabel("y points :")
        self.label_z_cube_points = QLabel("z points :")
        
        self.x_cube_points = QLineEdit()
        self.x_cube_points.textChanged.connect(self.udpdateEnable)
        self.y_cube_points = QLineEdit()
        self.y_cube_points.textChanged.connect(self.udpdateEnable)
        self.z_cube_points = QLineEdit()
        self.z_cube_points.textChanged.connect(self.udpdateEnable)
        
        self.unit_label = QLabel("Unit: mm")
        
        customCubeFieldWidgets = [self.label_x_cube_size, self.label_y_cube_size, self.label_z_cube_size, self.x_cube_size, self.y_cube_size, self.z_cube_size, self.label_x_cube_points, self.label_y_cube_points, self.label_z_cube_points, self.x_cube_points, self.y_cube_points, self.z_cube_points, self.unit_label]
        #cache les champs et garde l'espace
        for w in customCubeFieldWidgets:
            sizePolicy = w.sizePolicy()
            sizePolicy.setRetainSizeWhenHidden(True)
            w.setSizePolicy(sizePolicy)
            w.hide()
        
        self.figureGrid.addWidget(QFrame(), 0, 1)
        self.figureGrid.addWidget(self.unit_label, 0, 2)
        self.figureGrid.addWidget(self.label_x_cube_size, 1, 1)
        self.figureGrid.addWidget(self.label_y_cube_size, 2, 1)
        self.figureGrid.addWidget(self.label_z_cube_size, 3, 1)
        self.figureGrid.addWidget(self.label_x_cube_points, 1, 3)
        self.figureGrid.addWidget(self.label_y_cube_points, 2, 3)
        self.figureGrid.addWidget(self.label_z_cube_points, 3, 3)
        
        self.figureGrid.addWidget(self.x_cube_size, 1, 2)
        self.figureGrid.addWidget(self.y_cube_size, 2, 2)
        self.figureGrid.addWidget(self.z_cube_size, 3, 2)
        self.figureGrid.addWidget(self.x_cube_points, 1, 4)
        self.figureGrid.addWidget(self.y_cube_points, 2, 4)
        self.figureGrid.addWidget(self.z_cube_points, 3, 4)
        
        self.bottom_hbox.addLayout(self.figureGrid)
        self.ctrl_widget.setLayout(self.bottom_hbox) 
        
        #alias of the custom cube figure fields to reuse it with the custom cylindre figure
        self.label_radius = self.label_x_cube_size
        self.field_radius = self.x_cube_size
        self.label_height = self.label_y_cube_size
        self.field_height = self.y_cube_size
        self.label_z_points = self.label_x_cube_points
        self.field_z_points = self.x_cube_points
        self.label_theta_points = self.label_y_cube_points
        self.field_theta_points = self.y_cube_points
        
        #update the figure selected
        self.onFigureSelected()
        #------------------------Reference point coordinates----------------------
        text_widget = QWidget()
        self.layout_right_V.addWidget(text_widget)
        # text_widget.setGeometry(0, 560, 550, 80)
       
        grid = QGridLayout()
       
        text_widget.setLayout(grid)  # set grid QGridlayout to text_widget QWidget 
        
        
        # ptr_widget = QWidget(self.right_widget)
        # ptr_widget.setGeometry(0, 630, 400, 80)
        
        # clicked OK button then goto validateCoordPt function
        
        #------------------------Validatation Text Box----------------------
        validTextBox_widget = QWidget()
        self.layout_right_V.addWidget(validTextBox_widget)
        # validTextBox_widget.setGeometry(10, 710, 350, 100)
        validTextBox_widget.setStyleSheet("border: 1px solid black;")
    
        grid = QVBoxLayout()   #creat a QVBoxLayout--> organizes your widgets vertically in this window.no need use QGridLayout()
        global validationText
        validationText = QTextEdit()  # creat textbox
        validationText.setEnabled(False)  # False --> cant be edit, just display
        grid.addWidget(validationText)  # add validationText (test box) to the VBoxlayout named grid
        validTextBox_widget.setLayout(grid)
        
        #------------------------Emergency STOP Button----------------------
        STOP_widget = QWidget()
        self.layout_right_V.addWidget(STOP_widget)
        # STOP_widget.setGeometry(10, 810, 350, 200)
        STOP_btn = QPushButton('Emergency STOP')
        STOP_btn.clicked.connect(self.EmergencyStop)  # when click the button, go to function EmergencyStop in the class
        STOP_btn.setStyleSheet("background-color: red;")
        STOP_btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        STOP_grid = QVBoxLayout()
        STOP_grid.addWidget(STOP_btn)
        STOP_widget.setLayout(STOP_grid)
        self.udpdateEnable() #udpdateEnable() is a class defined in line 929

    # ------ MOCHE ---------
    def lecroy(self):
        self.buttonLecroy.setEnabled(False)
        self.buttonTektronix.setEnabled(True)
        self.buttonConnectOscilloscope.setEnabled(True)

    def tektronix(self):
        self.buttonTektronix.setEnabled(False)
        self.buttonLecroy.setEnabled(True)
        self.buttonConnectOscilloscope.setEnabled(True)
    # ----------------------
    """
     * @brief activate the buttons
    """ 
    def robotAxis5(self):
        self.RobID = "5axes"
        self.x_ptr = 0 
        self.y_ptr = 419   
        self.z_ptr = -170 
        self.buttonRobot5Axes.setEnabled(False)
        self.buttonRobot6Axes.setEnabled(True)
        self.buttonConnectRobot.setEnabled(True)
    
        
    """
     * @brief activate the buttons
    """ 
    def robotAxis6(self):
        self.RobID = "DENSO"
        self.x_ptr = 0 
        self.y_ptr = 100   
        self.z_ptr = -170 
        self.buttonRobot6Axes.setEnabled(False)
        self.buttonRobot5Axes.setEnabled(True)
        self.buttonConnectRobot.setEnabled(True)

 
    
    """
     * @brief Initialise the connection to the oscilloscope
    """ 
    def initOscilloscope(self):
        rm = oscilloscopeConnection(idOscilloscope)   # oscilloscopeConnection is a function in oscilloscopeAcquisition
        #display the action on the specific Text Box
        validationText.setText("Oscilloscope Connection: " + str(rm.list_resources()))  # settext is to write text in validationText
        
        #variable set to true in order to enable the setupOscilloscope fonction
        self.buttonSetupOscilloscope.setEnabled(True)
        self.buttonConnectOscilloscope.setEnabled(False)

            

    """
     * @brief Launch the measureOscillo function of the oscilloscope
    """
    def MeasureOscillo(self):
        current_time = datetime.datetime.now()
        log_dir = f"logAcquisition/measure_{current_time.strftime('%Y-%m-%d_%H-%M-%S')}"
        getAcquisition("(0 0 0)", 0, log_dir)
        validationText.setText("Autosetup Done")
   
    """
     * @brief Initialize the Robot position 
    """ 
    def initialization(self):
        import time
        self.robot = createRobot(self.RobID)
        err = self.robot.OpenComm(COM)
        if err != 0:
            print(f"Error failed to connect robot {err}")
            return
        
        self.moveHome()

        global initRobot
        initRobot = True
        self.buttonRobot5Axes.setEnabled(False)
        self.buttonRobot6Axes.setEnabled(False)
        self.buttonDisconnectRobot.setEnabled(True)
        self.buttonConnectRobot.setEnabled(False)
        self.canMove = True
        self.udpdateEnable()
            

    """
     * @brief Launch the NFC Acquisition 
    """  
    @asynchrone
    def nfc(self):
        current_time = datetime.datetime.now()
        log_dir = f"logAcquisition/nfc_{current_time.strftime('%Y-%m-%d_%H-%M-%S')}"#directory of the figure 
        #display the action on the specific Text Box
        def gui1():
            #update gui
            self.canMove = False
            self.buttonSetupOscilloscope.setEnabled(False)
            self.udpdateEnable()
            validationText.setText("Click on NFC Fonction")
        self.executeFunction.emit(gui1) #call gui1 function in main window thread
        
        dataset = []
        _x, _y, _z, dump = self.robot.GetPosition().split(", ", 3)
        x = float(_x)
        y = float(_y)
        z = float(_z)
        dataset = nfc(x, y, z)
        self.robot.SetSpeed(speed)
        Acquire_points(dataset, self.robot, x, y, z, log_dir)
        def gui2():
            self.canMove = True
            self.buttonSetupOscilloscope.setEnabled(True)
            self.udpdateEnable()
        self.executeFunction.emit(gui2) #call gui2 function in main window thread
            
    """
     * @brief Launch the EMVCO Acquisition 
    """   
    @asynchrone
    def emvco(self):
        current_time = datetime.datetime.now()
        log_dir = f"logAcquisition/emvco_{current_time.strftime('%Y-%m-%d_%H-%M-%S')}"
        
        #display the action on the specific Text Box
        def gui1():
            self.canMove = False
            self.buttonSetupOscilloscope.setEnabled(False)
            self.udpdateEnable()
            validationText.setText("Click on EMVCO Fonction")
        self.executeFunction.emit(gui1) #call gui1 function in main window thread
        
        dataset = []
        _x, _y, _z, dump = self.robot.GetPosition().split(", ", 3)
        x = float(_x)
        y = float(_y)
        z = float(_z)
        dataset = emvco(x, y, z)
        self.robot.SetSpeed(speed)
        Acquire_points(dataset, self.robot, x, y, z, log_dir)
        
        def gui2():
            self.canMove = True
            self.buttonSetupOscilloscope.setEnabled(True)
            self.udpdateEnable()
        self.executeFunction.emit(gui2) #call gui2 function in main window thread
    
    """
     * @brief launch the custom cube acquisition
    """
    @asynchrone
    def customCube(self):
        ready = False
        current_time = datetime.datetime.now()
        log_dir = f"logAcquisition/cutomcube_{current_time.strftime('%Y-%m-%d_%H-%M-%S')}"
        dataset = []
        x = self.x_ptr
        y = self.y_ptr
        z = self.z_ptr
        
        def gui1():
            #update gui and get fields
            nonlocal x
            nonlocal y
            nonlocal z
            nonlocal ready
            nonlocal dataset
            self.canMove = False
            self.buttonSetupOscilloscope.setEnabled(False)
            self.udpdateEnable()
            validationText.setText("Click on customCube Fonction")
            
            
            _x, _y, _z, dump = self.robot.GetPosition().split(", ", 3)
            x = float(_x)
            y = float(_y)
            z = float(_z)
            dataset = customCube(x, y, z, 
                float(self.x_cube_size.text()), float(self.y_cube_size.text()), float(self.z_cube_size.text()), 
                int(self.x_cube_points.text()), int(self.y_cube_points.text()), int(self.z_cube_points.text()))
            ready = True
            
            
        self.executeFunction.emit(gui1) #call gui1 function in main window thread*
        #wait for gui1 to finish
        while not ready:
            time.sleep(0.001)
        self.robot.SetSpeed(speed)
        
        Acquire_points(dataset, self.robot, x, y, z, log_dir)
        
        def gui2():
            self.canMove = True
            self.buttonSetupOscilloscope.setEnabled(True)
            self.udpdateEnable()
        self.executeFunction.emit(gui2) #call gui2 function in main window thread
    
    
    """
     * @brief launch the custom cylindre acquisition
    """
    @asynchrone 

        #multithread creation for customcylindre,because the cylindre measure will take long time. the main window cannot wait always until the measurement finish without any operation
        # (if operate during the measure, no response from main window). after adding multithread, you can operate main window at the same time when measure is doing.
        
    
    def customCylindre(self):
        ready = False
        current_time = datetime.datetime.now()
        log_dir = f"logAcquisition/cutomcylindre_{current_time.strftime('%Y-%m-%d_%H-%M-%S')}"
        dataset = []
        x = self.x_ptr
        y = self.y_ptr
        z = self.z_ptr
        
        def gui1():
            #update gui and get fields
            nonlocal x
            nonlocal y
            nonlocal z
            nonlocal ready
            nonlocal dataset
            self.canMove = False
            self.buttonSetupOscilloscope.setEnabled(False)
            self.udpdateEnable()
            validationText.setText("Click on customCylindre Function")
            
            
            _x, _y, _z, dump = self.robot.GetPosition().split(", ", 3)
            x = float(_x)
            y = float(_y)
            z = float(_z)
            dataset = customCylindre(x, y, z, 
                float(self.field_radius.text()), float(self.field_height.text()), 
                int(self.field_z_points.text()), int(self.field_theta_points.text()))
            ready = True
            
            
        self.executeFunction.emit(gui1)
        while not ready:
            time.sleep(0.001)
        self.robot.SetSpeed(speed)
        Acquire_points(dataset, self.robot, x, y, z, log_dir)
        def gui2():
            self.canMove = True
            self.buttonSetupOscilloscope.setEnabled(True)
            self.udpdateEnable()
        self.executeFunction.emit(gui2)
    
    @asynchrone#multithread creation for customcylindre
    def customSemisphere(self):
        ready = False
        current_time = datetime.datetime.now()
        log_dir = f"logAcquisition/cutomSemisphere_{current_time.strftime('%Y-%m-%d_%H-%M-%S')}"
        dataset = []
        x = self.x_ptr
        y = self.y_ptr
        z = self.z_ptr
        
        def gui1():
            nonlocal x
            nonlocal y
            nonlocal z
            nonlocal ready
            nonlocal dataset
            self.canMove = False
            self.buttonSetupOscilloscope.setEnabled(False)
            self.udpdateEnable()
            validationText.setText("Click on customSemisphere Function")
            
            
            _x, _y, _z, dump = self.robot.GetPosition().split(", ", 3)
            x = float(_x)
            y = float(_y)
            z = float(_z)
            dataset = customSemisphere(x, y, z, 
                float(self.field_radius.text()), 
                int(self.field_z_points.text()), int(self.field_theta_points.text()))
            ready = True
            
            
        self.executeFunction.emit(gui1)
        while not ready:
            time.sleep(0.001)
        self.robot.SetSpeed(speed)
        Acquire_points(dataset, self.robot, x, y, z, log_dir)
        def gui2():
            self.canMove = True
            self.buttonSetupOscilloscope.setEnabled(True)
            self.udpdateEnable()
        self.executeFunction.emit(gui2)

    """
     * @brief Launch the measureOscillo function of the oscilloscope
    """
    @asynchrone
    def MeasureOnePoint(self):
        _x, _y, _z, dump = self.robot.GetPosition().split(", ", 3)
        x = float(_x)
        y = float(_y)
        z = float(_z)
        if self.robot.type == "DENSO":
            err = self.robot.Energize(0)
        getAcquisition(f"({x-self.fixed_x} {y-self.fixed_y} {z-self.fixed_z})", 0, self.point_log_dir)
        if self.robot.type == "DENSO":
            err = self.robot.Energize(1)
        def gui():
            validationText.setText("Point done")
        self.executeFunction.emit(gui)
            
    """
     * @brief Disconnect the Robot Communication 
    """ 
    def disconnect(self): 
        self.robot.CloseComm()
        self.canMove = False
        validationText.setText("Robot Disconnection Done")
        self.buttonDisconnectRobot.setEnabled(False)
        self.buttonConnectRobot.setEnabled(True)
        if self.robot.type == "DENSO":
            self.buttonRobot5Axes.setEnabled(True)
            self.buttonRobot6Axes.setEnabled(False)
        else:
            self.buttonRobot5Axes.setEnabled(False)
            self.buttonRobot6Axes.setEnabled(True)
        self.udpdateEnable()
        
    """
     * @brief Check if the text entered in the text field is a number or not.
    """ 
    def validateCoordPt(self):       
        if not self.x.text() or not self.y.text() or not self.z.text() :   #The result of the .text() method is a string containing the combined text of all matched elements. 
                                                                          #self.x is Qlineedit class, self.x.text() is its value in string if it is not empty. if it's empty, return false. If certen value, return is True.
                validationText.setText("At least one of the coordinates is empty")
        else :
            if is_float(self.x.text()) and is_float(self.y.text()) and is_float(self.z.text()) :     #The isnumeric() method checks if all the characters in the string are numeric.
                self.x_ptr = float(self.x.text())    # convert the string value of self.x.text() to int
                self.y_ptr = float(self.y.text())
                self.z_ptr = float(self.z.text())
                validationText.setText("Coordinates Validate")
            else : validationText.setText("x, y or z is not a nomber")
    
    """
     * @brief Notify a change on the speed slider
    """ 
    def valueChanged(self):
        global speed
        slider = self.sender()              # self.sender() in a function connected to your button event to get the object that triggered the event.
                
        speed = slider.value()
        validationText.setText("speed : " + str(speed))    #validationText is defined in line 206
        print("speed : ", speed)  # To show the value in console
          
    """
     * @brief Notify a change on the COM list
    """   
    def valuePort(self):
        global COM
        port = self.sender()     # self.sender() in a function connected to your button event to get the object that triggered the event.
        
        print("port : ", port.currentIndex())  #  currentIndex() is the number the sender() get from port list (the line you chose with your mouse). 0-->first line or 1--> second line etc.
        COM = self.portlist[port.currentIndex()][0]  #if the first line chosen, currentIndex() = 0, so portlist[port.currentIndex][0] =  portlist[0][0]. Set FIRST LINE, first colone of portlist to com
        print("COM: ", COM)    # To show the value in console


    """
        * @brief Deenergize the robot and write the information in the console
    """
    def EmergencyStop(self):
        emergencyStop[0]=True
        self.robot.Energize(0)
        validationText.setText("Robot Stopped")
    
    def figureLaunch(self):
        self.getFigureList()[1][self.figureSelector.currentIndex()]()
    
    def getListStep(self):
        return [["10 mm", "5 mm", "4 mm", "3 mm", "2 mm", "1 mm","0,5 mm"],
                [10, 5, 4, 3, 2, 1, 0.5]]
    
    """
     * @brief update the step for the buttons to move robot
    """
    def changeStep(self):
        index = self.stepSelector.currentIndex()
        self.distanceMove = self.getListStep()[1][index]
    
    """
     * @brief return a list containing 3 lists for the figures [0:name of the figure, 1:function of the figure, 2:image of the figure
    """
    def getFigureList(self):
        return [
            ['Manual Measure', 'NFC', 'EMVCO', 'Custom Cube', 'One point', 'Custom Cylinder','Semi-sphere'],
            [self.MeasureOscillo, self.nfc, self.emvco, self.customCube, self.MeasureOnePoint, self.customCylindre, self.customSemisphere],
            ["image: url(./measure.png)", "image: url(./nfc.jpg)", "image: url(./emvco.jpg)", "image: url(./cube.png)", "image: url(./point.png)", "image: url(./point.png)", "image: url(./point.png)"]
        ]
    
    """
     * @brief update gui when a figure is selected
    """
    def onFigureSelected(self):
        index = self.figureSelector.currentIndex()
        customCubeFieldWidgets = [self.label_x_cube_size, self.label_y_cube_size, self.label_z_cube_size, self.x_cube_size, self.y_cube_size, self.z_cube_size, self.label_x_cube_points, self.label_y_cube_points, self.label_z_cube_points, self.x_cube_points, self.y_cube_points, self.z_cube_points, self.unit_label]
        customCylindreFieldWidgets = [self.label_radius,self.field_radius, self.label_height, self.field_height, self.label_z_points, self.field_z_points, self.label_theta_points, self.field_theta_points, self.unit_label]
        customSemispereFieldWidgets = [self.label_radius,self.field_radius, self.label_z_points, self.field_z_points, self.label_theta_points, self.field_theta_points, self.unit_label]
        
        if self.getFigureList()[1][self.figureSelector.currentIndex()] == self.customCube:
            for w in customCubeFieldWidgets:
                w.show()
            self.setCubeLabels()
        elif self.getFigureList()[1][self.figureSelector.currentIndex()] == self.customCylindre:
            for w in customCubeFieldWidgets:
                w.hide()
            for w in customCylindreFieldWidgets:
                w.show()
            self.setCylindreLabels()
        elif self.getFigureList()[1][self.figureSelector.currentIndex()] == self.customSemisphere:
            for w in customCubeFieldWidgets:
                w.hide()
            for w in customSemispereFieldWidgets:
                w.show()
            self.setSemisphereLabels()
        else:
            for w in customCubeFieldWidgets:
                w.hide()
        if self.getFigureList()[1][self.figureSelector.currentIndex()] == self.MeasureOnePoint:
            for w in [self.buttonFixPoint]:
                w.show()
        else:
            for w in [self.buttonFixPoint]:
                w.hide()
                
        
        self.left_widget.setStyleSheet(self.getFigureList()[2][self.figureSelector.currentIndex()])
        self.udpdateEnable()
    
    """
     * @brief update enable state of the buttons
    """
    def udpdateEnable(self):
        if self.canMeasure:
            if self.getFigureList()[1][self.figureSelector.currentIndex()] == self.MeasureOscillo:
                self.buttonLaunch.setEnabled(True)
        else:
            self.buttonLaunch.setEnabled(False)
            
        if self.canMove : #and self.RobID != "5axes"
            self.xPlus.setEnabled(True)
            self.yPlus.setEnabled(True)
            self.zPlus.setEnabled(True)
            self.xMinus.setEnabled(True)
            self.yMinus.setEnabled(True)
            self.zMinus.setEnabled(True)
            self.buttonSavePoint.setEnabled(True)
            self.buttonOpenPoint.setEnabled(True)
            self.buttonShowCoords.setEnabled(True)
            self.buttonHomePosition.setEnabled(True)
        elif self.canMove: #canMove and 6 axes
            self.buttonHomePosition.setEnabled(True)
        else: #not canMove
            self.xPlus.setEnabled(False)
            self.yPlus.setEnabled(False)
            self.zPlus.setEnabled(False)
            self.xMinus.setEnabled(False)
            self.yMinus.setEnabled(False)
            self.zMinus.setEnabled(False)
            self.buttonShowCoords.setEnabled(False)
            self.buttonHomePosition.setEnabled(False)
            self.buttonSavePoint.setEnabled(False)
            self.buttonOpenPoint.setEnabled(False)
            
        if self.canMeasure and self.canMove and self.getFigureList()[1][self.figureSelector.currentIndex()] != self.MeasureOscillo: #canMove and canMeasure except for manualMeasure
            if self.getFigureList()[1][self.figureSelector.currentIndex()] == self.customCube: 
                fields = [self.x_cube_size, self.y_cube_size, self.z_cube_size]
                fields_int = [self.x_cube_points, self.y_cube_points, self.z_cube_points]
                self.buttonLaunch.setEnabled(False)
                for field in fields:
                    if not field.text() or not is_positive_float(field.text()):
                        break
                else:
                    for f in fields_int:
                        if not f.text() or not f.text().isnumeric():
                            break
                    else:
                        self.buttonLaunch.setEnabled(True)
            elif self.getFigureList()[1][self.figureSelector.currentIndex()] == self.customCylindre:
                fields = [self.field_height, self.field_radius]
                fields_int = [self.field_theta_points, self.field_z_points]
                self.buttonLaunch.setEnabled(False)
                for field in fields:
                    if not field.text() or not is_positive_float(field.text()):
                        break
                else:
                    for f in fields_int:
                        if not f.text() or not f.text().isnumeric():
                            break
                    else:
                        self.buttonLaunch.setEnabled(True)
            elif self.getFigureList()[1][self.figureSelector.currentIndex()] == self.customSemisphere:
                fields = [self.field_radius]
                fields_int = [self.field_theta_points, self.field_z_points]
                self.buttonLaunch.setEnabled(False)
                for field in fields:
                    if not field.text() or not is_positive_float(field.text()):
                        break
                else:
                    for f in fields_int:
                        if not f.text() or not f.text().isnumeric():
                            break
                    else:
                        self.buttonLaunch.setEnabled(True)
            else:
                self.buttonLaunch.setEnabled(True)
        elif self.getFigureList()[1][self.figureSelector.currentIndex()] != self.MeasureOscillo:
            self.buttonLaunch.setEnabled(False)
        if self.canMeasure and self.canMove and self.getFigureList()[1][self.figureSelector.currentIndex()] == self.MeasureOnePoint:
            self.buttonFixPoint.setEnabled(True)
        else:
            self.buttonFixPoint.setEnabled(False)
    
    @asynchrone
    def moveXplus(self):
        def gui1():
            self.canMove = False
            self.udpdateEnable()
        self.executeFunction.emit(gui1)
        
        self.robot.SetSpeed(speed)
        self.robot.Move(self.distanceMove, 0, 0, 5000)
        
        def gui2():
            self.canMove = True
            self.udpdateEnable()
            self.showCoords()
        self.executeFunction.emit(gui2)
    
    @asynchrone
    def moveYplus(self):
        def gui1():
            self.canMove = False
            self.udpdateEnable()
        self.executeFunction.emit(gui1)
        
        self.robot.SetSpeed(speed)
        self.robot.Move(0, self.distanceMove, 0, 5000)
        
        def gui2():
            self.canMove = True
            self.udpdateEnable()
            self.showCoords()
        self.executeFunction.emit(gui2)
    
    @asynchrone    
    def moveZplus(self):
        def gui1():
            self.canMove = False
            self.udpdateEnable()
        self.executeFunction.emit(gui1)
        
        self.robot.SetSpeed(speed)
        self.robot.Move(0, 0, self.distanceMove, 5000)
        
        def gui2():
            self.canMove = True
            self.udpdateEnable()
            self.showCoords()
        self.executeFunction.emit(gui2)
    
    @asynchrone    
    def moveXminus(self):
        def gui1():
            self.canMove = False
            self.udpdateEnable()
        self.executeFunction.emit(gui1)
        
        self.robot.SetSpeed(speed)
        self.robot.Move(-self.distanceMove, 0, 0, 5000)
        
        def gui2():
            self.canMove = True
            self.udpdateEnable()
            self.showCoords()
        self.executeFunction.emit(gui2)
    
    @asynchrone
    def moveYminus(self):
        def gui1():
            self.canMove = False
            self.udpdateEnable()
        self.executeFunction.emit(gui1)
        
        self.robot.SetSpeed(speed)
        self.robot.Move(0, -self.distanceMove, 0, 5000)
        
        def gui2():
            self.canMove = True
            self.udpdateEnable()
            self.showCoords()
        self.executeFunction.emit(gui2)
        
    @asynchrone
    def moveZminus(self):
        def gui1():
            self.canMove = False
            self.udpdateEnable()
        self.executeFunction.emit(gui1)
        
        self.robot.SetSpeed(speed)
        self.robot.Move(0, 0, -self.distanceMove, 5000)
        
        def gui2():
            self.canMove = True
            self.udpdateEnable()
            self.showCoords()
        self.executeFunction.emit(gui2)
        
    @asynchrone
    def moveHome(self):
        def gui1():
            self.canMove = False
            self.udpdateEnable()
        self.executeFunction.emit(gui1)
        
        self.robot.SetSpeed(speed)
        if self.RobID == "DENSO":
            self.robot.GoToHome(6000)
        else:
            self.robot.GoToPosition( 0, 269, 27, 0,90, 6000)
        
        def gui2():
            self.canMove = True
            self.udpdateEnable()
            self.showCoords()
        self.executeFunction.emit(gui2)
        
    """
     * @brief affiche les coordonnées du robot 
    """
    def showCoords(self):
        pos = self.robot.GetPosition()
        print(pos)
        x, y, z, dump = pos.split(", ", 3)
        txtCoords = "coords : (x: {}, y: {}, z: {})".format(x, y, z)
        def gui():
            validationText.setText(txtCoords)
        self.executeFunction.emit(gui)
        
        print(txtCoords)
    
    """
     * @brief défini le point initial de la figure one point en tant que la position actuelle du robot. Permet aussi de définir l'emplacement des prochaines mesures de cette figure
    """
    def fixPoint(self):
        x, y, z, dump = self.robot.GetPosition().split(", ", 3)
        txtCoords = "point fixed : (x: {}, y: {}, z: {})".format(x, y, z)
        validationText.setText(txtCoords)
        print(txtCoords)
        self.fixed_x = float(x)
        self.fixed_y = float(y)
        self.fixed_z = float(z)
        current_time = datetime.datetime.now()
        self.point_log_dir = f"logAcquisition/points_{current_time.strftime('%Y-%m-%d_%H-%M-%S')}"
    
    """
     * @brief permet à l'utilisateur d'enregistrer la position du robot dans un fichier
    """
    def savePoint(self):
        _x, _y, _z, dump = self.robot.GetPosition().split(", ", 3)
        x = float(_x)
        y = float(_y)
        z = float(_z)
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self, 
            "Save File", "", "All Files(*);;Text Files(*.txt)", options = options)
        if fileName:
            with open(fileName, 'w') as f:
                f.write(f"{x},{y},{z}")
    
    
    """
     * @brief permet à l'utilisateur d'ouvrir un fichier contenant un point puis positionne le robot sur ce point
    """
    def openPoint(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, 
            "Open File", "", "All Files(*);;Text Files(*.txt)", options = options)
        if fileName:
            with open(fileName, 'r') as f:
                xyz = f.read()
                _x, _y, _z = xyz.split(",")
                x = float(_x)
                y = float(_y)
                z = float(_z)
                self.robot.SetSpeed(speed)
                self.canMove = False
                self.udpdateEnable()
                if self.robot.type == "DENSO":
                    self.robot.GoToPosition6x(x,y,z, 180,0, 180, 10000)
                else:
                    self.robot.GoToPosition(x, y, z, 0, 90, 6000)
                self.canMove = True
                self.udpdateEnable()
                self.showCoords()
                print(f"({x}, {y}, {z})")
    
    """
     * @brief permet d'ouvrir la popup pour choisir les mesures puis les envoient à l'oscilloscope
    """
    def setupMeasure(self):
        self.buttonSetupOscilloscope.setEnabled(False)
        self.centralWidget().setEnabled(False)
        self.popup = MeasureSetupPopup()
        self.popup.show()
        self.popup.onClose.connect(lambda: self.buttonSetupOscilloscope.setEnabled(True))
        def sendConfig(config):
            setOscilloscopeParameters(config)
            self.centralWidget().setEnabled(True)
            self.canMeasure = True
            self.udpdateEnable()
        def onCancel():
            self.centralWidget().setEnabled(True)
            self.canMeasure = True
            self.udpdateEnable()
        self.popup.onValidate.connect(sendConfig)
        self.popup.onClose.connect(onCancel)
    
    """
     * @brief change le texte des label pour la figure du cube
    """
    def setCubeLabels(self):
        self.label_x_cube_points.setText("x points :")
        self.label_y_cube_points.setText("y points :")
        self.label_z_cube_points.setText("z points :")
        
        self.label_x_cube_size.setText("x size :")
        self.label_y_cube_size.setText("y size :")
        self.label_z_cube_size.setText("z size :")
    
    """
     * @brief change le texte des label pour la figure du cylindre
    """
    def setCylindreLabels(self):
        self.label_height.setText("Height :")
        self.label_radius.setText("Radius :")
        self.label_theta_points.setText("Circle points :")
        self.label_z_points.setText("Layers :")

    """
     * @brief change le texte des label pour la figure du Semi-sphere
    """
    def setSemisphereLabels(self):
        self.label_radius.setText("Radius :")
        self.label_theta_points.setText("Circle points :")
        self.label_z_points.setText("Layers :")


        
if __name__ == '__main__':
    app = QApplication(sys.argv)  #
    window = MainWindow()
    window.show()
    sys.exit(app.exec_()) # for QT window need to begin with  QApplication(sys.argv) and finish with sys.exit(app.exec_())
    