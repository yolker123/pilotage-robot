# ------------ Project-specific functions ------------
"""
 * @file        Common.py
 * @brief       Contains the functions launching the specifics desired Acquisitions 
 * @author      Lisa Duterte | Romain Derrien | Clement Rouvier | Elsa Della Valle | Romeo Botuli-Bundol | Hamid Ajouaou | Samuel Decay | Baptiste Saby
 * @version     0.1
 * @date        2023
"""

from Robot import *
from oscilloscopeAcquisition import *
import time
from tektronix import *

timeToSleep = 2  # time before getting the oscilloscope acquisition when the robot move

"""
    * @brief Point class
    * @param name: name of the point
    * @param x: x coordinate of the point
    * @param y: y coordinate of the point
    * @param z: z coordinate of the point
    * @param w: rotation around the x axis
    * @param p: rotation around the y axis
    * @param timeout: timeout for the movement
"""


class Point:
    def __init__(self, name, x, y, z, w, p, zMin, timeout):
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        self.w = w
        self.p = p
        self.zMin = zMin
        self.timeout = timeout

    def setMinZ(self, z):
        self.zMin = z

    """
        * @brief Launch the acquisition of data on a given point
    """

    def Acquire_point(self, robot, mfa, log_dir="logAcquisition", oscilloName="", tektronix=None):
        print(f"Acquiring {self.name}")
        if robot.type == "DENSO":
            err = robot.Energize(1)
            err = robot.GoToPosition6x(self.x, self.y, self.z, 180, 0, 180, self.timeout)
            err = robot.Energize(0)
        else:
            err = robot.GoToPosition(self.x, self.y, self.z, self.w, self.p, self.timeout)
        if err != 0:
            print(f"Error moving to {self.name} : {err}")
        time.sleep(timeToSleep)

        # TODO : MOCHE
        if oscilloName == "lecroy":
            getAcquisition(mfa, f"{self.name}", 0, log_dir)
        if oscilloName == "tektronix":
            values = tektronix.get_measures(f"{self.name}")
            point_ = self.name.strip("()")
            val = point_.split()
            x, y, z = map(float, val)
            print(values)
            Hx = None
            Hy = None
            Hz = None
            for value in values:
                if value['meas_type'] == "MAXIMUM" and value['channel'] == 'CH1':
                    print("value CH :", value)
                    Bx = float(value['value']) / (mfa.simulation.S * mfa.simulation.omega)
                    Hx = Bx / mfa.simulation.mu_0
                if value['meas_type'] == "MAXIMUM" and value['channel'] == 'CH2':
                    By = float(value['value']) / (mfa.simulation.S * mfa.simulation.omega)
                    Hy = By / mfa.simulation.mu_0
                if value['meas_type'] == "MAXIMUM" and value['channel'] == 'CH3':
                    Bz = float(value['value']) / (mfa.simulation.S * mfa.simulation.omega)
                    Hz = Bz / mfa.simulation.mu_0
            print("measure du point :", self.x, self.y, self.z)
            # -> calcule Hx Hy Hz
            if Hx is not None or Hy is not None or Hz is not None:
                if Hx is None:
                    Hx= 0
                if Hy is None:
                    Hy= 0
                if Hz is None:
                    Hz= 0
                mfa.simulation.measuredPoints.append({
                    'x': x, 'y': y, 'z': z,
                    'Hx': Hx, 'Hy': Hy, 'Hz': Hz, 'display': True
                })
                mfa.update_all_graphs(False)
        print("aquired")


def createRobot(type):
    robot = RobotObject(type)
    return robot


"""
    * @brief Launch the acquisition of data on each points of a given dataset
    * @param points: list of points to acquire
"""


def Acquire_points(points, robot, mfa, form, x_ptr=None, y_ptr=None, z_ptr=None, log_dir="logAcquisition", oscilloName="",
                   tektronix=None):
    z_ptr_reel = z_ptr - 211.1
    mfa.simulation.hRobot = z_ptr_reel
    if oscilloName == "tektronix":
        tektronix.create_file(z_ptr_reel, form)
        tektronix.scope.write("FPANEL:PRESS AUTOset")
    for point in points:
        point.Acquire_point(robot, mfa, log_dir, oscilloName, tektronix)
    mfa.update_all_graphs()
    if oscilloName == "tektronix":
        base_waveform_path = "C:/Users/Public/Tektronix/TekScope/WaveForm/"
        base_screenshot_path = "C:/Users/Public/Tektronix/TekScope/Screenshots/"
        destination_base_path = "//PCROBOT/Users/isen/PycharmProjects/pilotage-robot/Code GUI version 202409/202409119/Measure/TekScope/"

        # Construction des chemins complets pour le waveform
        waveform_source = base_waveform_path + f"logWaveform_tektronix_{tektronix.current_time_str}/"
        waveform_dest = destination_base_path + 'WaveForm/' + f"logWaveform_tektronix_{tektronix.current_time_str}/"
        # Construction des chemins complets pour les screenshots
        screenshot_source = base_screenshot_path + f"logScreenshot_tektronix_{tektronix.current_time_str}/"
        screenshot_dest = destination_base_path + 'Screenshots/' + f"logScreenshot_tektronix_{tektronix.current_time_str}/"
        # Copier les fichiers de waveform
        tektronix.scope.write(f'FILESystem:COPy "{waveform_source}", "{waveform_dest}"')
        # Copier les fichiers de screenshot
        tektronix.scope.write(f'FILESystem:COPy "{screenshot_source}", "{screenshot_dest}"')

    if x_ptr == None or y_ptr == None or z_ptr == None:
        return
    timeout = 6000
    if robot.type == "DENSO":
        robot.Energize(1)
        robot.GoToPosition6x(x_ptr, y_ptr, z_ptr, 180, 0, 180, timeout)
    else:
        robot.GoToPosition(x_ptr, y_ptr, z_ptr, 0, 90, timeout)


"""
    * @brief Bring the card out of the rf field in order to unload it
    * @param timeout: timeout for the movement
"""


def goHorsChamp(timeout, robot):
    robot.GoToPosition(0, 386, 161, -77, 57, timeout)  # sort du champs d'alimentation


"""
    * @brief Do the acquisition loop on the default dataset NFC
    * @param x_ptr, y_ptr, z_ptr: position of the point (0,0,0) of the dataset
    * @return the mesh
"""


def nfc(x_ptr, y_ptr, z_ptr):
    h_nfc = 5  # hauteur  unit: mm
    rb_nfc = 5  # rayon bas
    rh_nfc = 10  # rayon haut
    w = 0
    p = 90
    timeout = 6000

    print('NFC Start')

    dataset = []
    for index, dx in enumerate([-rb_nfc, 0, rb_nfc]):
        dataset.append(Point(f"({index - 1} 0 0)", x_ptr + dx, y_ptr, z_ptr, w, p, z_ptr, timeout))
    for index, dy in enumerate([-rb_nfc, 0, rb_nfc]):
        dataset.append(Point(f"(0 {index - 1} 0)", x_ptr, y_ptr + dy, z_ptr, w, p, z_ptr, timeout))
    for index, dx in enumerate([-rh_nfc, 0, rh_nfc]):
        dataset.append(Point(f"({index - 1} 0 1)", x_ptr + dx, y_ptr, z_ptr + h_nfc, w, p, z_ptr, timeout))
    for index, dy in enumerate([-rh_nfc, 0, rh_nfc]):
        dataset.append(Point(f"(0 {index - 1} 1)", x_ptr, y_ptr + dy, z_ptr + h_nfc, w, p, z_ptr, timeout))
    return dataset


"""
    * @brief EMVCO Acquisition
    * @param x_ptr, y_ptr, z_ptr: position of the point (0,0,0) of the dataset
    * @return the mesh
"""


def emvco(x_ptr, y_ptr, z_ptr):
    h_emvco = 10
    rp_emvco = 15
    rg_emvco = 25

    Temps_limit = 6000

    print('emvco start')

    dataset = []
    for index, dx in enumerate([-rp_emvco, 0, rp_emvco]):
        dataset.append(Point(f"({index - 1} 0 0)", x_ptr + dx, y_ptr, z_ptr, 0, 90, z_ptr, Temps_limit))
    for index, dy in enumerate([-rp_emvco, 0, rp_emvco]):
        dataset.append(Point(f"(0 {index - 1} 0)", x_ptr, y_ptr + dy, z_ptr, 0, 90, z_ptr, Temps_limit))

    for index, dx in enumerate([-rg_emvco, 0, rg_emvco]):
        dataset.append(Point(f"({index - 1} 0 1)", x_ptr + dx, y_ptr, z_ptr + h_emvco, 0, 90, z_ptr, Temps_limit))
    for index, dy in enumerate([-rg_emvco, 0, rg_emvco]):
        dataset.append(Point(f"(0 {index - 1} 1)", x_ptr, y_ptr + dy, z_ptr + h_emvco, 0, 90, z_ptr, Temps_limit))

    for index, dx in enumerate([-rg_emvco, 0, rg_emvco]):
        dataset.append(Point(f"({index - 1} 0 2)", x_ptr + dx, y_ptr, z_ptr + 2 * h_emvco, 0, 90, z_ptr, Temps_limit))
    for index, dy in enumerate([-rg_emvco, 0, rg_emvco]):
        dataset.append(Point(f"(0 {index - 1} 2)", x_ptr, y_ptr + dy, z_ptr + 2 * h_emvco, 0, 90, z_ptr, Temps_limit))

    for index, dx in enumerate([-rg_emvco, 0, rg_emvco]):
        dataset.append(Point(f"({index - 1} 0 3)", x_ptr + dx, y_ptr, z_ptr + 3 * h_emvco, 0, 90, z_ptr, Temps_limit))
    for index, dy in enumerate([-rg_emvco, 0, rg_emvco]):
        dataset.append(Point(f"(0 {index - 1} 3)", x_ptr, y_ptr + dy, z_ptr + 3 * h_emvco, 0, 90, z_ptr, Temps_limit))

    for index, dx in enumerate([-rp_emvco, 0, rp_emvco]):
        dataset.append(Point(f"({index - 1} 0 4)", x_ptr + dx, y_ptr, z_ptr + 4 * h_emvco, 0, 90, z_ptr, Temps_limit))
    for index, dy in enumerate([-rp_emvco, 0, rp_emvco]):
        dataset.append(Point(f"(0 {index - 1} 4)", x_ptr, y_ptr + dy, z_ptr + 4 * h_emvco, 0, 90, z_ptr, Temps_limit))
    return dataset


"""
    * @brief Acquires a custom cube
    * @param x_ptr, y_ptr, z_ptr : position of the point (0,0,0) of the dataset
    * @param x_size, y_size, z_size : size of the edges of the cube
    * @param x_points, y_points, z_points : number of points for each edge
    * @return the mesh
"""


def customCube(x_ptr, y_ptr, z_ptr, x_size, y_size, z_size, x_points, y_points, z_points):
    print("cube start")
    time_limit = 6000
    dx = 0
    dy = 0
    dz = 0
    if x_points > 1:
        dx = x_size / (x_points - 1)  # size step on x
    if y_points > 1:
        dy = y_size / (y_points - 1)  # size step on y
    if z_points > 1:
        dz = z_size / (z_points - 1)  # size step on z
    dataset = []
    z = z_ptr  # current z position
    for i in range(z_points):
        y = y_ptr - y_size / 2  # current y position
        for j in range(y_points):
            x = x_ptr - x_size / 2  # current x position
            for k in range(x_points):
                dataset.append(Point(f"({x - x_ptr} {y - y_ptr} {z - z_ptr})", x, y, z, 0, 90, z_ptr, time_limit))
                x += dx
            y += dy
        z += dz
    return dataset


"""
    * @brief Acquires a custom cylindre
    * @param x_ptr, y_ptr, z_ptr : position of the point (0,0,0) of the dataset
    * @param r : radius
    * @param h : height
    * @param h_points : number of points on height
    * @param circle_points : number of points on the circles
    * @param x_points, y_points, z_points : number of points for each edge
    * @return the mesh
"""


def customCylindre(x_ptr, y_ptr, z_ptr, r, h, h_points, circle_points):
    import math
    time_limit = 6000
    dz = 0  # distance between z points
    if h_points > 1:
        dz = h / (h_points - 1)
    dataset = []
    dtheta = 0  # distance between circle points

    if circle_points > 1:
        dtheta = 2 * math.pi / circle_points
    z = z_ptr
    for i in range(h_points):
        theta = 0
        dataset.append(Point(f"(0 0 {z - z_ptr})", x_ptr, y_ptr, z, 0, 90, z_ptr, time_limit))
        for j in range(circle_points):
            x = r * math.cos(theta) + x_ptr
            y = r * math.sin(theta) + y_ptr
            dataset.append(Point(f"({x - x_ptr} {y - y_ptr} {z - z_ptr})", x, y, z, 0, 90, z_ptr, time_limit))
            theta += dtheta
        z += dz
    return dataset


def customSemisphere(x_ptr, y_ptr, z_ptr, r, z_points, circle_points):
    import math
    time_limit = 6000
    dz = 0
    if z_points > 1:
        dz = r / (z_points - 1)
    dataset = []
    dtheta = 0

    if circle_points > 1:
        dtheta = 2 * math.pi / circle_points
    z = z_ptr

    for i in range(z_points):
        theta = 0
        dataset.append(Point(f"(0 0 {z - z_ptr})", x_ptr, y_ptr, z, 0, 90, z_ptr, time_limit))
        if i < z_points - 1:
            h = i * dz  # h is the height in z direction between intial position and measure point#
            l = math.sqrt(r * r - h * h)  # l is the radius of each layer

            for j in range(circle_points):
                x = l * math.cos(theta) + x_ptr
                y = l * math.sin(theta) + y_ptr
                dataset.append(Point(f"({x - x_ptr} {y - y_ptr} {z - z_ptr})", x, y, z, 0, 90, z_ptr, time_limit))
                theta += dtheta
        z += dz
    return dataset