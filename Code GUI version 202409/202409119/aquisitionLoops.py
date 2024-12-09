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

timeToSleep = 2 #time before getting the oscilloscope acquisition when the robot move




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
    def Acquire_point(self, robot, log_dir="logAcquisition"):
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
        getAcquisition(f"{self.name}", 0, log_dir)
        print("aquired")
    

def createRobot(type):
    robot = RobotObject(type)
    return robot

"""
    * @brief Launch the acquisition of data on each points of a given dataset
    * @param points: list of points to acquire
"""
def Acquire_points(points, robot, x_ptr = None, y_ptr = None, z_ptr = None, log_dir="logAcquisition"):
    for point in points:
        point.Acquire_point(robot, log_dir)
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
    robot.GoToPosition(0, 386, 161, -77, 57, timeout) #sort du champs d'alimentation
    
"""
    * @brief Do the acquisition loop on the default dataset NFC
    * @param x_ptr, y_ptr, z_ptr: position of the point (0,0,0) of the dataset
    * @return the mesh
"""
def nfc(x_ptr, y_ptr, z_ptr):
    h_nfc = 5    #hauteur  unit: mm
    rb_nfc = 5   #rayon bas
    rh_nfc = 10  #rayon haut
    w=0
    p=90
    timeout=6000
    
    print('NFC Start')

    dataset = []
    for index, dx in enumerate([-rb_nfc, 0, rb_nfc]):
        dataset.append(Point(f"({index-1} 0 0)",x_ptr+dx, y_ptr, z_ptr, w, p, z_ptr, timeout))
    for index, dy in enumerate([-rb_nfc, 0, rb_nfc]):
        dataset.append(Point(f"(0 {index-1} 0)",x_ptr, y_ptr+dy, z_ptr, w, p, z_ptr, timeout))
    for index, dx in enumerate([-rh_nfc, 0, rh_nfc]):
        dataset.append(Point(f"({index-1} 0 1)",x_ptr+dx, y_ptr, z_ptr+h_nfc, w, p, z_ptr, timeout))
    for index, dy in enumerate([-rh_nfc, 0, rh_nfc]):
        dataset.append(Point(f"(0 {index-1} 1)",x_ptr, y_ptr+dy, z_ptr+h_nfc, w, p, z_ptr, timeout))
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
        dataset.append(Point(f"({index-1} 0 0)",x_ptr+dx, y_ptr, z_ptr, 0, 90, z_ptr, Temps_limit))
    for index, dy in enumerate([-rp_emvco, 0, rp_emvco]):
        dataset.append(Point(f"(0 {index-1} 0)",x_ptr, y_ptr+dy, z_ptr, 0, 90, z_ptr, Temps_limit))
    
    for index, dx in enumerate([-rg_emvco, 0, rg_emvco]):
        dataset.append(Point(f"({index-1} 0 1)",x_ptr+dx, y_ptr, z_ptr+h_emvco, 0, 90, z_ptr, Temps_limit))
    for index, dy in enumerate([-rg_emvco, 0, rg_emvco]):
        dataset.append(Point(f"(0 {index-1} 1)",x_ptr, y_ptr+dy, z_ptr+h_emvco, 0, 90, z_ptr, Temps_limit))

    for index, dx in enumerate([-rg_emvco, 0, rg_emvco]):
        dataset.append(Point(f"({index-1} 0 2)",x_ptr+dx, y_ptr, z_ptr+2*h_emvco, 0, 90, z_ptr, Temps_limit))
    for index, dy in enumerate([-rg_emvco, 0, rg_emvco]):
        dataset.append(Point(f"(0 {index-1} 2)",x_ptr, y_ptr+dy, z_ptr+2*h_emvco, 0, 90, z_ptr, Temps_limit))

    for index, dx in enumerate([-rg_emvco, 0, rg_emvco]):
        dataset.append(Point(f"({index-1} 0 3)",x_ptr+dx, y_ptr, z_ptr+3*h_emvco, 0, 90, z_ptr, Temps_limit))
    for index, dy in enumerate([-rg_emvco, 0, rg_emvco]):
        dataset.append(Point(f"(0 {index-1} 3)",x_ptr, y_ptr+dy, z_ptr+3*h_emvco, 0, 90, z_ptr, Temps_limit))

    for index, dx in enumerate([-rp_emvco, 0, rp_emvco]):
        dataset.append(Point(f"({index-1} 0 4)",x_ptr+dx, y_ptr, z_ptr+4*h_emvco, 0, 90, z_ptr, Temps_limit))
    for index, dy in enumerate([-rp_emvco, 0, rp_emvco]):
        dataset.append(Point(f"(0 {index-1} 4)",x_ptr, y_ptr+dy, z_ptr+4*h_emvco, 0, 90, z_ptr, Temps_limit))
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
        dx = x_size/(x_points-1) #size step on x
    if y_points > 1:
        dy = y_size/(y_points-1) #size step on y 
    if z_points > 1:
        dz = z_size/(z_points-1) #size step on z
    dataset = []
    z = z_ptr #current z position
    for i in range(z_points):
        y = y_ptr - y_size/2 #current y position
        for j in range(y_points):
            x = x_ptr - x_size/2 #current x position
            for k in range(x_points):
                dataset.append(Point(f"({x-x_ptr} {y-y_ptr} {z-z_ptr})", x, y, z, 0, 90, z_ptr, time_limit))
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
    dz = 0 #distance between z points
    if h_points > 1:
        dz = h/(h_points-1)
    dataset = []
    dtheta = 0 #distance between circle points
    
    if circle_points > 1:
        dtheta = 2*math.pi/circle_points
    z = z_ptr
    for i in range(h_points):
        theta = 0
        dataset.append(Point(f"(0 0 {z-z_ptr})", x_ptr, y_ptr, z, 0, 90, z_ptr, time_limit))
        for j in range(circle_points):
            x = r*math.cos(theta)+x_ptr
            y = r*math.sin(theta)+y_ptr
            dataset.append(Point(f"({x-x_ptr} {y-y_ptr} {z-z_ptr})", x, y, z, 0, 90, z_ptr, time_limit))
            theta += dtheta
        z += dz
    return dataset


def customSemisphere(x_ptr, y_ptr, z_ptr, r, z_points, circle_points):
    import math
    time_limit = 6000
    dz = 0
    if z_points > 1:
        dz = r/(z_points-1)
    dataset = []
    dtheta = 0
    
    if circle_points > 1:
        dtheta = 2*math.pi/circle_points
    z = z_ptr
    
    for i in range(z_points):
        theta = 0
        dataset.append(Point(f"(0 0 {z-z_ptr})", x_ptr, y_ptr, z, 0, 90, z_ptr, time_limit))
        if i < z_points -1:
            h=i*dz  # h is the height in z direction between intial position and measure point# 
            l=math.sqrt(r*r-h*h)  # l is the radius of each layer
            
            for j in range(circle_points):     
                x = l*math.cos(theta)+x_ptr
                y = l*math.sin(theta)+y_ptr
                dataset.append(Point(f"({x-x_ptr} {y-y_ptr} {z-z_ptr})", x, y, z, 0, 90, z_ptr, time_limit))
                theta += dtheta
        z += dz
    return dataset