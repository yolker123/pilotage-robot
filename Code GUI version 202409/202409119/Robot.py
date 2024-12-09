
"""
 * @file        Robot.py
 * @brief       Contains the functions to control the robot and get its position and limits
 * @author      Lisa Duterte | Romain Derrien | Clement Rouvier | Elsa Della Valle 
 * @version     0.1
 * @date        2023
"""

"""
Code based on Common.py example code from keolabs
code is availabe in C:/Program Files (x86)/KEOLABS/RobX/Examples/Python/Common.py
"""

import win32com.client
import ctypes
import winreg
import threading

"""
 * @brief RobotObject class
 * @param type: type of the robot
 * @param RobX: object to control the robot

"""

MessageBox = ctypes.windll.user32.MessageBoxA
MB_OK = 0
MB_SYSTEMMODAL = 4096
MB_ICONERROR = 16
MB_ICONWARNING = 48

cart_mode = 1
scale = 1+cart_mode*9 # 10 if cart_mode = 1, 1 if cart_mode = 0



#debug traces
#0 = no
#1 = yes
debug = 0
mutex = threading.Lock()
th_lock = None
global workers
workers = {}
global emergencyStop
emergencyStop=[False]  #variable globale permettant d'arreter le robot en cas d'urgence


"""
 * @brief annonation permettant à chaque threads d'obtenir une communication avec le programme KEOLABS.
"""
def reDispatch(fn):
    def wrapper(self, *args, **kwargs):
        result = None
        has_lock = (th_lock != threading.current_thread())
        if has_lock:
            mutex.acquire()
        if workers.get(threading.current_thread().ident) == None:
            workers[threading.current_thread().ident] = [win32com.client.Dispatch("KEOLABS.RobX"), threading.current_thread()]
        else:
            if workers.get(threading.current_thread().ident)[0] is None:
                workers[threading.current_thread().ident][0] = win32com.client.Dispatch("KEOLABS.RobX")
            self.RobX = workers[threading.current_thread().ident][0]
        
        try:
            result =  fn(self, *args, **kwargs)
        except Exception as e:
            print("error : ", e)
        if has_lock:
            mutex.release()
        return result
    return wrapper

class RobotObject:
    def __init__(self, type):
        self.type = type
        self.RobX = win32com.client.Dispatch("KEOLABS.RobX")
        

    def LogError(msg):
        MessageBox(0,str(msg),"Error", MB_OK | MB_ICONERROR | MB_SYSTEMMODAL)

    def LogDebug(msg):
        if (debug == 1) :
            MessageBox(0,str(msg),"Debug", MB_OK | MB_ICONWARNING | MB_SYSTEMMODAL)
    
    @reDispatch
    def OpenComm(self, COM=None):
        err=0
        if(self.type == "DENSO"):
            err = self.RobX.RobotSelect(6, "DENSO", "VS-AV6")
            err = self.RobX.OpenComm("192.168.30.34")
        elif(self.type == "5axes"):
            self.RobX.RobotSelect(5)
            print(COM)
            err = self.RobX.OpenComm(COM)
        if err != 0:
            print(f"Error opening communication : {err}")
            return err
        
        #Cartesian mode 0 -> unit in mm
        #Cartesian mode 1 -> unit in 1/10 mm
        if(self.type == "DENSO"):
            err = self.RobX.SetCartesianExtMode(cart_mode)
            err = self.RobX.InitRobot()
        elif (self.type == "5axes") :
            err = self.RobX.InitRobot()
            if err != 0: return err
            #Making sure Robot is correctly placed
            err = self.RobX.Energize(0)
            if err != 0: return err
            MessageBox(0,"Make sure Robot is approximately in Home Position","User action required", MB_OK | MB_SYSTEMMODAL) 
            err = self.RobX.Energize(1)
            if err != 0: return err
            err = self.RobX.SetCartesianExtMode(cart_mode)
            if err != 0: return err
            try:
                err = self.RobX.SetToolParam(220*scale)
                if err != 0: return err
                #err = self.RobX.Align(1)
                if err != 0: return err
                err = self.RobX.Calibrate(1)
                if err != 0: return err
            except Exception as detail:
                print(detail)
            err = self.RobX.SetSpeed(2000)
            if err != 0: return err
        if err != 0:
            print(f"Error initializing robot : {err}")
            return err
        return 0
    
    @reDispatch
    def Energize(self, state):
        if emergencyStop[0]:
            state = 0
        err = self.RobX.Energize(state)
        return err
    
    @reDispatch
    def GoToPosition(self, x, y, z, w, p, timeout):
        err = 0
        if not emergencyStop[0]:
            err = self.RobX.GoToPosition(x*scale, y*scale, z*scale, w*scale, p*scale, timeout)
        return err
    
    @reDispatch
    def GoToPosition6x(self, x, y, z, rX, rY, rZ, timeout):
        err = 0
        if not emergencyStop[0]:
            err = self.RobX.RestoreDefaultLimits()
            err = self.RobX.GoToPosition6x(x*scale, y*scale, z*scale, rX*scale, rY*scale, rZ*scale, timeout)
        return err
    
    @reDispatch
    def GoToHome(self, timeout):
        err = 0
        if not emergencyStop[0]:
            err = self.RobX.GoToHome(timeout)
            #if self.type == "5axes":
                #err = self.GoToPosition(0.0, 398.0, 58.0, 0, 90, 6000)
                
        return err
    @reDispatch
    def CloseComm(self):
        self.RobX.Energize(0)
        self.RobX.CloseComm()
        return 0
    
    @reDispatch
    def SetSpeed(self, speed, rob=None):
        err = 0
        err = self.RobX.SetSpeed(speed)
        return err
    
    @reDispatch
    def GetLimits(self):
        err = self.RobX.GetLimitsStatus()
        if err != 0:
            return err
        else:
            return self.RobX.GetLimitsStr()
    
    @reDispatch
    def GetPosition(self):
        err = self.RobX.GetCurrentPositionStatus()
        if err != 0:
            #if self.type == "5axes":
            #    return "0.0, 398.0, 58.0, 0, 90"
            #else:
            return ""
        else:
            _x,_y,_z, _w, _p = self.RobX.GetCurrentPositionStr().split(", ", 5)
            x=int(_x)/scale
            y=int(_y)/scale
            z=int(_z)/scale
            w=int(_w)/scale
            p=int(_p)/scale
            scale
            return "{}, {}, {}, {}, {}".format(x,y,z,w,p)
    
    @reDispatch
    def Move(self, x, y, z, timeout):
        err = 0
        if not emergencyStop[0]:
            err = self.RobX.Move(x*scale, y*scale, z*scale, timeout)
        if err != 0:
            print(f"Error moving : {err}")
        return err




