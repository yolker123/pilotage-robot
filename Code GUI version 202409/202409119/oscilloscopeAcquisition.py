import win32com.client
import pyvisa
#import visa
import time
from lecroy import *
from bddSetupOscilloscope import *
import os
import datetime
import csv

timeout_commands = 20000

global measure_config
global wf_img_config

"""
* @brief Enable the connection to the oscilloscope
* @param The Id of the Oscilloscope
* @exception range_error If the oscilloscope is not connected 
"""
def oscilloscopeConnection(id_oscillo):
    global lecroy

    scope = win32com.client.Dispatch("LeCroy.ActiveDSOCtrl.1") #creates instances of the ActiveDSO Control
    scope.MakeConnection("USBTMC:" + id_oscillo)  #Connects to the oscilloscope
    rm = pyvisa.ResourceManager()   #pyvisa.ResourceManager()  get return value of connection, 
                                    #'USB0::0x05FF::0x1023::3561N16324::INSTR' is the series of oscilloscope. if it's in the return, oscilloscope connected.
    print(rm.list_resources())
    lecroy = lecroy(pyvisa_instr=rm.open_resource(id_oscillo, timeout=timeout_commands), dso=scope)
    lecroy.scope.write("TRMD NORM")
    return rm


"""
* @brief Set the parameters of the oscilloscope depending on the desired acquisition 
"""
def setOscilloscopeParameters(config = None) :
    (ver_offset,analog_channels,measurement_channels) = acquisition1()  # get dictionary value from acquisition1() in bddSetupOscilloscope
    if config is None:
        config = measure_config
    measurement_channels = formatMeasureConfig(config)
    lecroy.measurement_setup(measurement_channels,numberOfAnalogChannels=8)
    for i in range(1, 5):
        view=False
        for _, it in measurement_channels.items():
            if it[0] == f"C{i}":
                view = True
                break
        if wf_img_config.get(f"C{i}").get("img") or wf_img_config.get(f"C{i}").get("wf"):
            view = True
        lecroy.scope.write(f'VBS app.Acquisition.C{i}.View={view}')
  
    
"""
* @brief Get the configuration of the oscilloscope
"""
def getOscillocopeConfiguration():
    sampleRate = lecroy.getSampleRate()
    horizontal_scale = lecroy.get_horizontal_scale()
    C1VerticalScale = lecroy.getChannelVerticalScale("C1")
    C2VerticalScale = lecroy.getChannelVerticalScale("C2")

    print(sampleRate, horizontal_scale, C1VerticalScale, C2VerticalScale)
    # In this way we delete the measurement of the C3 channel
    
"""
* @brief Get six specifics variables mesured by the oscilloscope and record
          the measures and the screenshot of the window.
"""
def getAcquisition(nomPoint, err=0, acquisition_directory="logAcquisition", config = None):
    if config is None:
        config = measure_config
        
    if not os.path.exists(acquisition_directory):
        os.makedirs(acquisition_directory)
    current_time = datetime.datetime.now()
    
    
    #________MEASURE RECORD__________

    # Get the date information
    
    # Create a file thanks to the data time
    log_filename = f"logMeasureOscillo_{nomPoint}_{current_time.strftime('%Y-%m-%d_%H-%M-%S')}.txt"
    # Specify the full path
    log_directory = acquisition_directory+"/logMeasure"
    # Check if the screenshot folder exists, create it if it doesn't
    if not os.path.exists(log_directory):
        os.makedirs(log_directory)
    # Create the full path of the directory
    log_filepath = os.path.join(log_directory, log_filename)

    if err != 0: 
        print(f"Error moving to {nomPoint} : {err}")
        with open(log_filepath, 'a') as f:
            f.write(f"{nomPoint}, time:{current_time.strftime('%H-%M-%S')}\nError moving to {nomPoint} : {err}\n")
        return

    conf = formatMeasureConfig(config)
    instant = current_time.strftime('%Y-%m-%d_%H-%M-%S')
    if len(conf) > 0:
        captureMeasures(conf, current_time, instant, acquisition_directory, nomPoint, log_filepath)

    for kt, it in wf_img_config.items():
        if it.get("wf"):
            captureWF(acquisition_directory, instant, kt, nomPoint)
        if it.get("img"):
            captureScreen(acquisition_directory, instant, kt, nomPoint)

"""
 * @brief capture la waveform
"""
def captureWF(acquisition_directory, instant, chan, nomPoint):
    # Enregistrer la forme d'onde dans un fichier Excel
    waveform_filename = f"logWaveform_{chan}_{nomPoint}_{instant}"
    waveform_directory = acquisition_directory+"/logWaveform"
    if not os.path.exists(waveform_directory):
        os.makedirs(waveform_directory)
    waveform_filepath = os.path.join(waveform_directory, waveform_filename+".txt")

    # Enregistrement et tranfert de la Waveform
    wf_dir = "D:\\Waveforms"
    wf_name = waveform_filename
    wf_path = wf_dir+"\\"+wf_name
    lecroy.saveWaveform(wf_dir, chan, wf_name)
    wf_reelpath = wf_dir+"\\"+chan+wf_name+"00000.csv"
    lecroy.transferFile(wf_reelpath, waveform_filepath)

"""
 * @brief capture l'image
"""
def captureScreen(acquisition_directory, instant, chan, nomPoint):
    screenshot_filename = f"logScreenshot_{chan}_{nomPoint}_{instant}.png"
    # We can put the %S for the seconds to avoid rewrite the screenshot if we take two measures in the same minute
    # Specify the full path to the screenshot folder
    screenshot_directory = acquisition_directory+"/logScreenshot"
    # Check if the screenshot folder exists, create it if it doesn't
    if not os.path.exists(screenshot_directory):
        os.makedirs(screenshot_directory)
    # Create the full path of the screenshot file
    screenshot_filepath = os.path.join(screenshot_directory, screenshot_filename)
    # Capture and save the screenshot
    lecroy.get_screen_image(path_with_filename=screenshot_filepath)
    print(f"Screenshot captured and saved to {screenshot_filepath}")


"""
 * @brief prend les P 
"""
def captureMeasures(conf, current_time, instant, acquisition_directory, nomPoint, log_filepath):
    # Get value of measurement channels
    measures={}
    
    
    for kp, p in conf.items():
        measures[f"{p[0]}_{p[1]}"] = "{:.4f}".format(lecroy.getValueOnChannel(f'P{kp}', 'value')) if isinstance(lecroy.getValueOnChannel(f'P{kp}', 'value'), (float, int)) else "{}".format(lecroy.getValueOnChannel(f'P{kp}', 'value'))
  
    with open(log_filepath, 'a') as f:  # Utiliser 'a' pour ajouter au fichier existant
        txt = f"{nomPoint}, time:{current_time.strftime('%H-%M-%S')}\n"
        for km, m in measures.items():
            txt += f"{km}:{m}\n"
        f.write(txt)
    # Enregistrer les mesures dans un fichier CSV
    csv_filename = f"logMeasure_{nomPoint}_{instant}.csv"
    csv_directory = acquisition_directory+"/logCSV"
    if not os.path.exists(csv_directory):
        os.makedirs(csv_directory)
    csv_filepath = os.path.join(csv_directory, csv_filename)

    with open(csv_filepath, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        row1 = ['Nom du point', 'Heure']
        row2 = [nomPoint, current_time.strftime('%H-%M-%S')]
        for km, m in measures.items():
            row1.append(km)
            row2.append(m)
        writer.writerow(row1)
        writer.writerow(row2)






