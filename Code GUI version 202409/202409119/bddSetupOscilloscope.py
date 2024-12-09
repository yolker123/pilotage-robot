"""
 * @file        bddSetupOscilloscope.py
 * @brief       Contains differents possible setup for the oscilloscope to set up
 * @author      Lisa Duterte modified by Hamid Ajouaou & Romeo Botuli
 * @version     0.2
 * @date        2023 & 2024
"""

"""
* @brief  This code make the setup of the Oscilloscope Window (Grid, Channels, Scales, Offset, Measures)
          with the record of the parameters, all in the acquisition1() function.
"""
import os
import datetime

global measure_config
measure_config = list()
global wf_img_config
wf_img_config = {
    'C1' : {
        "img" : False,
        "wf" : False
    },
    'C2' : {
        "img" : False,
        "wf" : False
    },
    'C3' : {
        "img" : False,
        "wf" : False
    },
    'C4' : {
        "img" : False,
        "wf" : False
    }
}
"""[
    
    {
        "ch" : "C1",
        "info" : "Frequency"
    },
    {
        "ch" : "C2",
        "info" : "RMS"
    },
    {
        "ch" : "C1",
        "info" : "pkpk"
    },
    {
        "ch" : "C2",
        "info" : "pkpk"
    }
]"""

def acquisition1():
    #Set axes parameters
    ver_scale =  {
        'Channel_1' : 0.2,
        'Channel_2' : 0.1,
        #'Channel_3' : 0.005
    }

    ver_offset =  {
        #'Channel_1' : -0.005,
        'Channel_1' : 1,
        'Channel_2' : 0.0,
        #'Channel_3' : 0.05
    }

    # set analog channels parameters (label, ver_scale, ver_offset, bandwidth, coupling)
    analog_channels = {
        1: ('Channel_1',  ver_scale['Channel_1'],   ver_offset['Channel_1'],  '20kHz', 'DC1M'),
        2: ('Channel_2',  ver_scale['Channel_2'],   ver_offset['Channel_2'],  '20kHz', 'DC1M'),
        #3: ('Channel_3',  ver_scale['Channel_3'],   ver_offset['Channel_3'],  '20kHz', 'DC1M')
    }
    # set measurement channels options :
    #   - Frequency - MAX - MIN
    #   - RMS (root mean square: moyenne quadratique)
    #   - pkpk (pic-à-pic),
    #   - Threshold (dépassement de seuil)
    #   - RiseTime - FallTime
    #   - PropagationDelay
    measurement_channels = { 
        1:  ('C1', 'Frequency'),
        2:  ('C2', 'RMS'),
        #3:  ('C3', 'RMS'),
        4:  ('C1', 'pkpk'),
        5:  ('C2', 'pkpk'),
        #6:  ('C3', 'pkpk')
    }
    return (ver_offset,analog_channels,measurement_channels)

def formatMeasureConfig(conf):
    from functools import cmp_to_key
    def compare(a,b) :
        ch_a = a.get("ch")
        ch_b = b.get("ch")
        info_a = a.get("info")
        info_b = b.get("info")
        if ch_a == ch_b and info_a == info_b: return 0
        if ch_a > ch_b or (info_a > info_b and ch_a == ch_b): return 1
        return -1
    config = sorted(conf, key=cmp_to_key(compare))
    measurement_channels = {}
    for i in range(len(config)):
        measurement_channels[i+1]=(config[i]["ch"],config[i]["info"])
    return measurement_channels

