"""Save a screenshot on the device and copy it to the local machine/environment."""

from tm_devices import DeviceManager
from tm_devices.drivers import MSO6B

with DeviceManager(verbose=True) as dm:
    # Add a scope
    scope: MSO6B = dm.add_scope("172.16.115.218")
    #cwd_command = 'FILESystem:CWD "C:/Users/Public/Tektronix/TekScope/Screenshots"'

    scope.write('FILESystem:COPy "C:/Users/Public/Tektronix/TekScope/Screenshots/logScreenshot_tektronix_20250124_150234/C1","//PCROBOT/Users/isen/PycharmProjects/pilotage-robot/Code GUI version 202409/202409119/Measure/TekScope/Screenshots/logScreenshot_tektronix_20250124_150234/"')