"""Save a screenshot on the device and copy it to the local machine/environment."""

from tm_devices import DeviceManager
from tm_devices.drivers import MSO6B

with DeviceManager(verbose=True) as dm:
    # Add a scope
    scope: MSO6B = dm.add_scope("172.16.115.218")
    #cwd_command = 'FILESystem:CWD "C:/Users/Public/Tektronix/TekScope/Screenshots"'


    waveform_command = 'SAVEON:WAVEFORM:DEST "C:\\Users\\Public\\Tektronix\\TekScope"'
    scope.write(waveform_command)
    scope.commands.save.waveform.write('CH1, "test.csv"')

    #waveform_command = 'SAVEON:FILE:DEST "C:\\Users\\Public\\Tektronix\\TekScope\\Screenshots"'
    #scope.write(waveform_command)
