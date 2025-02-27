"""Save a screenshot on the device and copy it to the local machine/environment."""

from tm_devices import DeviceManager
from tm_devices.drivers import MSO6B

from tm_devices.helpers import (
    DMConfigOptions,
    SYSTEM_DEFAULT_VISA_BACKEND,
)
with DeviceManager(verbose=True) as dm:
    # Add a scope
    CONFIG_OPTIONS = DMConfigOptions(
        setup_cleanup=True,  # update the value for this option, all other options will remain untouched
        teardown_cleanup=True,
    )
    device_manager = DeviceManager(verbose=False, config_options=CONFIG_OPTIONS)
    device_manager.visa_library = SYSTEM_DEFAULT_VISA_BACKEND
    scope: MSO6B = device_manager.add_scope("USB0::0x0699::0x0530::C071483::INSTR")

    screenshot_filepath = "TEK00000.SET"
    write_screenshot = 'RECALL:SETUP \"' + screenshot_filepath + '"'
    scope.write(write_screenshot)