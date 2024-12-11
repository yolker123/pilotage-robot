def find_USB_device():
    myports = [tuple(p) for p in list(serial.tools.list_ports.comports())]
    usb_port_list = [p[0:2] for p in
                     myports]  # Prendre p[0], p[1], no need p[2], because myports has 3 string parts.use print to see the value of myports.
    # for ex. myports[0]=('COM4', 'Périphérique série USB (COM4)', 'USB VID:PID=8087:0ACA SER=05022016'), only take first 2 parts.
    # p[0:2] take only p[0] ET P[1] . no p[2].   0<=i<2
    return usb_port_list