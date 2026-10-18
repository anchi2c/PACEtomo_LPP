#!Python
# need sem
import sys
import platform
if platform.system() == 'Windows':
    is_simu = False
    required_sys_paths = [    'C:\Program Files\SerialEM\PythonModules',
                    'C:\\Users\\VALUEDGATANCUSTOMER\\Desktop\\anchi\\PACEtomo_LPP',
                ]
    for p in required_sys_paths:
        if p not in sys.path:
            sys.path.insert(0,p)
    import serialem as sem
else:
    is_simu = True
    print('testing on Mac/Linux with simulator')
    import sem_simulator as sem

def breakpoint():
    """Breakpoint for debugging in SerialEM."""

    while not sem.KeyBreak():
        sem.Delay(0.1, "s")
    for i in range(5):
        if sem.KeyBreak("d"):
            dumpVars()
            break
        sem.Delay(0.1, "s")

def log(text, color=0, style=0):
    if text.startswith("DEBUG:") and not debug:
        return
    if text.startswith("NOTE:"):
        color = 4
    elif text.startswith("WARNING:"):
        color = 5
    elif text.startswith("ERROR:"):
        color = 2
        style = 1 
    elif text.startswith("DEBUG:"):
        color = 1
        if breakpoints:
            breakpoint()

    if sem.IsVersionAtLeast("40200", "20240205"):
        sem.SetNextLogOutputStyle(style, color)
    sem.EchoBreakLines(text)
