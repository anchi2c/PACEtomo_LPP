#!Python
nav_id = 42 
tgts_path = 'Position_22_tgts.txt'

##### END INPUT ####
realignToItem   = False     # Use SerialEM's RealignToItem routine instead of simple image realignment (was default in PACEtomo <=v1.9.1)
beamTiltComp    = False      # use beam tilt compensation (uses coma vs image shift calibrations)
ronchiC3Offset     = -20          # added to ReportImageDistanceOffset before Trial shot
debug = False
### settings from cal_util
import sys
sys.path.insert(0, 'C:\\Users\\VALUEDGATANCUSTOMER\\Desktop\\anchi\\PACEtomo_LPP')
sys.path.insert(0, 'C:\Program Files\SerialEM\PythonModules')
import os
import serialem as sem
import ronchi_sem_lib
import cal_util
import display_util

sem.GoToLowDoseArea("T")
working_dir = sem.ReportDirectory()
cal_dir, session_name = cal_util.getCalibrationsDir(working_dir)
ronchi_sem_lib.ronchiC3Offset = ronchiC3Offset
optics = ronchi_sem_lib.getOpticsKey(ronchiC3Offset)
ref_correct_ks, ronchi_c3_offset = cal_util.readCalibration('ronchi_ref_ks', cal_dir, optics)
ref_phases, ronchi_c3_offset = cal_util.readCalibration('ronchi_ref_phase', cal_dir, optics)

def log(text, color=0, style=0):
    pass

def doLafis(is_x,is_y):
    pass

def restoreLafis():
    pass

def ronchi_before_preview_align(label):
    print(label)

def alignTo(buffer, debug=False):
    sem.AlignTo(buffer, 0, 0, 0, int(debug))
    if debug:
        try:
            sem.AddBufToStackWindow("A", 0, 0, 0, 0, "CC") #M #S [#B] [#O] [title]
        except AttributeError:
            # Show CC briefly, then switch back to aligned buffer for buffer shift
            sem.Delay(1, "s")
        sem.Copy("B", "A")
        sem.AlignTo(buffer)

def realignTo(nav_id=None, target=None):
    if target is not None and not realignToItem:
        # Move stage to target position
        sem.MoveStageTo(float(target["stageX"]), float(target["stageY"]))
        if "viewfile" in target.keys():
            print(f'reading {target["viewfile"]}')
            sem.ReadOtherFile(0, "O", target["viewfile"]) # reads view file for first AlignTo instead
            display_util.addImage(sem.bufferImage('O'))
            is_x, is_y, *_ = sem.ReportImageShift()
            print('before goto_dose_area_V', is_x, is_y)
            sem.GoToLowDoseArea("V")
            is_v_x, is_v_y, *_ = sem.ReportImageShift()
            print('goto_dose_area_V', is_v_x, is_v_y)
            sem.SetImageShift(0, 0)
            sem.SetImageShift(is_x, is_y)
            sem.V()
            display_util.addImage(sem.bufferImage('A'))
            alignTo("O", debug)
            ASX, ASY = sem.ReportAlignShift()[4:6]
            log(f"Alignment (View) error in X | Y: {round(ASX, 0)} nm | {round(ASY, 0)} nm")
            print(f"Alignment (View) error in X | Y: {round(ASX, 0)} nm | {round(ASY, 0)} nm")
            print('paused at viewfile realignment. press enter to continue')
        if "tgtfile" in target.keys():                
            print(f'reading {target["tgtfile"]}')
            sem.ReadOtherFile(0, "O", target["tgtfile"]) # reads tgt file for first AlignTo instead
            display_util.addImage(sem.bufferImage('O'))
            is_x, is_y, *_ = sem.ReportImageShift()
            print('before ronchi_before_preview align', is_x, is_y)
            sem.GoToLowDoseArea("R")
            sem.SetImageShift(0, 0)
            sem.SetImageShift(is_x, is_y)
            if beamTiltComp:
                doLafis(is_x,is_y)
            ronchi_before_preview_align("initial realign preview (tgtfile)")
            sem.L()
            display_util.addImage(sem.bufferImage('A'))
            print('showing preview after ronchi_before preview')
            if beamTiltComp:
                restoreLafis()
            alignTo("O", debug)
            AISX, AISY, ASX, ASY = sem.ReportAlignShift()[2:6]
            print(is_x, is_y)
            log(f"Alignment (Prev) error in X | Y: {round(ASX, 0)} nm | {round(ASY, 0)} nm")
            print(f"Alignment (Prev) error in X | Y: {round(ASX, 0)} nm | {round(ASY, 0)} nm")
        elif "viewfile" in target.keys():
            # Use align between mags to align preview image to view image
            # If View image was already aligned, take new centered View image at startTilt and use as reference instead
            is_x, is_y, *_ = sem.ReportImageShift()
            sem.GoToLowDoseArea("V")
            sem.SetImageShift(0, 0)
            sem.SetImageShift(is_x, is_y)
            sem.V()
            print('V at start of align between mags',is_x, is_y)
            print('paused at align between mag after V. press enter to continue')
            sem.Copy("A", "O")
            # Check defocus offset
            is_x, is_y, *_ = sem.ReportImageShift()
            sem.GoToLowDoseArea("R") # Switch to R before applying defocus offset to not mess with potential mP/nP offsets between View and Rec
            sem.SetImageShift(0, 0)
            defocus_offset = max(-10, sem.ReportLDDefocusOffset("V"))
            if defocus_offset != 0:
                sem.ChangeFocus(defocus_offset) # Higher defocus for better correlation, but max at 10 to avoid major distortions
            sem.SetImageShift(is_x, is_y)
            print('R IS set to V image shift here', is_x, is_y)
            print('image shift set to R is from V. press enter to continue')
            if beamTiltComp:
                doLafis(is_x,is_y)
            ronchi_before_preview_align("initial realign preview (view to Record)")
            temp_is_x, temp_is_y, *_ = sem.ReportImageShift()
            print(temp_is_x, temp_is_y)
            print('L image taken with Lafis restored. press enter to continue')
            sem.L()

def parseTargets(file_path):
    """Reads targets file."""

    with open(file_path) as f:                                                                         # open last tgts or tgts_run file
        targetFile = f.readlines()

    targets = []
    geoPoints = []
    savedRun = []
    branch = None
    resume = {"sec": 0, "pos": 0}
    for line in targetFile:
        col = line.strip(os.linesep).split(" ")
        if col[0] == "": continue
        if line.startswith("_set") and len(col) == 4:
            if col[1] in globals():
                log(f"WARNING: Read setting from tgts file and overwrite: {col[1]} = {col[3]}")
                globals()[col[1]] = float(col[3])
            else:
                log(f"WARNING: Attempted to overwrite {col[1]} but variable does not exist!")
        elif line.startswith("_bset") and len(col) == 4:
            if col[1] in globals():
                val = True if col[3].lower() in ["true", "yes", "y", "on"] else False
                log(f"WARNING: Read setting from tgts file and overwrite: {col[1]} = {val}")
                globals()[col[1]] = val
            else:
                log(f"WARNING: Attempted to overwrite {col[1]} but variable does not exist!")                
        elif line.startswith("_spos"):
            resume["sec"] = int(col[2].split(",")[0])
            resume["pos"] = int(col[2].split(",")[1])
        elif line.startswith("_tgt"):
            targets.append({})
            branch = None
        elif line.startswith("_pbr"):
            savedRun.append([{},{}])
            branch = 0
        elif line.startswith("_nbr"):
            branch = 1
        elif line.startswith("_geo"):
            geoPoints.append({})
            branch = "geo"
        else:
            if branch is None:
                targets[-1][col[0]] = col[2]
            elif branch == "geo":
                geoPoints[-1][col[0]] = float(col[2])
            else:
                savedRun[-1][branch][col[0]] = col[2]
    if savedRun == []: savedRun = False
    return targets, savedRun, resume, geoPoints

def convertValueType(target):
    float_keys = ['SSX','SSY','stageX','stageY','SPACEscore']
    bool_keys = ['skip']

    for k in float_keys:
        try:
            target[k] = float(target[k])
        except:
            if target[k] == 'None':
                target[k] = None
            else:
                raise
    for k in bool_keys:
        if target[k]=='False':
            target[k] = False
        else:
            target[k] = True
    return target

targets = parseTargets(tgts_path)

target = targets[0][0]
target = convertValueType(target)

sem.SetImageShift(0,0)

realignTo(nav_id=20, target=target)
display_util.showImages(ncols=4, panel_width=4, cmap='gray',title=tgts_path, savefig=True)
