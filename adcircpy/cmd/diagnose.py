import pathlib
import numpy as np


def parse(logfile):
    logfile = pathlib.Path(logfile).resolve()
    with open(logfile, 'r') as f:
        lines = "".join(f.readlines())
    elmax = list()
    speedmax = list()
    index = list()
    _lines = lines.split(
        '** ERROR: Elevation.gt.ErrorElev, ADCIRC stopping. **\n')
    line0 = "".join(_lines[0]).split('\n')
    for line in line0:
        if "** WARNING: Elevation.gt.WarnElev **" in line:
            elmax.append(float(line.split("AT NODE")[0].split("=")[-1]))
            speedmax.append(
                float(line.split("SPEEDMAX =")[0].split("AT NODE")[-1]))
            index.append(
                np.abs(int(
                    line.split("AT NODE")[-1].split("ON MYPROC")[0].strip())
                )-1)
    return elmax, speedmax, index
