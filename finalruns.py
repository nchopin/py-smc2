#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import subprocess

#nruns = 5

#for i in range(nruns):
#    try:
#        subprocess.call(["python", "launch.py", "userfiles/BSMCSVonesynth.py"])
#    except:
#        pass
#    try:
#        subprocess.call(["python", "launch.py", "userfiles/PMCMCSVonesynth.py"])
#    except:
#        pass

for i in range(1):
    try:
        subprocess.call(["python", "launch.py", "userfiles/BSMCSVonesynthotherh1.py"])
    except:
        pass
    try:
        subprocess.call(["python", "launch.py", "userfiles/BSMCSVonesynthotherh2.py"])
    except:
        pass

#nruns = 3
#processes = []
#for j in range(nruns):
#    processes.append(subprocess.Popen(["python", "launch.py", "userfiles/SMC2SVmultiSP500.py"]))
#for j in range(nruns):
#     processes[j].wait()
#
#processes = []
#for j in range(nruns):
#    processes.append(subprocess.Popen(["python", "launch.py", "userfiles/SMC2SVoneSP500.py"]))
#for j in range(nruns):
#     processes[j].wait()
#
#processes = []
#for j in range(nruns):
#    processes.append(subprocess.Popen(["python", "launch.py", "userfiles/SMC2SVfixedrhoSP500.py"]))
#for j in range(nruns):
#     processes[j].wait()
#
