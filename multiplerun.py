#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import subprocess

nruns = 10
ncores = 1
for i in range(nruns):
    processes = []
    for j in range(ncores):
        try:
            processes.append(subprocess.Popen(["python", "launch.py", "userfiles/athletics.py"]))
        except:
            pass
    for j in range(ncores):
         processes[j].wait()
#
#for i in range(nruns):
#    processes = []
#    for j in range(ncores):
#        try:
#            processes.append(subprocess.Popen(["python", "launch.py", "userfiles/PMCMCSVonesynth.py"]))
#        except:
#            pass
#    for j in range(ncores):
#         processes[j].wait()
