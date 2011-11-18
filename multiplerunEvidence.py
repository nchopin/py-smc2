#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import subprocess

nruns = 3
for i in range(nruns):
    processes = []
    processes.append(subprocess.Popen(["python", "launchEvidenceSMC2.py"]))

for i in range(nruns):
    processes = []
    processes.append(subprocess.Popen(["python", "launchEvidenceBSMC.py"]))




