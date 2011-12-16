#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import subprocess

nruns = 2

for i in range(nruns):
    try:
        subprocess.call(["python", "launch.py", "userfiles/BSMCSVonesynth.py"])
    except:
        pass

for i in range(nruns):
    try:
        subprocess.call(["python", "launch.py", "userfiles/BSMCSVonesynth.py"])
    except:
        pass
