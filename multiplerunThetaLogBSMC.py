#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import subprocess

nbruns = 10
for i in range(nbruns):
    try:
        subprocess.call(["python", "launchThetaLogBSMC.py"])
    except:
        pass

