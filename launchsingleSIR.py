#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import os, os.path, sys, imp
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32, cov, load, isnan, isinf, newaxis
from numpy import max as numpymax
from numpy import min as numpymin
from numpy import sum as numpysum
from scipy.stats import norm
from snippets.localfolder import get_path
from src.singleSIR import SingleSIR


MODEL = "thetalogistic"


THISPATH = get_path()

xmodulename = MODEL + "x"
thetamodulename = MODEL + "theta"
modelfolder = os.path.join(THISPATH, "models")
sys.path.append(modelfolder)
f, filename, description = imp.find_module(xmodulename)
xmodule = imp.load_module("xmodule", f, filename, description)
f, filename, description = imp.find_module(thetamodulename)
thetamodule = imp.load_module("thetamodule", f, filename, description)

print "Creating data set..."
xmodule.modelx.generateData(10000, xmodule.modelx.parameters, savefilename = "/tmp/txt.txt")

nbparameters = thetamodule.modeltheta.parameterdimension
Nx = 50000
T = min(300, xmodule.modelx.model_obs.shape[0])
model_obs = xmodule.modelx.model_obs[0:T, :]
model_states = xmodule.modelx.model_states[0:T, :]

theta = xmodule.modelx.parameters.copy()
print theta

import cProfile
cProfile.run("""
singleSIR = SingleSIR(Nx, theta, model_obs, xmodule.modelx, verbose = True)
""", "prof")
import pstats
p = pstats.Stats('prof')
p.sort_stats("time").print_stats(3)


import rpy2.robjects as robjects
r = robjects.r
mstates = robjects.FloatVector(model_states[0:T,0])
mobs = robjects.FloatVector(model_obs[0:T,0])
yrange = r.range(mstates)
yrange[0] -= 1
yrange[1] += 1
myrange = robjects.FloatVector([-40, 40])
r("""par(mfrow = c(1, 1))""")
r.plot(x = robjects.FloatVector(array(range(T))), y = mstates, ylim = yrange, xlab = "index", ylab = "x", type = "b", col = "black", lwd = 2.5)
#r.lines(x = robjects.FloatVector(array(range(T))), y = mobs, col = "blue", lwd = 2.5, lty = 2)
#r.plot(my, ylab = "x", type = "l", ylim = myrange, lwd = 2)
approx = robjects.FloatVector(singleSIR.path[:, 0])
#r.points(mx, col = "black", lwd = 2.5)
#r.points(mx, col = "white", lwd = 1.5)
r.lines(x = robjects.FloatVector(array(range(T))), y = approx, col = "black", lwd = 3.5, lty = 4)
r.lines(x = robjects.FloatVector(array(range(T))), y = approx, col = "yellow", lwd = 2, lty = 4)
#for index in range(Nx):
#    print index
#    onetraj = singleSIR.retrieveTrajectory(index)[1:(T+1), :]
#    Ronetraj = robjects.FloatVector(onetraj[:, 0])
#    r.lines(x = robjects.FloatVector(array(range(T))), y = Ronetraj, col = "red", lwd = 0.4)

raw_input("appuyez sur une touche pour fermer la fenetre")







