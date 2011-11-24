###################################################
#    This file is part of py-smc2.
#    http://code.google.com/p/py-smc2/
#
#    py-smc2 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    py-smc2 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with py-smc2.  If not, see <http://www.gnu.org/licenses/>.
###################################################

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
from src.SIR import SIR
from src.plot.plotresultsSMC import PlotResultsSMC
import subprocess

#random.seed(17)
MODEL = "simplestmodel"

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
xmodule.modelx.generateData(1000, xmodule.modelx.parameters, savefilename = "/tmp/txt.txt")

nbparameters = thetamodule.modeltheta.parameterdimension
Nx = 10**2
T = min(100, xmodule.modelx.model_obs.shape[0])
model_obs = xmodule.modelx.model_obs[0:T, :]
model_states = xmodule.modelx.model_states[0:T, :]

theta = xmodule.modelx.parameters.copy()
print theta

SAVINGTIMES = [T / 2, T]
import cProfile
cProfile.run("""
singleSIR = SIR(Nx, theta, model_obs, xmodule.modelx, savingtimes = SAVINGTIMES, verbose = False)
""", "prof")
import pstats
p = pstats.Stats('prof')
p.sort_stats("time").print_stats(10)
results = singleSIR.getResults()
#results.update({"truestates": xmodule.modelx.model_states[0:T,:]})
#RDatafile = "/tmp/SMCfiltering" + ".RData"
#try:
#    import rpy2
#    from snippets.pickletoRdata import dictionary2RDataWithRPY 
#    dictionary2RDataWithRPY(results, RDatafilename = RDatafile)
#except ImportError:
#    print "I'd recommend installing rpy2 for faster saving of the results in RData..."
#    from snippets.pickletoRdata import dictionary2RDataWithoutRPY
#    dictionary2RDataWithoutRPY(results, RDatafilename = RDatafile)
#resultsfolder = os.path.join(THISPATH, "results")
#plotter = PlotResultsSMC(RDatafile, resultsfolder)
#if hasattr(xmodule.modelx, "RLinearGaussian"):
#    plotter.setDLM(xmodule.modelx.RLinearGaussian)
#if hasattr(xmodule.modelx, "RLinearGaussian"):
#    plotter.addKalmanComparison()
#plotter.close()
#print "R code to plot results in %s" % plotter.plotresultsfile
#subprocess.call(["R", "CMD", "BATCH", "--vanilla", plotter.plotresultsfile, "/dev/null"])



