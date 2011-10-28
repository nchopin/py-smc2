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
import os, os.path, imp, sys, shutil
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, var, average, prod, log, sum, repeat, \
        array, float32, int32, cov, newaxis, genfromtxt
from snippets.localfolder import get_path
from src.singleSIR import SingleSIR
import userfileSMC as userfile

userfilefile = open("userfileSMC.py", "r")
userfilecontent = userfilefile.read()
userfilefile.close()

if not(userfile.RANDOMSEED):
    random.seed(127)

THISPATH = get_path()

xmodulename = userfile.MODEL + "x"
thetamodulename = userfile.MODEL + "theta"
modelfolder = os.path.join(THISPATH, "models")
sys.path.append(modelfolder)
f, filename, description = imp.find_module(xmodulename)
xmodule = imp.load_module("xmodule", f, filename, description)
f, filename, description = imp.find_module(thetamodulename)
thetamodule = imp.load_module("thetamodule", f, filename, description)

if userfile.DATASET == "synthetic":
    datasetname = userfile.MODEL.replace("CUDA", "")
    syntheticdatasetpath = os.path.join(THISPATH, "data/%s-syntheticdata.R" % datasetname)
    if os.path.isfile(syntheticdatasetpath):
        xmodule.modelx.loadData(syntheticdatasetpath)
        states = genfromtxt(syntheticdatasetpath.replace(".R", "-states.R"))
        if len(states.shape) == 1:
            states = states[:, newaxis]
        xmodule.modelx.model_states = states
    else:
        print "Synthetic data set not found, creating it..."
        xmodule.modelx.generateData(20000, xmodule.modelx.parameters, \
                savefilename = syntheticdatasetpath)
elif not(userfile.DATASET == "synthetic"):
    try:
        xmodule.modelx.loadData(os.path.join(THISPATH, "data/%s.R" % userfile.DATASET))
    except:
        raise ImportError("ERROR: file %s.R not found in the data/ subfolder")


### Parameters, Models
#if userfile.METHOD == "SMC2":
#    algorithmparameters = {"Ntheta": userfile.NTHETA, "Nx": userfile.NX, \
#            "rwvariance": userfile.RWVARIANCE, \
#            "ESSthreshold": userfile.ESSTHRESHOLD, "nbmoves": \
#            userfile.NBMOVES, "proposalkernel": userfile.PROPOSALKERNEL, \
#            "dynamicNx": userfile.DYNAMICNX, "dynamicNxThreshold": \
#            userfile.DYNAMICNXTHRESHOLD, "NxLimit": userfile.NXLIMIT, \
#            "filtering": userfile.FILTERING, "smoothing": userfile.SMOOOTHING, \
#            "smoothingtimes": userfile.SMOOTHINGTIMES, "storesmoothingtime": userfile.STORESMOOTHINGTIME}
#elif userfile.METHOD == "SOPF":
#    algorithmparameters = {"N": userfile.NSOPF}
#elif userfile.METHOD == "BSMC":
#    algorithmparameters = {"N": userfile.NBSMC, "smooth": userfile.BSMCsmooth, \
#            "ESSthreshold": userfile.BSMCESSTHRESHOLD}
#elif userfile.METHOD == "adPMCMC":
#    algorithmparameters = {"N": userfile.NPMCMC, "nbiterations": userfile.TPMCMC}
#else:
#    raise ValueError("ERROR: METHOD should be SMC2, SOPF, BSMC or adPMCMC.")
#

userfile.T = min(userfile.T, xmodule.modelx.model_obs.shape[0])
print "using %i observations out of %i..." % (userfile.T, xmodule.modelx.model_obs.shape[0])
model_obs = xmodule.modelx.model_obs[0:userfile.T]

### hyperparameters might depend on the data
#thetamodule.modeltheta.hyperparameters = thetamodule.modeltheta.updateHyperParam(thetamodule.modeltheta.hyperparameters, model_obs)

model = {"modelx": xmodule.modelx, "modeltheta": thetamodule.modeltheta, "observations": model_obs}

theta = xmodule.modelx.parameters.copy()
print "using the following parameters:"
print theta

#if len(userfile.SAVINGTIMES) > 0:
#    if userfile.METHOD == "adPMCMC":
#        print "Saving times are provided but the adaptive PMCMC algorithm cannot perform sequential analysis."
#    elif (sum(array(userfile.SAVINGTIMES) >= userfile.T) > 0):
#            raise ValueError("ERROR: cant have a saving time bigger than T. Fix the SAVINGTIMES list.")

# launching algorithm
if userfile.PROFILING:
    import cProfile
    counter = 0
    tempproffile = "/tmp/prof(%i)" % counter
    while os.path.isfile(tempproffile):
        counter += 1
        tempproffile = "/tmp/prof(%i)" % counter
    cProfile.run("""\
SMCresults = SingleSIR(userfile.NX, theta, model_obs, xmodule.modelx, verbose = True)
    """, tempproffile)
    import pstats
    p = pstats.Stats(tempproffile)
    p.sort_stats("cumulative").print_stats(10)
    p.sort_stats("time").print_stats(10)
else:
    singleSIR = SingleSIR(userfile.NX, theta, model_obs, xmodule.modelx, verbose = True)

if userfile.USESUBFOLDERS:
    resultsfolder = os.path.join(THISPATH, "results", userfile.MODEL, userfile.DATASET)
    if not(os.path.isdir(resultsfolder)):
        print "creating folder %s" % resultsfolder
        modelsubfolder = os.path.join(THISPATH, "results", userfile.MODEL)
        if os.path.isdir(modelsubfolder):
            datasetsubfolder = os.path.join(modelsubfolder, userfile.DATASET)
            if not(os.path.isdir(datasetsubfolder)):
                os.mkdir(datasetsubfolder)
        else:
            os.mkdir(modelsubfolder)
            datasetsubfolder = os.path.join(modelsubfolder, userfile.DATASET)
            os.mkdir(datasetsubfolder)
else:
    resultsfolder = os.path.join(THISPATH, "results")
if len(userfile.RESULTSFILENAME) > 0:
    basename = userfile.RESULTSFILENAME
else:
    if not(userfile.USESUBFOLDERS):
        prefix = "%s-%s-%s" % ("SMC", userfile.MODEL, userfile.DATASET)
    else:
        prefix = "SMC"
    basename = "%s-T%i-N%i" % (prefix, userfile.T, userfile.NX)
if userfile.REPLACEFILE:
    resultsfile = os.path.join(resultsfolder, basename)
else:
    counter = 0
    currentname = "%s(%i)" % (basename, counter)
    filelist = os.listdir(resultsfolder)
    while len([x for x in filelist if x.startswith(currentname)]) > 0:
        counter += 1
        currentname = "%s(%i)" % (basename, counter)
    resultsfile = os.path.join(resultsfolder, currentname)
print "Saving results in %s*..." % resultsfile

copyuserfile = resultsfile + "-userfile.py"
print "Copying userfile to %s..." % copyuserfile
userfilefile = open(copyuserfile, "w")
userfilefile.write(userfilecontent)
userfilefile.close()

if userfile.PROFILING:
    profilepath = resultsfile + "-profile.txt"
    print "Copying profile to %s..." % profilepath
    profilefile = open(profilepath, "w")
    p = pstats.Stats(tempproffile, stream = profilefile)
    p.sort_stats("cumulative").print_stats(50)
    p.sort_stats("time").print_stats(50)
    profilefile.close()
    os.remove(tempproffile)

resultsDict = SMCresults.getResults()

if "cpickle" in userfile.RESULTSFILETYPE:
    import cPickle
    f = open(resultsfile + ".cpickle", "w")
    cPickle.dump(resultsDict, f)
    f.close()

if "RData" in userfile.RESULTSFILETYPE:
    RDatafile = resultsfile + ".RData"
    try:
        import rpy2
        from snippets.pickletoRdata import dictionary2RDataWithRPY 
        dictionary2RDataWithRPY(resultsDict, RDatafilename = RDatafile)
    except ImportError:
        print "I'd recommend installing rpy2 for faster saving of the results in RData..."
        from snippets.pickletoRdata import dictionary2RDataWithoutRPY
        dictionary2RDataWithoutRPY(resultsDict, RDatafilename = RDatafile)

if userfile.GENERATERFILE:
    if not("RData" in userfile.RESULTSFILETYPE):
        raise ImportError('no RData file to plot, you should add "RData" to the RESULTSFILETYPE list.')
    pdffile = RDatafile.replace(".RData", ".pdf")
    print "plotting results in pdf file: %s ..." % pdffile
    from src.plotresultsSMC import PlotResultsSMC
    plotter = PlotResultsSMC(resultsfolder, RDatafile)
    plotter.close()
    if userfile.PLOT:
        import subprocess
        subprocess.call(["R", "CMD", "BATCH", "--vanilla", plotter.plotresultsfile, "/dev/null"])

print "...done! Bye."



