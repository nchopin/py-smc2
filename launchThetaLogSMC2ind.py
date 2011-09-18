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
        array, float32, int32, cov, isnan, load, savez
from src.smcsquare import SMCsquare
from src.sopf import SOPF 
from src.bsmc import BSMC
from src.adpmcmc import AdaptivePMCMC
from snippets.localfolder import get_path
import userfileThetaLogSMC2ind as userfile

userfilefile = open("userfile.py", "r")
userfilecontent = userfilefile.read()
userfilefile.close()

#random.seed(127)

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
    syntheticdatasetpath = os.path.join(THISPATH, "data/%s-syntheticdata.R" % userfile.MODEL)
    if os.path.isfile(syntheticdatasetpath):
        xmodule.modelx.loadData(syntheticdatasetpath)
    else:
        print "Synthetic data set not found, creating it..."
        xmodule.modelx.generateData(20000, xmodule.modelx.parameters, \
                savefilename = syntheticdatasetpath)
elif not(userfile.DATASET == "synthetic"):
    try:
        xmodule.modelx.loadData(os.path.join(THISPATH, "data/%s.R" % userfile.DATASET))
    except:
        raise ImportError("ERROR: file %s.R not found in the data/ subfolder")


## Parameters, Models
if userfile.METHOD == "SMC2":
    algorithmparameters = {"Ntheta": userfile.NTHETA, "Nx": userfile.NX, \
            "rwvariance": userfile.RWVARIANCE, \
            "ESSthreshold": userfile.ESSTHRESHOLD, "nbmoves": \
            userfile.NBMOVES, "proposalkernel": userfile.PROPOSALKERNEL, \
            "dynamicNx": userfile.DYNAMICNX, "dynamicNxThreshold": \
            userfile.DYNAMICNXTHRESHOLD, "NxLimit": userfile.NXLIMIT, \
            "filtering": userfile.FILTERING, "smoothing": userfile.SMOOOTHING, \
            "smoothingtimes": userfile.SMOOTHINGTIMES, "storesmoothingtime": userfile.STORESMOOTHINGTIME}
elif userfile.METHOD == "SOPF":
    algorithmparameters = {"N": userfile.NSOPF}
elif userfile.METHOD == "BSMC":
    algorithmparameters = {"N": userfile.NBSMC, "smooth": userfile.BSMCsmooth, \
            "ESSthreshold": userfile.BSMCESSTHRESHOLD}
elif userfile.METHOD == "adPMCMC":
    algorithmparameters = {"N": userfile.NPMCMC, "nbiterations": userfile.TPMCMC}
else:
    raise ValueError("ERROR: METHOD should be SMC2, SOPF, BSMC or adPMCMC.")

userfile.T = min(userfile.T, xmodule.modelx.model_obs.shape[0])
print "using %i observations out of %i..." % (userfile.T, xmodule.modelx.model_obs.shape[0])
model_obs = xmodule.modelx.model_obs[0:userfile.T]

### hyperparameters might depend on the data
thetamodule.modeltheta.hyperparameters = thetamodule.modeltheta.updateHyperParam(thetamodule.modeltheta.hyperparameters, model_obs)

model = {"modelx": xmodule.modelx, "modeltheta": thetamodule.modeltheta, "observations": model_obs}

if len(userfile.SAVINGTIMES) > 0:
    if userfile.METHOD == "adPMCMC":
        print "Saving times are provided but the adaptive PMCMC algorithm cannot perform sequential analysis."
    elif (sum(array(userfile.SAVINGTIMES) >= userfile.T) > 0):
            raise ValueError("ERROR: cant have a saving time bigger than T. Fix the SAVINGTIMES list.")

# launching algorithm
if userfile.PROFILING:
    import cProfile
    if userfile.METHOD == "SMC2":
        cProfile.run("""\
algo = SMCsquare(model, algorithmparameters, \
dynamicNx = userfile.DYNAMICNX, savingtimes = userfile.SAVINGTIMES)\
    """, "prof")
    elif userfile.METHOD == "SOPF":
        cProfile.run("""\
algo = SOPF(model, algorithmparameters, savingtimes = userfile.SAVINGTIMES)\
    """, "prof")
    elif userfile.METHOD == "BSMC":
        cProfile.run("""\
algo = BSMC(model, algorithmparameters, savingtimes = userfile.SAVINGTIMES)\
    """, "prof")
    elif userfile.METHOD == "adPMCMC":
        cProfile.run("""\
algo = AdaptivePMCMC(model, algorithmparameters)\
    """, "prof")
    import pstats
    p = pstats.Stats('prof')
    p.sort_stats("cumulative").print_stats(10)
    p.sort_stats("time").print_stats(10)
else:
    if userfile.METHOD == "SMC2":
        algo = SMCsquare(model, algorithmparameters, \
                dynamicNx = userfile.DYNAMICNX, savingtimes = userfile.SAVINGTIMES)
    elif userfile.METHOD == "SOPF":
        algo = SOPF(model, algorithmparameters, savingtimes = userfile.SAVINGTIMES)
    elif userfile.METHOD == "BSMC":
        algo = BSMC(model, algorithmparameters, savingtimes = userfile.SAVINGTIMES)
    elif userfile.METHOD == "adPMCMC":
        algo = AdaptivePMCMC(model, algorithmparameters) 

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
        prefix = "%s-%s-%s" % (userfile.METHOD, userfile.MODEL, userfile.DATASET)
    else:
        prefix = "%s" % userfile.METHOD
    if userfile.DYNAMICNX and userfile.METHOD == "SMC2":
        basename = "%s-T%i-dynamicNx%i-Nth%i" % (prefix, userfile.T, userfile.NX, userfile.NTHETA)
    elif not(userfile.DYNAMICNX) and userfile.METHOD == "SMC2":
        basename = "%s-T%i-Nx%i-Nth%i" % (prefix, userfile.T, userfile.NX, userfile.NTHETA)
    elif userfile.METHOD == "SOPF":
        basename = "%s-T%i-N%i" % (prefix, userfile.T, userfile.NSOPF)
    elif userfile.METHOD == "BSMC":
        basename = "%s-T%i-N%i" % (prefix, userfile.T, userfile.NBSMC)
    else:
        basename = "%s-T%i" % (prefix, userfile.T)
if userfile.REPLACEFILE:
    resultsfile = os.path.join(resultsfolder, "%s.cpickle" % basename)
else:
    checkcpickle = os.path.join(resultsfolder, "%s(0).cpickle" % basename)
    checkRData = os.path.join(resultsfolder, "%s(0).RData" % basename)
    counter = 0
    while os.path.isfile(checkcpickle) or os.path.isfile(checkRData):
        counter += 1
        checkcpickle = os.path.join(resultsfolder, "%s(%i).cpickle" % (basename, counter))
        checkRData = os.path.join(resultsfolder, "%s(%i).RData" % (basename, counter))
    resultsfile = checkcpickle
print "Saving results in %s..." % resultsfile

copyuserfile = resultsfile.replace(".cpickle", "-userfile.py")
print "Copying userfile to %s..." % copyuserfile
userfilefile = open(copyuserfile, "w")
userfilefile.write(userfilecontent)
userfilefile.close()

if userfile.PROFILING:
    profilepath = resultsfile.replace(".cpickle", "-profile.txt")
    print "Copying profile to %s..." % profilepath
    profilefile = open(profilepath, "w")
    p = pstats.Stats('prof', stream = profilefile)
    p.sort_stats("cumulative").print_stats(50)
    p.sort_stats("time").print_stats(50)
    profilefile.close()
    os.remove("prof")

resultsDict = algo.getResults()

import cPickle
f = open(resultsfile, "w")
cPickle.dump(resultsDict, f)
f.close()

if "RData" in userfile.RESULTSFILETYPE:
    from snippets.pickletoRdata import pickle2RData
    ## wait RDatafile looks kinda useless now
    RDatafile = resultsfile.replace(".cpickle", ".RData")
    print "...and saving results in %s..." % RDatafile
    pickle2RData(resultsfile)
if not("cpickle" in userfile.RESULTSFILETYPE):
    os.remove(resultsfile)

if userfile.GENERATERFILE:
    if not("RData" in userfile.RESULTSFILETYPE):
        raise ImportError('no RData file to plot, you should add "RData" to the RESULTSFILETYPE list.')
    pdffile = RDatafile.replace(".RData", ".pdf")
    print "plotting results in pdf file: %s ..." % pdffile
    if userfile.METHOD == "SMC2":
        from src.plotresultsSMC2 import PlotResultsSMC2
        plotter = PlotResultsSMC2(resultsfolder, RDatafile)
    elif userfile.METHOD == "adPMCMC":
        from src.plotresultsPMMH import PlotResultsPMMH
        plotter = PlotResultsPMMH(resultsfolder, RDatafile)
    elif userfile.METHOD == "SOPF":
        from src.plotresultsSOPF import PlotResultsSOPF
        plotter = PlotResultsSOPF(resultsfolder, RDatafile)
    elif userfile.METHOD == "BSMC":
        from src.plotresultsBSMC import PlotResultsBSMC
        plotter = PlotResultsBSMC(resultsfolder, RDatafile)
    plotter.setParameters(thetamodule.modeltheta.parameternames)
    plotter.setModelTheta(thetamodule.modeltheta)
    plotter.setModelX(xmodule.modelx)
    plotter.everything()
    if userfile.PLOT:
        import subprocess
        subprocess.call(["R", "CMD", "BATCH", "--vanilla", \
                plotter.plotresultsfile, os.path.join(resultsfolder, "/tmp/R-output.out")])

print "...done! Bye."


