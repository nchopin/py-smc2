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
from src.adpmcmc import AdaptivePMCMC
from src.bsmc import BSMC
from snippets.localfolder import get_path

from models.locallevelx import modelx
from models.localleveltheta import modeltheta

random.seed(17)

THISPATH = get_path()
syntheticdatasetpath = os.path.join(THISPATH, "testdata.R")

T = 100
modelx.generateData(T, modelx.parameters, \
    savefilename = syntheticdatasetpath)

model_obs = modelx.model_obs[0:T]

### hyperparameters might depend on the data
modeltheta.hyperparameters = modeltheta.updateHyperParam(modeltheta.hyperparameters, model_obs)
model = {"modelx": modelx, "modeltheta": modeltheta, "observations": model_obs}

algorithmparameters = {"N": 10000}
algo = BSMC(model, algorithmparameters, autoinit = False, savingtimes = [])
algo.first_step()
algo.next_steps()
results = algo.getResults()
#print results["thetahistory"]
#print algo.savingtimes
#print algo.__dict__
#print algo.__dict__["thetaparticles"]
#print algo.__dict__["transformedthetaparticles"]
#print modeltheta.transform(algo.__dict__["thetaparticles"])
#print modeltheta.untransform(algo.__dict__["transformedthetaparticles"])


resultsfile = os.path.join(THISPATH, "results.cpickle")
import cPickle
f = open(resultsfile, "w")
cPickle.dump(results, f)
f.close()

from snippets.pickletoRdata import pickle2RData
RDatafile = resultsfile.replace(".cpickle", ".RData")
print "...and saving results in %s..." % RDatafile
pickle2RData(resultsfile)
os.remove(resultsfile)

