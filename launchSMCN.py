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
        array, float32, int32, cov, load, isnan, isinf, newaxis, \
        transpose
from numpy import max as numpymax
from numpy import min as numpymin
from numpy import sum as numpysum
from scipy.stats import norm
from snippets.localfolder import get_path
from src.parallelSIRs import ParallelSIRs
import subprocess

random.seed(17)
MODEL = "locallevel"

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
#xmodule.modelx.generateData(5000, xmodule.modelx.parameters, savefilename = "/tmp/txt.txt")
syntheticdatasetpath = os.path.join(THISPATH, "data/%s-syntheticdata.R" % MODEL)
if os.path.isfile(syntheticdatasetpath):
    xmodule.modelx.loadData(syntheticdatasetpath)
xmodule.modelx.model_states = xmodule.modelx.model_states[:, newaxis]
print xmodule.modelx.model_obs.shape
nbparameters = thetamodule.modeltheta.parameterdimension
Nx = 100
Ntheta = 1000
T = min(500, xmodule.modelx.model_obs.shape[0])
model_obs = xmodule.modelx.model_obs[0:T, :]
model_states = xmodule.modelx.model_states[0:T, :]

theta = xmodule.modelx.parameters.copy()
theta = theta.reshape(nbparameters, 1)
thetaparticles = repeat(theta, Ntheta).reshape(nbparameters, Ntheta)

print "parallel SIRs starting..."
import cProfile
cProfile.run("""
parallelSIRs = ParallelSIRs(Nx, thetaparticles, model_obs, xmodule.modelx, saveLL = True, verbose = True)
parallelSIRs.first_step()
parallelSIRs.next_steps()
""", "prof")
print "done."
import pstats
p = pstats.Stats('prof')
p.sort_stats("time").print_stats(15)

SMCLogLikelihoods = parallelSIRs.getTotalLogLike()
results = {"SMCLL": SMCLogLikelihoods,
           "observations": model_obs,
           "allLL": parallelSIRs.allLL}
resultsfolder = os.path.join(THISPATH, "results")
basename = "parallelSMC"
basename = basename + "T%iNx%iNtheta%i" % (T, Nx, Ntheta)
tryname = basename
RDatafile = os.path.join(resultsfolder, basename)
counter = 0
#while os.path.isfile(os.path.join(resultsfolder, tryname + ".RData")):
#    counter += 1
#    tryname = basename + "(%i)" % counter
print "results in %s" % resultsfolder
RDatafile = os.path.join(resultsfolder, tryname + ".RData")
pdffile = os.path.basename(RDatafile.replace(".RData", ".pdf"))
plotresultsfile = os.path.join(resultsfolder, RDatafile.replace(".RData", "-plots.R"))
try:
    import rpy2
    from snippets.pickletoRdata import dictionary2RDataWithRPY 
    dictionary2RDataWithRPY(results, RDatafilename = RDatafile)
except ImportError:
    print "I'd recommend installing rpy2 for faster saving of the results in RData..."
    from snippets.pickletoRdata import dictionary2RDataWithoutRPY
    dictionary2RDataWithoutRPY(results, RDatafilename = RDatafile)
Rcode = \
"""
rm(list = ls())
gc()
graphics.off()
resultsfolder <- "%(resultspath)s/"
resultsfile <- "%(RDatafile)s"
pdffile <- "%(pdffile)s"
library(ggplot2)
setwd(resultsfolder)
cat("working with file", resultsfile, "...\n")
load(file = resultsfile)
pdf(file = pdffile, useDingbats = FALSE, title = "parallel SMC results")
""" % {"RDatafile": RDatafile, "resultspath": resultsfolder, \
        "pdffile": pdffile}
Rcode += xmodule.modelx.RLinearGaussian + "\n"
Rcode += \
"""
KF <- function(observations, somedlm){
  # the notation comes from the package dlm
  T <- length(observations)
  m <- rep(0, T + 1); C <- rep(1, T + 1)
  a <- rep(0, T); R <- rep(0, T)
  f <- rep(0, T); Q <- rep(0, T)
  m[1] <- somedlm$m0; C[1] <- somedlm$C0
  for (t in 1:T){
    a[t] <- somedlm$GG * m[t]
    R[t] <- somedlm$GG * C[t] * somedlm$GG + somedlm$W
    f[t] <- somedlm$FF * a[t]
    Q[t] <- somedlm$FF * R[t] * somedlm$FF + somedlm$V
    m[t+1] <- a[t] + R[t] * somedlm$FF * (1 / Q[t]) * (observations[t] - f[t])
    C[t+1] <- R[t] - R[t] * somedlm$FF * (1 / Q[t]) * somedlm$FF * R[t]
  }
  return(list(observations = observations, NextObsMean = f, NextObsVar = Q,
              NextStateMean = a, NextStatevar = R,
              FiltStateMean = m[2:(T+1)], FiltStateVar = C[2:(T+1)]))
}
getLoglikelihood <- function(KFresults){
  IncrLogLike <- log(dnorm(KFresults$observations, 
            mean = KFresults$NextObsMean, 
            sd = sqrt(KFresults$NextObsVar)))
  loglikelihood <- sum(IncrLogLike)
  return(list(IncrLogLike = IncrLogLike, loglikelihood = loglikelihood))
}
KFLL <- function(observations, dlm){
  KFres <- KF(observations, dlm)
  return(getLoglikelihood(KFres)$loglikelihood)
}
trueLogLikelihood <- KFLL(observations, dlm)
trueIncrlogLikelihood <- getLoglikelihood(KF(observations, dlm))$IncrLogLike
trueCumLogLikelihood <- cumsum(trueIncrlogLikelihood)
"""
Rcode += \
"""
g <- qplot(x = exp(SMCLL), geom = "blank") + geom_histogram(aes(y = ..density..)) + geom_density(alpha = 0.2, fill = "blue")
g <- g + geom_vline(xintercept = exp(trueLogLikelihood), size = 2, colour = "green", linetype = 2)
g <- g + geom_vline(xintercept = mean(exp(SMCLL)), size = 2, colour = "red", linetype = 3)
g <- g + xlab("likelihoods (log scale)") + scale_x_log()
print(g)
SMCcumloglikelihoods <- apply(allLL, 2, cumsum)
SMCcumlikelihoods <- exp(SMCcumloglikelihoods)
SMCcumlikelihoods <- SMCcumlikelihoods / exp(trueCumLogLikelihood)
means <- apply(SMCcumlikelihoods, 1, mean)
vars <- apply(SMCcumlikelihoods, 1, var)
T <- length(observations)
g <- qplot(x = 1:T, y = vars, geom = "line")
g <- g + xlab("time") + ylab("normalised variance")
print(g)
upper <- means + sqrt(vars)
lower <- means - sqrt(vars)
g <- qplot(x = 1:T, y = means, geom = "line")
g <- g + geom_line(aes(y = upper), colour = "green") +
    geom_line(aes(y = lower), colour = "green")
g <- g + xlab("time") + ylab("likelihood estimates")
print(g)
"""
Rcode += \
"""
dev.off()
""" + "\n"
print "R code to plot results in %s" % os.path.join(resultsfolder, plotresultsfile)
f = open(plotresultsfile, "w")
f.write(Rcode)
f.close()
subprocess.call(["R", "CMD", "BATCH", "--vanilla", plotresultsfile, "/dev/null"])


#if hasattr(xmodule.modelx, "RLinearGaussian"):
#    plotter.setDLM(xmodule.modelx.RLinearGaussian)
#RDatafile = "/tmp/testSMCres" + ".RData"
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
#subprocess.call(["R", "CMD", "BATCH", "--vanilla", plotter.plotresultsfile, "/dev/null"])




