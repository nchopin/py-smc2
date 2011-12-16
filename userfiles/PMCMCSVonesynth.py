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

# Requirements: python, numpy, a functional scipy.weave # (hence a C compiler installed and configured)
# like gcc (mac os users should install xcode)
# Optional but useful: R, ggplot2 (for the plots), rpy2

# Once you've chosen the paramerers, launch the script "launch.py" by typing
# :~$ python launch.py pathtotheuserfile

RANDOMSEED = True 
##
MODEL = "SVonefactor"
T = 1000
DATASET = "synthetic"
METHOD = "adPMCMC"
##
GENERATERFILE = True
PLOT = True
##
NTHETA = 1000
NX = 100
DYNAMICNX = True
PROPOSALKERNEL = "independent"
ESSTHRESHOLD = 0.5
##
DYNAMICNXTHRESHOLD = 0.2
NXLIMIT = 5000
RWVARIANCE = 0.1
NBMOVES = 1
##
NSOPF = 100000
##
NBSMC = 20000
BSMCsmooth = 0.1
##
NPMCMC = 500
TPMCMC = 50000
PMCMCBURNIN = TPMCMC / 10
##
#SAVINGTIMES = [250]
SAVINGTIMES = [250, 500, 750]
#SAVINGTIMES = []
SMOOTHING = False
FILTERING = True
PREDICTION = True
SMOOTHINGTIMES = []
STORESMOOTHINGTIME = 0
##
PROFILING = True
RESULTSFILENAME = ""
REPLACEFILE = False
USESUBFOLDERS = True
RESULTSFILETYPE = ["RData"]

##################################################################
########################## EXPLANATIONS ##########################
##################################################################
# If RANDOMSEED is True, a new seed is used at each run, otherwise it's fixed (for debugging purposes)
##########################
# METHOD could be:
# SMC2: SOPF, BSMC, adPMCMC
##########################
### Models and datasets options
# MODEL is a string, there has to be files named MODELx.py and MODELtheta.py
# in the models/ subfolder. Provided models are: 
# - hiddenAR: most basic linear gaussian SSM
# - locallevel: basic linear gaussian SSM
# - thetalogistic: population model with non linear hidden process (Polansky's parameterization)
# - periodic: highly non linear model from Gordon Salmond and Smith
# - SVonefactor: SV stands for stochastic volatility. This model is the one factor model.
# - SVfixedrho: this one is the two factor model without leverage (rho_1 = rho_2 = 0)
# - SVmultifactor: this one is the two factor model with leverage
# - athletics: SSM with GEV observation function, to model the best times in athletics.
# T: number of observations to be taken into account
# DATASET could either be a string corresponding to a file named
# DATASET.R in the data/ subfolder, or "synthetic", in which case
# data points are simulated from the model.
# Provided datasets are:
# - "SP500recent": 753 data points from the S&P500 index. You can test the
# nonlinear SV model on this data set.
# - "anisus" (18 data points), or "nutria" (120 data points). You can test these datasets
# with the thetalogistic model.
# - "athletics-best-two": (for the athletics records model), taken from http://www.alltime-athletics.com/
##########################
### SMC2 algorithm parameters
# NTHETA: number of theta-particles
# NX: number of x-particles 
# DYNAMICNX: automatic increase of NX:
# PROPOSALKERNEL can be either "randomwalk" or "independent".
# ESSTHRESHOLD: resample-move steps are triggered when the ESS goes below ESSTHRESHOLD.
##########################
### Advanced parameters for SMC2
## DYNAMICNXTHRESHOLD: If DYNAMICNX, acceptance rate limit
## NXLIMIT: If DYNAMICNX, maximum number of x-particles allowed
## If PROPOSALKERNEL == "randomwalk",
## the random walk has variance RWVARIANCE * Cov(theta-particles).
# NBMOVES: Number of moves to perform at each resample-move step
##########################
### SOPF algorithm parameters
# NSOPF: If you want to try the SOPF algorithm to compare the results,
# specify the number of particles here
##########################
### BSMC algorithm parameters
# NBSMC: If you want to try the BSMC algorithm to compare the results,
# specify the number of particles here
# BSMCsmooth: specify the "h" factor used in the kernel density approximation
# of the distribution of the parameters given the data
# should be "slowly decreasing with N", N being the number of particles
##########################
### adaptive PMCMC algorithm parameters
# If you want to try the adaptive PMCMC algorithm to compare the results:
# NPMCMC: Number of x-particles
# TPMCMC: Number of iterations
# PMCMCBURNIN: burn-in (for the plots only)
##########################
### Plot options
# GENERATERFILE: Generate R file to plot results (does not require R)
## PLOT: Generate pdf figures using the R file (requires R and package "ggplot2")
##########################
### Results options
# Specify time at which the theta-particles should be saved
# (useful when doing sequential analysis).
# The final time T = number of observations is automatically saved,
# no need to add it to the list.
### : Should the program report profiling? (Slows the program a little bit).
# RESULTSFILENAME: Specify results file name (without extensions),
# leave empty if an automatic name is preferred, based on the algorithm parameters.
# REPLACEFILE: Replace already existing files; if False,
# a counter will be added to the result file name, so that no existing files are erased.
# USESUBFOLDERS: Use subfolders for each model and each dataset to organize the results.
# If False, use long names and store the results at the root of the results/ folder.
# RESULTSFILETYPE could include "RData" and "cpickle".
# This list is used to save the results either in RData format 
# or in cPickle format or in both formats. RData is required to use the graph programs
# of subfolder rgraphs/.
#RESULTSFILETYPE = ["RData", "cpickle"]
##########################



