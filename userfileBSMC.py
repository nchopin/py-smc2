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
# Optional but useful: R, ggplot2 (for the plots), rpy2
# Very optional: CUDA (for the model files with their names finishing by "CUDA")

# If RANDOMSEED is True, a new seed is used at each run, otherwise it's fixed (for debugging purposes)
RANDOMSEED = True

##########################
####  USER PARAMETERS ####
##########################

# Only this file should be modified by most users.
# Once this file is modified and saved, launch.py can be launched with:
# python launch.py
# and the results will
# be saved in the results/ folder. The graph programs in rgraphs/ can
# then be used on the RData files to visualize the results.

##########################
### Model and dataset:
# (see below for possible values)
MODEL = "SVonefactor"
T = 500
DATASET = "synthetic"

##########################
### METHOD
# METHOD could be:
# SMC2: SMC square as described in Chopin, Jacob and Papaspiliopoulos
# SOPF: Self Organizing Particle Filter, as described in Kitagawa
# BSMC: Self Organizing Particle Filter, as described in Kitagawa
# adPMCMC: Adaptive PMCMC as described in Peters, Hosack and Hayes
### (see below for explanations)
##########################

METHOD = "BSMC"

##########################
### Plot options
# Generate R file to plot results (does not require R)
GENERATERFILE = True
## Generate pdf figures using the R file (requires R and package "ggplot2")
PLOT = True

###
##########################
### SMC2 algorithm parameters
# NTHETA: number of theta-particles
# NX: number of x-particles 
# DYNAMICNX: automatic increase of NX:
# PROPOSALKERNEL can be either "randomwalk" or "independent".
# ESSTHRESHOLD: resample-move steps are triggered when the ESS goes below ESSTHRESHOLD.
###
NTHETA = 1000
NX = 125
DYNAMICNX = True
PROPOSALKERNEL = "independent"
ESSTHRESHOLD = 0.5
##########################
### Advancer parameters for SMC2
## If DYNAMICNX, acceptance rate limit:
DYNAMICNXTHRESHOLD = 0.2
## If DYNAMICNX, maximum number of x-particles allowed:
NXLIMIT = 5000
## If PROPOSALKERNEL == "randomwalk",
## the random walk has variance RWVARIANCE * Cov(theta-particles).
RWVARIANCE = 0.1
# Number of moves to perform at each resample-move step:
NBMOVES = 1

### 
##########################
### SOPF algorithm parameters
# If you want to try the SOPF algorithm to compare the results,
# specify the number of particles here:
NSOPF = 100000

##########################
### BSMC algorithm parameters
# If you want to try the BSMC algorithm to compare the results,
# specify the number of particles here:
NBSMC = 250000
# specify the "h" factor used in the kernel density approximation
# of the distribution of the parameters given the data
# should be "slowly decreasing with N", N being the number of particles
BSMCsmooth = 0.1
# specify the ESS threshold (between 0 and 1). 1 means resampling
# steps occur at each steps, while 0 means resampling never occurs
BSMCESSTHRESHOLD = 0.99

### 
##########################
### adaptive PMCMC algorithm parameters
# If you want to try the adaptive PMCMC algorithm to compare the results:
# Number of x-particles:
NPMCMC = 500
# Number of iterations:
TPMCMC = 1000

### 
##########################
### Models and datasets options
# MODEL is a string, there has to be files named MODELx.py and MODELtheta.py
# in the models/ subfolder. Provided models are: 
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


### 
##########################
### Results options

# Specify time at which the theta-particles should be saved
# (useful when doing sequential analysis).
# The final time T = number of observations is automatically saved,
# no need to add it to the list.
#SAVINGTIMES = [250, 500, 750]
SAVINGTIMES = []
SMOOOTHING = False
FILTERING = False
PREDICTION = True
SMOOTHINGTIMES = []
STORESMOOTHINGTIME = 0

### Should the program report profiling? (Slows the program a little bit).
PROFILING = True 

# Specify results file name (without extensions),
# leave empty if an automatic name is preferred, based on the algorithm parameters.
RESULTSFILENAME = "tempBSMC"
# Replace already existing files; if False,
# a counter will be added to the result file name, so that no existing files are erased.
REPLACEFILE = True
# Use subfolders for each model and each dataset to organize the results.
# If False, use long names and store the results at the root of the results/ folder.
USESUBFOLDERS = False
# RESULTSFILETYPE could include "RData" and "cpickle".
# This list is used to save the results either in RData format 
# or in cPickle format or in both formats. RData is required to use the graph programs
# of subfolder rgraphs/.
#RESULTSFILETYPE = ["RData", "cpickle"]
RESULTSFILETYPE = ["RData"]
##########################



