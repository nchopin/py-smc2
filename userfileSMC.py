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

# If RANDOMSEED is True, a new seed is used at each run, otherwise it's fixed (for debugging purposes)
RANDOMSEED = True

MODEL = "locallevel"
T = 100
DATASET = "synthetic"

##########################
### Plot options
# Generate R file to plot results (does not require R)
GENERATERFILE = True
## Generate pdf figures using the R file (requires R and package "ggplot2")
PLOT = True

###
##########################
### SMC algorithm parameters
# NX: number of x-particles 
###
NX = 25000

### 
##########################
### Results options

TIMES = True
### Should the program report profiling? (Slows the program a little bit).
PROFILING = True 

# Specify results file name (without extensions),
# leave empty if an automatic name is preferred, based on the algorithm parameters.
RESULTSFILENAME = "temp"
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
RESULTSFILETYPE = ["RData"]
##########################




