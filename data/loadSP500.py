#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32, zeros_like, savez, \
        newaxis, load, genfromtxt, savetxt
from scipy.stats import norm, truncnorm, gamma
import os

filename = "SP500recent.csv"

data = genfromtxt(filename, delimiter = ",")
print type(data)
print data.shape

SP500 = array(data[0:(data.shape[0]),data.shape[1]-1])
SP500 = SP500[::-1]
SP500 = SP500[0:(SP500.shape[0] - 1)]
print SP500.shape
#rescale
SP500 = 10**(5/2)* SP500
savetxt(filename.replace(".csv", ".R"), SP500)
#### Draw the data
import rpy2.robjects as robjects
r = robjects.r
y = robjects.FloatVector(SP500)
r(""" observations <- %s """ % y.r_repr())
#r(""" dump("observations", file = "littleSP500.R")""" )
r("""pdf(file = "%s")""" % os.path.join(filename.replace(".csv", ".pdf")))
r.plot(y, ylab = "yields", xlab = "days", type = "l", col = "lightgreen")
r("dev.off()")
#raw_input("appuyez sur une touche pour fermer la fenetre")



