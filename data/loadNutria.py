#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32, zeros_like, savez, \
        newaxis, load, genfromtxt, loadtxt, savetxt
from scipy.stats import norm, truncnorm, gamma
import os

filename = "Data_nutria.csv"
data = loadtxt(filename, usecols = (2, 3), delimiter = ",", dtype = 'str', skiprows = 1)
#print data
years = data[:,0]
population = [float(x) for x in data[:,1]]
print years
print population
savetxt("nutria.R", population)

#import rpy2.robjects as robjects
#r = robjects.r
#y = robjects.FloatVector(population)
#r(""" observations <- %s """ % y.r_repr())
#r(""" dump("observations", file = "nutria.R")""" )
#
#f = open("nutria.R", "r")
#txt = f.read()
#f.close()
#txt = txt.replace("\n", "")
#f = open("nutria.R", "w")
#f.write(txt)
#f.close()
#





