#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32, zeros_like, savez, \
        newaxis, load, genfromtxt, loadtxt, savetxt
from scipy.stats import norm, truncnorm, gamma
import os

filename = "Data_Anisus.csv"
data = loadtxt(filename, usecols = (1, 2), delimiter = ",", dtype = 'str', skiprows = 1)
years = [x.replace('"', '') for x in data[:,0]]
population = [int(x) for x in data[:,1]]
print years
print population
#import rpy2.robjects as robjects
savetxt("anisus.R", population)
#r = robjects.r
#y = robjects.FloatVector(population)
#r(""" observations <- %s """ % y.r_repr())
#r(""" dump("observations", file = "anisus.R")""" )
#
#f = open("anisus.R", "r")
#txt = f.read()
#f.close()
#txt = txt.replace("\n", "")
#f = open("anisus.R", "w")
#f.write(txt)
#f.close()
#

