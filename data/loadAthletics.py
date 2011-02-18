#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32, zeros_like, savez, \
        newaxis, load, genfromtxt, loadtxt, savetxt, vstack, hstack, transpose
from scipy.stats import norm, truncnorm, gamma
import os
import re

filename = "data/athletics.txt"
#data = loadtxt(filename, usecols = (2, 3), delimiter = ",", dtype = 'str', skiprows = 1)
data = loadtxt(filename, delimiter = r"\s+", dtype = 'str', skiprows = 6)
years = []
records = []
times = []
for i in range(data.shape[0]):
    splittedlist = [x for x in data[i].split(r" ") if len(x) > 0]
    years.append(splittedlist[len(splittedlist) - 1])
    records.append(splittedlist[1])

years = [x[(len(x) - 4):len(x)] for x in years]
for i, r in enumerate(records):
    #print r
    r = r.replace("+", "")
    r = re.split(r"[:]", r)
    time = float(r[0]) * 60 + float(r[1])
    times.append(time)
#print array(times)
#print array(years)
dataset = transpose(vstack([times, years]))
#print dataset.shape
#print dataset

uniqueyears = list(set(years))
uniqueyears.sort()
print uniqueyears
nbobsperyear = []
besttwo = zeros((len(uniqueyears), 2))
for ind, year in enumerate(uniqueyears):
    #print year
    timeperyear = [dataset[i,0] for i in range(dataset.shape[0]) if dataset[i,1] == year]
    timeperyear.sort()
    #print timeperyear
    besttwo[ind,:] = timeperyear[0:2]
    nbobsperyear.append(len(timeperyear))
#print dataset
print uniqueyears
print nbobsperyear
print besttwo



savetxt("athletics-best-two.R", besttwo)





