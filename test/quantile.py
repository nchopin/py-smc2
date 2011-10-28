
#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32, zeros_like, newaxis, \
        arange
from scipy.stats import norm, truncnorm, gamma
import scipy.weave as weave
import os
import math



print norm.pdf(0)
print arange(10)
print repeat(1, 10)
#print norm.ppf(repeat(2.5 / 100, 10), loc = arange(10), scale = repeat(1, 10))
def pppp(location):
    return(norm.ppf(0.025, location, 1))
#print map(pppp, arange(10))

print norm.ppf(0.025, 2, 23)
print 2 + 23 * norm.ppf(0.025, 0, 1)
print norm.ppf(0.025)

