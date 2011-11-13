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
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, average, prod, log, sum, repeat, \
        array, zeros_like, newaxis
from numpy import sum as numpysum
from numpy import max as numpymax
from scipy.stats import norm, truncnorm, gamma
import scipy.weave as weave
import os
import math

random.seed(127)
allK = random.poisson(lam = 0, size = 1)

print len(allK)
allK[allK > 10000] = -1000
print allK
x = array([1, 2, 3, 5])
ind = x > 2.4
print ind
x[1 - ind] = 99
print x


a = array([True, False, False, False])
b = array([False, False, True, False])
print 1 - (a | b)




