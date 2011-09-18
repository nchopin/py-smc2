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
import os, os.path, sys
from numpy import random, power, sqrt, exp, zeros, zeros_like,\
        ones, mean, average, prod, log, sum, repeat, newaxis, \
        array, float32, int32, cov, load, isinf, isnan, zeros_like, \
        var, linalg, pi, dot, argmax, transpose, diag
from numpy import max as numpymax
from numpy import min as numpymin
from numpy import sum as numpysum
import scipy.weave as weave

def ESSfunction(weights):
    """
    Computes the ESS, given unnormalized weights.
    """
    norm_weights = weights / sum(weights)
    sqweights = power(norm_weights, 2)
    return 1 / sum(sqweights)

def fastWeightedCov(X, unnormalizedw):
    weights = unnormalizedw / numpysum(unnormalizedw)
    Xbar = average(X, weights = weights, axis = 0)
    code = \
    """
    int row,col;
    for (row = 0; row < d(0); row++)
    {
        for(col = 0; col < d(0); col++)
        {
          for (int index = 0; index < N(0); index ++){
            covariance(row, col) += weights(index) * (X(index, row) - Xbar(row)) * (X(index, col) - Xbar(col));
          }
        }
    }
    """
    d = X.shape[1]
    covariance = zeros((d, d))
    d = array([d])
    N = array([X.shape[0]])
    weave.inline(code,['covariance', 'd', 'N', 'Xbar', 'X', 'weights'], \
        type_converters=weave.converters.blitz, libraries = ["m"])
    weightedcovariance = covariance / (1 - numpysum(power(weights, 2)))
    return {"mean": Xbar, "cov": weightedcovariance}

def progressbar(ratio, text=None, ticks=50):
   progress = int(ticks * ratio)
   s = '%.1f%%' % (100.0 * ratio)
   length = len(s)
   if progress > ticks / 2 - length:
       sys.stdout.write('\r[' + int(ticks / 2 - length) * '-' + s
                        + int(progress - ticks / 2) * '-' + int(min(ticks -
                            progress, ticks / 2)) * ' ' + ']')
   else:
       sys.stdout.write('\r[' + progress * '-' + int(ticks / 2 - length - progress) * ' ' + s
                        + int(ticks / 2) * ' ' + ']')
   if not text is None: sys.stdout.write(text)
   sys.stdout.flush()

