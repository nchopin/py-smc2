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
import os
import scipy.weave as weave
from numpy import random, array, sqrt, empty, double, log, exp, \
        sum, power, cumsum, zeros, zeros_like
from scipy.stats import norm

def IndResample(weights, nparticles):
    """ Indexed of resampled particles
    (deterministic scheme) """
    code = \
    """
    int j = 0;
    double csw = weights(0);
    for(int k = 0; k < N; k++)
    {
        while(csw < u)
        {
        j++;
        csw += weights(j);
        }
        Ind(k) = j;
        u = u + 1.;
    }
    """
    u = float(random.uniform(size = 1, low = 0, high = 1)[0])
    N = nparticles
    weights = double(weights * N / sum(weights))
    Ind = zeros(N, dtype='int')
    weave.inline(code,['u','N','weights','Ind'], type_converters=weave.converters.blitz)
    return Ind

def IndResample2D(weights, Nx, Ntheta):
    """ Indexed of resampled particles
    (deterministic scheme) """
    code = \
    """
    for(int i = 0; i < Ntheta; i++)
    {
        int j = 0;
        double csw = weights(0, i);
        for(int k = 0; k < Nx; k++)
        {
            while(csw < u(i))
            {
            j++;
            csw += weights(j, i);
            }
            Ind(k, i) = j;
            u(i) = u(i) + 1.;
        }
    }
    """
    u = random.uniform(size = Ntheta, low = 0, high = 1).astype(float)
    weights = double(weights * Nx / sum(weights, axis = 0))
    Ind = zeros((Nx, Ntheta), dtype='int')
    weave.inline(code,['u','Nx', 'Ntheta','weights','Ind'], type_converters=weave.converters.blitz)
    return Ind

def resample2D(array, weights, Nx, xdim, Ntheta):
    """ Indexed of resampled particles
    (deterministic scheme) """
    code = \
    """
    for(int i = 0; i < Ntheta; i++)
    {
        int j = 0;
        double csw = weights(0, i);
        for(int k = 0; k < Nx; k++)
        {
            while(csw < u(i))
            {
                j++;
                csw += weights(j, i);
            }
            for (int l = 0; l < xdim; l++){
                newarray(k, l, i) = array(j, l, i);
            }
            u(i) = u(i) + 1.;
        }
    }
    """
    u = random.uniform(size = Ntheta, low = 0, high = 1).astype(float)
    weights = double(weights * Nx / sum(weights, axis = 0))
    newarray = zeros_like(array)
    weave.inline(code,['u','Nx', 'Ntheta','weights', 'array', 'newarray', 'xdim'], type_converters=weave.converters.blitz)
    return newarray




