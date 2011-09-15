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
import os, os.path
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32, cov, load, isnan, isinf, newaxis
from numpy import max as numpymax
from numpy import min as numpymin
from numpy import sum as numpysum
from scipy.stats import norm
from resampling import resample2D

def ESSfunction(weights):
    """
    Computes the ESS, given unnormalized weights.
    """
    norm_weights = weights / sum(weights)
    sqweights = power(norm_weights, 2)
    return 1 / sum(sqweights)

class ParallelSIRs:
    """
    This class launches a SMC filter with Nx particles
    for each theta in the matrix thetaparticles.
    It computes the log likelihood of the data
    in order to use it in a PMCMC acceptance ratio.
    Resampling is made at each iteration, and the 
    particles are moved according to the hidden state
    dynamic.
    """
    def __init__(self, Nx, thetaparticles, observations, modelx, \
            savepath = False, verbose = False, saveproposals = False):
        self.modelx = modelx
        self.observations = observations
        self.statedimension = modelx.xdimension
        self.obsdimension = modelx.ydimension
        self.Nx = Nx
        self.T = observations.shape[0]
        self.excludedobservations = self.modelx.excludedobservations
        self.thetaparticles = thetaparticles
        self.Ntheta = thetaparticles.shape[1]
        self.xparticles = zeros((self.Nx, self.statedimension, self.Ntheta))
        self.xweights = zeros((self.Nx, self.Ntheta))
        self.logxweights = zeros((self.Nx, self.Ntheta))
        self.constants = zeros((self.T, self.Ntheta))
        self.totalLogLike = zeros(self.Ntheta)
        self.verbose = verbose
        self.saveproposals = saveproposals
        self.savepath = savepath
        if self.saveproposals:
            self.allxparticles = zeros((self.T, self.Nx, self.statedimension, self.Ntheta))
            self.allproposals = zeros((self.T, self.Nx, self.statedimension, self.Ntheta))
        if self.savepath:
            self.paths = zeros((self.T, self.statedimension, self.Ntheta))
        if self.verbose:
            print "Parallel SIRs with %i theta-particles, %i x-particles, %i observations" % \
                    (self.Ntheta, self.Nx, self.observations.shape[0])
    def xresample(self):
        self.xparticles[...] = resample2D(self.xparticles, self.xweights, self.Nx, self.statedimension, self.Ntheta)
    def first_step(self):
        for i in range(self.Ntheta):
            self.xparticles[:, :, i] = self.modelx.firstStateGenerator(self.thetaparticles[:, i], size = self.Nx)
    def next_steps(self):
        for t in range(self.T):
            excluded = t in self.excludedobservations
            if self.verbose:
                print "time %i" % t
                if excluded:
                    print "observations", self.observations[t,:], "set to be excluded"
            TandWresults = self.modelx.transitionAndWeight(self.xparticles, \
                    self.observations[t,:], self.thetaparticles, t + 1)
            self.xparticles[...] = TandWresults["states"]
            if not(excluded):
                self.logxweights[...] = TandWresults["weights"]
                self.logxweights[isnan(self.logxweights)] = -(10**150)
                self.logxweights[isinf(self.logxweights)] = -(10**150)
                self.constants[t, :] = numpymax(self.logxweights, axis = 0)
                self.logxweights[...] -= self.constants[t, :]
            else:
                self.logxweights = zeros((self.Nx, self.Ntheta))
            self.xweights[...] = exp(self.logxweights)
            if self.saveproposals:
                self.allproposals[t, ...] = self.xparticles.copy()
            if self.savepath:
                for xdim in range(self.xparticles.shape[1]):
                    self.paths[t, xdim, :] = average(self.xparticles[:, xdim, :], weights = self.xweights, axis = 0)
            logLike = log(mean(self.xweights, axis = 0))
            logLike[isnan(logLike)] = -(10**150)
            logLike[isinf(logLike)] = -(10**150)
            self.totalLogLike += logLike
            if not(excluded):
                self.xresample()
            if self.saveproposals:
                self.allxparticles[t, ...] = self.xparticles.copy()
    def getTotalLogLike(self):
        csts = numpysum(self.constants, axis = 0)
        csts[isnan(csts)] = -(10**150)
        csts[isinf(csts)] = -(10**150)
        return self.totalLogLike + csts




