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
import os, os.path, sys, imp
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32, cov, load, isnan, isinf, newaxis
from numpy import max as numpymax
from numpy import min as numpymin
from numpy import sum as numpysum
from scipy.stats import norm
from resampling import IndResample
import time

class SingleSIR:
    """
    This class launches a SMC filter with Nx particles
    for a single theta parameter.
    The goal is to sample a x-trajectory, so the 
    x-trajectories are stored.
    """
    def __init__(self, Nx, theta, observations, modelx, \
            savingtimes = [], storall = False, \
            verbose = False, autoinit = True):
        self.modelx = modelx
        self.observations = observations
        self.statedimension = modelx.xdimension
        self.obsdimension = modelx.ydimension
        self.Nx = Nx
        self.T = observations.shape[0]
        self.excludedobservations = self.modelx.excludedobservations
        self.theta = theta
        self.savingtimes = savingtimes
        self.storall = storall
        if self.storall:
            self.savingtimes = range(T)
            self.lineage = zeros((Nx, self.T), dtype = int32)
        if len(self.savingtimes) > 0:
            self.alreadystored = 0
            if not (0 in self.savingtimes):
                self.savingtimes.append(0)
            self.savingtimes.sort()
            self.xhistory = zeros((Nx, self.statedimension, len(self.savingtimes)))
        self.xparticles = zeros((Nx, self.statedimension))
        self.xweights = zeros(Nx)
        self.logxweights = zeros(Nx)
        self.constants = zeros(self.T)
        self.verbose = verbose
        self.meanpath = zeros((self.T, self.statedimension))
        self.computingtimes = zeros(self.T)
        if self.verbose:
            print "Single SIRs with %i x-particles, %i observations" % \
                    (self.Nx, self.observations.shape[0])
        if autoinit:
            self.first_step()
            self.next_steps()
    def xresample(self, t):
        parentsindices = IndResample(self.xweights, self.Nx)
        self.xparticles[...] = self.xparticles[parentsindices, :]
        if self.storall:
            self.lineage[:, t] = parentsindices
    def first_step(self):
        self.xparticles[...] = self.modelx.firstStateGenerator(self.theta, size = self.Nx)
        if 0 in self.savingtimes:
            self.xhistory[..., 0] = self.xparticles
            self.alreadystored += 1
    def next_steps(self):
        last_tic = time.time()
        for t in range(self.T):
            excluded = t in self.excludedobservations
            if self.verbose:
                print "time %i" % t
                if excluded:
                    print "observations", self.observations[t,:], "set to be excluded"
            TandWresults = self.modelx.transitionAndWeight(self.xparticles[..., newaxis], \
                    self.observations[t,:], self.theta[:, newaxis], t + 1)
            self.xparticles[...] = TandWresults["states"][..., 0]
            if not(excluded):
                self.logxweights[...] = TandWresults["weights"][..., 0]
                self.logxweights[isnan(self.logxweights)] = -(10**150)
                self.logxweights[isinf(self.logxweights)] = -(10**150)
                self.constants[t] = numpymax(self.logxweights)
                self.logxweights[...] -= self.constants[t]
            else:
                self.logxweights = zeros(self.Nx)
                self.constants[t] = numpymax(self.logxweights)
            self.xweights[...] = exp(self.logxweights)
            self.xresample(t)
            if ((t+1) in self.savingtimes):
                self.xhistory[..., self.alreadystored] = self.xparticles
                self.alreadystored += 1
            self.meanpath[t,:] = mean(self.xparticles, axis = 0)
            new_tic = time.time()
            self.computingtimes[t] = new_tic - last_tic
            last_tic = new_tic
    def retrieveTrajectory(self, particleindex):
        """ 
        return a complete trajectory (starting from t = 0)
        so the length is T + 1
        """
        trajectory = zeros((self.T + 1, self.statedimension))
        trajectory[self.T, :] = self.xhistory[particleindex, :, self.T]
        parentindex = particleindex
        for t in reversed(xrange(self.T)):
            parentindex = self.lineage[parentindex, t]
            trajectory[t, :] = self.xhistory[parentindex, :, t]
        return trajectory
    def getResults(self):
        """
        Return a dictionary with vectors of interest.
        """
        model_obs = self.modelx.model_obs[0:self.T, :]
        resultsDict = {"parameters": self.theta, \
                "statedimension": self.statedimension, \
                "Nx" : self.Nx, "T": self.T, \
                "meanpath": self.meanpath, \
                "observations": model_obs, \
                "computingtimes": self.computingtimes}
        if self.storall:
            ntraj = 50
            trajectories = zeros((self.T + 1, self.statedimension, ntraj))
            for index in range(ntraj):
                trajectories[..., index] = self.retrieveTrajectory(index)
            resultsDict.update({"trajectories": trajectories})
        if len(self.savingtimes) > 0:
            resultsDict.update({"xhistory": self.xhistory, "savingtimes":
                self.savingtimes})
        return resultsDict

