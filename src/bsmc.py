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
from numpy import random, power, sqrt, exp, zeros, zeros_like,\
        ones, mean, average, prod, log, sum, repeat, newaxis, \
        array, float32, int32, cov, load, isinf, isnan, zeros_like, \
        var, linalg, pi, dot, argmax, transpose, diag
from numpy import max as numpymax
from numpy import min as numpymin
from numpy import sum as numpysum
import scipy.weave as weave
from resampling import IndResample
from various import fastWeightedCov, ESSfunction

class BSMC:
    """
    Liu & West's SMC for parameter estimation, 
    using a random walk to move the parameters, 
    and a shrinkage parameter for the variance not 
    to explode too quickly.
    """
    def __init__(self, model, algorithmparameters, savingtimes = [], autoinit = True):
        # models... 
        self.modelx = model["modelx"]
        self.modeltheta = model["modeltheta"]
        self.observations = model["observations"]
        # parameters...
        self.N = algorithmparameters["N"]
        self.smooth = algorithmparameters["smooth"]
        self.ESSthreshold = algorithmparameters["ESSthreshold"]
        # useful quantities that we want to access easily 
        self.statedimension = self.modelx.xdimension
        self.obsdimension = self.modelx.ydimension
        self.T = self.observations.shape[0]
        self.thetaparticles = zeros((self.modeltheta.parameterdimension, self.N))
        self.transformedthetaparticles = zeros((self.modeltheta.parameterdimension, self.N))
        self.xparticles = zeros((self.statedimension, self.N))
        self.xweights = zeros(self.N)
        self.logxweights = zeros(self.N)
        self.constants = zeros(self.T)
        self.savingtimes = savingtimes
        self.savingtimes.append(self.T)
        self.hsq = power(self.smooth, 2)
        self.shrink = sqrt(1 - self.hsq)
        self.ESS = zeros(self.T)
        self.resamplingindices = []
        # number of already past saving times
        self.alreadystored = 0
        self.thetahistory = zeros((len(self.savingtimes), self.modeltheta.parameterdimension, self.N))
        self.weighthistory = zeros((len(self.savingtimes), self.N))
        print "------------------"
        print "launching Liu and West's SMC, with algorithm parameters:"
        for key, element in algorithmparameters.items():
            print key, ":", element
        print "------------------"
        if autoinit:
            self.first_step()
            self.next_steps()
    def resample(self):
        """
        Resample (x, theta) particles according to their weights.
        """
        indices = IndResample(self.xweights, self.N)
        self.xparticles[...] = self.xparticles[:, indices]
        self.thetaparticles[...] = self.thetaparticles[:, indices]
        self.transformedthetaparticles[...] = self.transformedthetaparticles[:, indices]
    def first_step(self):
        """
        First step: generate N theta-particles from the prior, and then
        for each theta_i, simulate one x-particle from the initial distribution
        p(x_0 | theta_i)
        """
        if self.modeltheta.hasInitDistribution:
            print "init distribution is specified, using it..."
            self.thetaparticles[...] = self.modeltheta.rinit(self.N)
            self.transformedthetaparticles[...] = self.modeltheta.transform(self.thetaparticles)
            for i in range(self.N):
                self.logxweights[i] = self.modeltheta.priorlogdensity(self.transformedthetaparticles[:, i]) - \
                        self.modeltheta.dinit(self.transformedthetaparticles[:, i])
        else:
            print "no init distribution is specified, using prior distribution instead..."
            self.thetaparticles[...] = self.modeltheta.priorgenerator(self.N)
            self.transformedthetaparticles[...] = self.modeltheta.transform(self.thetaparticles)
        for i in range(self.N):
            self.xparticles[:, i] = self.modelx.firstStateGenerator(self.thetaparticles[:, i], size = 1)
    def next_steps(self):
        """
        Perform all the iterations until time T == number of observations.
        """
        for t in range(0, self.T):
            print "time %i" % t
            TandWresults = self.modelx.transitionAndWeight(self.xparticles[newaxis, ...], \
                    self.observations[t], self.thetaparticles, t + 1)
            self.xparticles[...] = TandWresults["states"][0, ...]
            self.logxweights[:] = TandWresults["weights"][0, :]
            self.logxweights[isnan(self.logxweights)] = -(10**150)
            self.logxweights[isinf(self.logxweights)] = -(10**150)
            self.constants[t] = numpymax(self.logxweights)
            self.logxweights[:] -= self.constants[t]
            self.xweights[:] = exp(self.logxweights)
            covmean = self.computeCovarianceAndMean()
            m = (self.shrink) * self.transformedthetaparticles + \
                (1 - self.shrink) * transpose(covmean["mean"][newaxis])
            noise = transpose(random.multivariate_normal(repeat(0, self.modeltheta.parameterdimension), \
                self.hsq * covmean["cov"], size = self.N))
            self.transformedthetaparticles[...] = m + noise
            self.thetaparticles[...] = self.modeltheta.untransform(self.transformedthetaparticles)
            self.ESS[t] = ESSfunction(self.xweights[:])
            if self.ESS[t] < (self.ESSthreshold * self.N):
                self.resample()
                self.resamplingindices.append(t)
            if t in self.savingtimes or t == self.T - 1:
                print "saving particles at time %i" % t
                self.thetahistory[self.alreadystored, ...] = self.thetaparticles.copy()
                self.weighthistory[self.alreadystored, ...] = self.xweights.copy()
                self.alreadystored += 1
    def computeCovarianceAndMean(self):
        X = transpose(self.transformedthetaparticles)
        res = fastWeightedCov(X, self.xweights[:])
        res["cov"] += diag(zeros(self.modeltheta.parameterdimension) + \
            10**(-4)/self.modeltheta.parameterdimension)
        return res
    def getResults(self):
        resultsDict = {"trueparameters" : self.modelx.parameters,\
                "N" : self.N, "T" : self.T, "nbparameters" : self.modeltheta.parameterdimension, \
                "observations": self.observations, \
                "savingtimes" : self.savingtimes, \
                "thetahistory": self.thetahistory, \
                "weighthistory": self.weighthistory, \
                "ESS": self.ESS, \
                "resamplingindices": self.resamplingindices}
        return resultsDict
            




