###
# warning: this should be modified so that the init distribution
# of the parameters
# is used when it is specified (instead of the prior)
#
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
        array, float32, int32, cov, load, isinf, isnan
from numpy import max as numpymax
from numpy import min as numpymin
from numpy import sum as numpysum
from resampling import IndResample

def ESSfunction(weights):
    """
    Computes the ESS, given unnormalized weights.
    """
    norm_weights = weights / sum(weights)
    sqweights = power(norm_weights, 2)
    return 1 / sum(sqweights)

class SOPF:
    """
    Kitagawa's Self Organizing Particle Filter
    """
    def __init__(self, model, algorithmparameters, savingtimes = [], autoinit = True):
        self.modelx = model["modelx"]
        self.modeltheta = model["modeltheta"]
        self.observations = model["observations"]
        self.statedimension = self.modelx.xdimension
        self.obsdimension = self.modelx.ydimension
        self.N = algorithmparameters["N"]
        self.T = self.observations.shape[0]
        self.thetaparticles = zeros((self.modeltheta.parameterdimension, self.N))
        self.xparticles = zeros((self.statedimension, self.N))
        self.xweights = zeros(self.N)
        self.logxweights = zeros(self.N)
        self.constants = zeros(self.T)
        self.savingtimes = savingtimes
        self.savingtimes.append(self.T)
        self.allReducedParticles = []
        self.allCounts = []
        print "------------------"
        print "launching Kitagawa's SOPF, with algorithm parameters:"
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
    def first_step(self):
        """
        First step: generate N theta-particles from the prior, and then
        for each theta_i, simulate one x-particle from the initial distribution
        p(x_0 | theta_i)
        """
        self.thetaparticles[...] = self.modeltheta.priorgenerator(self.N)
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
            self.resample()
            if t in self.savingtimes or t == self.T - 1:
                reducedParticles, counts = self.reduceParticles(self.thetaparticles)
                self.allReducedParticles.append(reducedParticles)
                self.allCounts.append(counts)
    def reduceParticles(self, thetaparticles):
        """
        Since a lot of theta-particles are identical, we save
        a lot of memory by counting the number of occurrences of each particles
        and storing the particles with the number of counts. 
        This takes a long time but allows to store the results in a more
        efficient way, and it might be necessary if N is large (> 10**7).
        """
        print "reducing particles..."
        onedimension = list(thetaparticles[0, :])
        uniqueelements = list(set(onedimension))
        nbuniqueparticles = len(uniqueelements)
        print "# of unique particles", nbuniqueparticles
        proportions = zeros(nbuniqueparticles, dtype = int32)
        reducedparticles = zeros((nbuniqueparticles, thetaparticles.shape[0]))
        for indextheta in xrange(self.N):
            index = uniqueelements.index(thetaparticles[0, indextheta])
            proportions[index] += 1
            reducedparticles[index, :] = thetaparticles[:, indextheta]
        proportions = proportions / self.N
        return reducedparticles, proportions 
    def getResults(self):
        resultsDict = {"trueparameters" : self.modelx.parameters,\
                "N" : self.N, "T" : self.T, "nbparameters" : self.modeltheta.parameterdimension, \
                "allreducedparticles" : self.allReducedParticles, "observations": self.observations, \
                "allcounts" : self.allCounts, "savingtimes" : self.savingtimes}
        return resultsDict
            



