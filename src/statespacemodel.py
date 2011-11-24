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
from numpy import random, exp, zeros, \
        ones, mean, log, repeat, array, zeros_like, \
        transpose, newaxis, savetxt, genfromtxt, minimum, maximum, \
        argsort, var, cumsum, searchsorted, average
from scipy.stats import norm, truncnorm, gamma
import os

class SSM:
    def __init__(self, name = "", xdimension = 1, ydimension = 1):
        print 'creating state space model "%s"' % name
        self.name = name
        self.xdimension = xdimension
        self.ydimension = ydimension
        self.excludedobservations = []
        self.filteringlist = []
        self.predictionlist = []
    def setParameters(self, parameters):
        self.parameters = parameters
    def setFirstStateGenerator(self, function):
        self.firstStateGenerator = function
    def setVectorTransition(self, function):
        self.vectorTransition2D = function
    def setTransitionAndWeight(self, function):
        self.transitionAndWeight = function
    def setObservationGenerator(self, function):
        self.observationGenerator = function
    def addStateFiltering(self):
        for idim in range(self.xdimension):
            func_str = \
            """
def filttmp%(idim)i(xparticles, thetaweights, thetaparticles, t):
    meanpertheta = mean(xparticles[:, %(idim)i, :], axis = 0)
    xmean = average(meanpertheta, weights = thetaweights)
    return xmean
            """ % {"idim": idim}
            exec(func_str)
            f = locals()["filttmp%i" % idim]
            self.filteringlist.append({"function": f, "dimension": 1, "name": "state%i" % (idim + 1)})
    def addStatePrediction(self):
        for idim in range(self.xdimension):
            func_str = \
            """
def predtmp%(idim)i(xparticles, thetaweights, thetaparticles, t):
    Nx = xparticles.shape[0]
    Ntheta = xparticles.shape[2]
    predictedstate = zeros(Nx * Ntheta)
    weight = zeros(Nx * Ntheta)
    for j in range(Ntheta):
        predictedstate[(Nx * j):(Nx * (j+1))] = xparticles[..., %(idim)i, j]
        weight[(Nx * j):(Nx * (j+1))] = repeat(thetaweights[j], repeats = Nx)
    weight = weight / sum(weight)
    xmean = average(predictedstate, weights = weight)
    ind = argsort(predictedstate)
    predictedstate = predictedstate[ind]
    weight = weight[ind]
    cumweight = cumsum(weight)
    quantile5 = predictedstate[searchsorted(cumweight, 0.05)]
    quantile95 = predictedstate[searchsorted(cumweight, 0.95)]
    return array([xmean, quantile5, quantile95])
            """ % {"idim": idim}
            exec(func_str)
            f = locals()["predtmp%i" % idim]
            self.predictionlist.append({"function": f, "dimension": 3, "name": "state%i" % (idim + 1)})
    def addObsPrediction(self):
        def predictionObservations(xparticles, thetaweights, thetaparticles, t):
            Nx = xparticles.shape[0]
            Ntheta = xparticles.shape[2]
            result = zeros(3)
            observations = zeros(Nx * Ntheta)
            weightobs = zeros(Nx * Ntheta)
            for j in range(Ntheta):
                observations[(Nx * j):(Nx * (j+1))] = \
                        self.observationGenerator(xparticles[..., j], thetaparticles[:, j]).reshape(Nx)
                weightobs[(Nx * j):(Nx * (j+1))] = repeat(thetaweights[j], repeats = Nx)
            weightobs = weightobs / sum(weightobs)
            obsmean = average(observations, weights = weightobs)
            ind = argsort(observations)
            observations = observations[ind]
            weightobs = weightobs[ind]
            cumweightobs = cumsum(weightobs)
            quantile5 = observations[searchsorted(cumweightobs, 0.05)]
            quantile95 = observations[searchsorted(cumweightobs, 0.95)]
            return array([obsmean, quantile5, quantile95])
        self.predictionlist.append({"function": predictionObservations, "dimension": 3, "name": "observations"})
    def addPredictionList(self, l):
        self.predictionlist.extend(l)
    def generateData(self, size, parameters, savefilename):
        print "generating data"
        self.parameters = parameters
        self.model_states = zeros((size + 1, self.xdimension))
        self.model_obs = zeros((size, self.ydimension))
        self.model_states[0, :] = self.firstStateGenerator(parameters, 1)
        for t in xrange(0, size):
            previousstate = self.model_states[t, :]
            TaWresults = self.transitionAndWeight(previousstate[newaxis, :, newaxis], y = 0., \
                    parameters = parameters[:, newaxis], t = t)
            self.model_states[t + 1, :] = TaWresults["states"][0, :, 0]
        self.model_states = self.model_states[1:(size + 1), :]
        self.model_obs[...] = self.observationGenerator(self.model_states, parameters)
        savetxt(savefilename, self.model_obs)
        savefilenamestates = savefilename.replace(".R", "-states.R")
        savetxt(savefilenamestates, self.model_states)
    def loadData(self, filename):
        print "loading data from %s" % filename
        self.model_obs =  genfromtxt(filename)
        if len(self.model_obs.shape) == 1:
            self.model_obs = self.model_obs[:, newaxis]
        print "loaded: %i observations" % self.model_obs.shape[0]
        hiddenstatefile = filename.replace(".R", "-states.R")
        if os.path.exists(hiddenstatefile):
            self.model_states = genfromtxt(hiddenstatefile)

