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
        transpose, newaxis, savetxt, genfromtxt, minimum, maximum
from scipy.stats import norm, truncnorm, gamma
import scipy.weave as weave
import os

class SSM:
    def __init__(self, name = "", xdimension = 1, ydimension = 1):
        print 'creating state space model "%s"' % name
        self.name = name
        self.xdimension = xdimension
        self.ydimension = ydimension
        self.excludedobservations = []
        self.functionals = {}
    def setFirstStateGenerator(self, function):
        self.firstStateGenerator = function
    def setVectorTransition(self, function):
        self.vectorTransition2D = function
    def setTransitionAndWeight(self, function):
        self.transitionAndWeight = function
    def setObservationGenerator(self, function):
        self.observationGenerator = function
    def generateData(self, size, parameters, savefilename):
        print "generating data"
        self.parameters = parameters
        self.model_states = zeros((size + 1, self.xdimension))
        self.model_obs = zeros((size, self.ydimension))
        self.model_states[0, :] = self.firstStateGenerator(parameters, 1)
        for t in xrange(0, size):
            previousstate = self.model_states[t, :]
            TaWresults = self.transitionAndWeight(previousstate[newaxis, :, newaxis], y = 0., parameters = parameters[:, newaxis], t = t)
            self.model_states[t + 1, :] = TaWresults["states"][0, :, 0]
        self.model_states = self.model_states[1:(size + 1), :]
        self.model_obs[...] = self.observationGenerator(self.model_states, parameters)
        savetxt(savefilename, self.model_obs)
        savefilenamestates = savefilename.replace(".R", "-states.R")
        savetxt(savefilenamestates, self.model_states)
    def loadData(self, filename):
        """
        """
        print "loading data from %s" % filename
        self.model_obs =  genfromtxt(filename)
        if len(self.model_obs.shape) == 1:
            self.model_obs = self.model_obs[:, newaxis]
        print "loaded: %i observations" % self.model_obs.shape[0]
        
class ParameterModel:
    def __init__(self, name, dimension, priorandtruevaluesspecified = False):
        print 'creating parameter model "%s"' % name
        self.name = name
        self.parameterdimension = dimension 
        self.priorandtruevaluesspecified = priorandtruevaluesspecified
        self.hasInitDistribution = False
        self.plottingInstructions = []
        def update(hyperparameters, observations):
            return hyperparameters
        self.updateHyperParam = update
        self.RpriorAvailable = False
        self.additionalPlots = ""
        self.truevaluesAvailable = False
    def setRtruevalues(self, truevalues):
        self.truevalues = truevalues
        self.truevaluesAvailable = True
    def setRprior(self, Rfunctionlist):
        self.Rfunctionlist = Rfunctionlist
        self.RpriorAvailable = True
    def setParameterNames(self, names):
        self.parameternames = names
    def setHyperparameters(self, hyperparameters):
        self.hyperparameters = hyperparameters
    def setUpdateHyperparametersWithData(self, function):
        self.updateHyperParam = function
    def setPriorlogdensity(self, function):
        self.priorlogdensityfunction = function
    def priorlogdensity(self, parameters):
        return self.priorlogdensityfunction(parameters, self.hyperparameters)
    def setPriorgenerator(self, function):
        self.priorgeneratorfunction = function
    def priorgenerator(self, size):
        return self.priorgeneratorfunction(size, self.hyperparameters)
    def proposal(self, currentvalue, proposalcovmatrix, hyperparameters = {}, \
            proposalkernel = "randomwalk", proposalmean = None):
        """ takes and returns transformed parameters """
        nbparameters = currentvalue.shape[0]
        proposedparameters = zeros_like(currentvalue)
        noise = transpose(random.multivariate_normal(repeat(0, nbparameters), \
                proposalcovmatrix, size = currentvalue.shape[1]))
        if proposalkernel == "randomwalk":
            proposedparameters = currentvalue + noise
        elif proposalkernel == "independent":
            proposedparameters = proposalmean[:, newaxis] + noise
        return proposedparameters
    def setInitDistribution(self, rfunction, dfunction):
        self.hasInitDistribution = True
        self.rinit = rfunction
        self.dinit = dfunction
    def setTransformation(self, transformations):
        self.whichAreLog = [transfo == "log" for transfo in transformations]
        self.whichAreLogit = [transfo == "logit" for transfo in transformations]
    def transform(self, parameters):
        transformedparameters = zeros_like(parameters)
        for index in range(parameters.shape[0]):
            if self.whichAreLog[index]:
                transformedparameters[index, :] = log(parameters[index, :])
            elif self.whichAreLogit[index]:
                transformedparameters[index, :] = log(parameters[index, :]/(1 - parameters[index, :]))
            else:
                transformedparameters[index, :] = parameters[index, :]
        return transformedparameters
    def untransform(self, transformedparameters):
        parameters = zeros_like(transformedparameters)
        for index in range(transformedparameters.shape[0]):
            if self.whichAreLog[index]:
                parameters[index, :] = exp(transformedparameters[index, :])
            elif self.whichAreLogit[index]:
                parameters[index, :] = maximum(10**(-10), minimum((1 - 10**(-10)), exp(transformedparameters[index, :]) / (1 + exp(transformedparameters[index, :]))))
            else:
                parameters[index, :] = transformedparameters[index, :]
        return parameters

