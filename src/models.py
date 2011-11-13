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
from numpy import min as numpymin
from numpy import max as numpymax
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
        self.predictionfunctionals = {}
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
        # load hidden states if it is a synthetic data set
        # (for plotting purpose)
        hiddenstatefile = filename.replace(".R", "-states.R")
        if os.path.exists(hiddenstatefile):
            self.model_states = genfromtxt(hiddenstatefile)
        else:
            self.model_states = "unknown"
    def setRLinearGaussian(self, string):
        self.RLinearGaussian = string
    def setRlikelihood(self, alist):
        self.Rlikelihood = alist
class ParameterModel:
    def __init__(self, name, dimension):
        print 'creating parameter model "%s"' % name
        self.name = name
        self.parameterdimension = dimension 
        #self.priorandtruevaluesspecified = priorandtruevaluesspecified
#        self.hasInitDistribution = False
        self.plottingInstructions = []
        def update(hyperparameters, observations):
            return hyperparameters
        self.updateHyperParam = update
        self.additionalPlots = ""
    def setRtruevalues(self, truevalues):
        self.truevalues = truevalues
    def setRprior(self, Rfunctionlist):
        self.Rprior = Rfunctionlist
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
        #print "\n cov matrix"
        #print proposalcovmatrix
        proposedparameters = zeros_like(currentvalue)
        if (nbparameters == 1):
            noise = transpose(random.normal(0, \
                proposalcovmatrix, size = currentvalue.shape[1]))
        else:
            noise = transpose(random.multivariate_normal(repeat(0, nbparameters), \
                proposalcovmatrix, size = currentvalue.shape[1]))
        if proposalkernel == "randomwalk":
            proposedparameters = currentvalue + noise
        elif proposalkernel == "independent":
            proposedparameters = proposalmean[:, newaxis] + noise
        return proposedparameters
#    def setInitDistribution(self, rfunction, dfunction):
#        self.hasInitDistribution = True
#        self.rinit = rfunction
#        self.dinit = dfunction
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
#                indTooBig = (transformedparameters[index, :] > 10)
#                indTooSmall = (transformedparameters[index, :] < -10)
#                indProblem = indTooBig | indTooSmall
#                trparam = zeros_like(transformedparameters[index, :])
#                trparam[indOK] = exp(transformedparameters[index, :][indOK]) / (1 + exp(transformedparameters[index, :][indOK]))
#                trparam[1 - indProblem] = exp(transformedparameters[index, 1 - indProblem]) / \
#                                          (1 + exp(transformedparameters[index, 1 - indProblem]))
#                trparam[indTooBig] = 1 - exp(-transformedparameters[index, indTooBig])
#                trparam[indTooSmall] = exp(transformedparameters[index, indTooSmall])
#                print "\ninside untransform"
#                print sum(trparam == 1)
#                print "the problem comes from parameters:"
#                print transformedparameters[index,:][trparam == 1]
#                print numpymin(trparam), numpymax(trparam)
                parameters[index, :] = exp(transformedparameters[index, :]) / (1 + exp(transformedparameters[index, :]))
                #print "range of the untransformed particles"
                #print numpymin(parameters[index, :]), numpymax(parameters[index, :])
            else:
                parameters[index, :] = transformedparameters[index, :]
        return parameters

