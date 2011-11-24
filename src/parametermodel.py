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
from numpy import min as numpymin
from numpy import max as numpymax
import os

class ParameterModel:
    def __init__(self, name, dimension):
        print 'creating parameter model "%s"' % name
        self.name = name
        self.parameterdimension = dimension 
        def update(hyperparameters, observations):
            return hyperparameters
        self.updateHyperParam = update
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
                parameters[index, :] = exp(transformedparameters[index, :]) /\
                        (1 + exp(transformedparameters[index, :]))
            else:
                parameters[index, :] = transformedparameters[index, :]
        return parameters











