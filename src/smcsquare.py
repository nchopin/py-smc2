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
        array, float32, int32, cov, isnan, zeros_like, \
        var, isinf, linalg, pi, dot, argmax, transpose, diag
from numpy import max as numpymax
from numpy import min as numpymin
from numpy import sum as numpysum
from scipy.stats import norm
from resampling import IndResample, resample2D
from parallelSIRs import *

class SMCsquare:
    """
    """
    def __init__(self, model, algorithmparameters, \
            dynamicNx = True, savingtimes = [], autoinit = True):
        ## essential things
        # get the model:
        self.modelx = model["modelx"]
        self.modeltheta = model["modeltheta"]
        self.observations = model["observations"]
        self.statedimension = self.modelx.xdimension
        self.obsdimension = self.modelx.ydimension
        # get the basic algorithmic parameters
        self.AP = algorithmparameters
        self.Nx = algorithmparameters["Nx"]
        self.Ntheta = algorithmparameters["Ntheta"]
        self.T = self.observations.shape[0]
        self.excludedobservations = self.modelx.excludedobservations
        # initialize huge matrices and vectors
        self.thetaparticles = zeros((self.modeltheta.parameterdimension, self.Ntheta))
        self.transformedthetaparticles = zeros((self.modeltheta.parameterdimension, self.Ntheta))
        self.thetalogweights = zeros((self.T, self.Ntheta))
        self.xparticles = zeros((self.Nx, self.statedimension, self.Ntheta))
        self.xweights = zeros((self.Nx, self.Ntheta))
        self.logxweights = zeros((self.Nx, self.Ntheta))
        self.constants = zeros(self.Ntheta)
        self.logLike = zeros(self.Ntheta)
        self.totalLogLike = zeros(self.Ntheta)
        self.evidences = zeros(self.T)
        ## Filtering and Smoothing
        self.filteringEnable = algorithmparameters["filtering"]
        self.filtered = {}
        if self.filteringEnable:
            for functionalnames in self.modelx.functionals.keys():
                self.filtered[functionalnames] = zeros(self.T)
        self.smoothingEnable = algorithmparameters["smoothing"]
        self.smoothedmeans = {}
        self.smoothedvalues= {}
        if self.smoothingEnable:
            self.smoothingtimes = algorithmparameters["smoothingtimes"]
            self.storesmoothingtime = algorithmparameters["storesmoothingtime"]
        ## other things:
        # store ESS at each time
        self.ESS = zeros(self.T)
        # store iterations where resample-moves are performed
        self.resamplingindices = []
        # store the acceptance ratios of each move step
        self.acceptratios = []
        # store all Nx 
        # (in case it is automatically increasing along the iterations)
        self.Nxlist = [self.Nx]
        # store iterations at each Nx is increased
        self.increaseindices = [0]
        # iterations at which all the theta-particles are saved
        self.savingtimes = savingtimes
        self.savingtimes.append(self.T)
        # initialize matrices to store the weighted theta-particles
        self.thetahistory = zeros((len(self.savingtimes), self.modeltheta.parameterdimension, self.Ntheta))
        self.weighthistory = zeros((len(self.savingtimes), self.Ntheta))
        # number of already past saving times
        self.alreadystored = 0
        print "------------------"
        print "launching SMC^2, with algorithm parameters:"
        for key, element in algorithmparameters.items():
            print key, ":", element
        print "------------------"
        if autoinit:
            self.first_step()
            self.next_steps()

    def increaseParticlesNb(self, t):
        """
        Double the number of x-particles.
        """
        print "increasing Nx: from %i" % self.Nx
        self.Nx = 2 * self.Nx
        print "to %i" % self.Nx
        biggerSIRs = ParallelSIRs(self.Nx, self.thetaparticles, self.observations[0:(t+1),:], self.modelx)
        biggerSIRs.first_step()
        biggerSIRs.next_steps()
        biggerTotalLogLike = biggerSIRs.getTotalLogLike()
        self.thetalogweights[t, :] = self.thetalogweights[t, :] +  biggerTotalLogLike - self.totalLogLike
        self.totalLogLike = biggerTotalLogLike.copy()
        self.xparticles = biggerSIRs.xparticles.copy()
        self.xweights = zeros_like(biggerSIRs.xweights)
        self.logxweights = zeros_like(biggerSIRs.xweights)
        self.Nxlist.append(self.Nx)
        self.increaseindices.append(t)

    def xresample(self):
        """
        Resample all the x-particles according to their weights.
        """
        self.xparticles[...] = resample2D(self.xparticles, self.xweights, self.Nx, self.statedimension, self.Ntheta)
    def thetaresample(self, t):
        """
        Resample the theta-particles according to their weights.
        """
        indices = IndResample(exp(self.thetalogweights[t, :]), self.Ntheta)
        self.thetaparticles[...] = self.thetaparticles[:, indices]
        self.transformedthetaparticles[...] = self.transformedthetaparticles[:, indices]
        self.xparticles[...] = self.xparticles[:, :, indices]
        self.totalLogLike[:] = self.totalLogLike[indices]
        self.thetalogweights[t, :] = 0.

    def PMCMCstep(self, t):
        """
        Perform a PMMH move step on each theta-particle.
        """
        transformedthetastar = self.modeltheta.proposal(self.transformedthetaparticles, \
                self.proposalcovmatrix, hyperparameters = self.modeltheta.hyperparameters,\
                proposalmean = self.proposalmean, proposalkernel = self.AP["proposalkernel"])
        thetastar = self.modeltheta.untransform(transformedthetastar)
        proposedSIRs = ParallelSIRs(self.Nx, thetastar, self.observations[0:(t+1)], self.modelx)
        proposedSIRs.first_step()
        proposedSIRs.next_steps()
        proposedTotalLogLike = proposedSIRs.getTotalLogLike()
        acceptations = zeros(self.Ntheta)
        if self.AP["proposalkernel"] == "randomwalk":
            for i in range(self.Ntheta):
                proposedlogomega = proposedTotalLogLike[i] + self.modeltheta.priorlogdensity(transformedthetastar[:, i])
                currentlogomega = self.totalLogLike[i] + self.modeltheta.priorlogdensity(self.transformedthetaparticles[:, i])
                acceptations[i]  = (log(random.uniform(size = 1)) < (proposedlogomega - currentlogomega))
                if acceptations[i]:
                    self.transformedthetaparticles[:, i] = transformedthetastar[:, i].copy()
                    self.thetaparticles[:, i] = thetastar[:, i].copy()
                    self.totalLogLike[i] = proposedTotalLogLike[i].copy()
                    self.xparticles[:, :, i] = proposedSIRs.xparticles[:, :, i].copy()
        elif self.AP["proposalkernel"] == "independent":
            invSigma = linalg.inv(self.proposalcovmatrix)
            def multinorm_logpdf(x):
                centeredx = (x - self.proposalmean)
                return -0.5 * dot(dot(centeredx, invSigma), centeredx)
            for i in range(self.Ntheta):
                proposedlogomega = proposedTotalLogLike[i] + self.modeltheta.priorlogdensity(transformedthetastar[:, i])
                proposedlogomega -= multinorm_logpdf(transformedthetastar[:, i])
                currentlogomega = self.totalLogLike[i] + self.modeltheta.priorlogdensity(self.transformedthetaparticles[:, i])
                currentlogomega -= multinorm_logpdf(self.transformedthetaparticles[:, i])
                acceptations[i]  = (log(random.uniform(size = 1)) < (proposedlogomega - currentlogomega))
                if acceptations[i]:
                    self.transformedthetaparticles[:, i] = transformedthetastar[:, i].copy()
                    self.thetaparticles[:, i] = thetastar[:, i].copy()
                    self.totalLogLike[i] = proposedTotalLogLike[i].copy()
                    self.xparticles[:, :, i] = proposedSIRs.xparticles[:, :, i].copy()
        acceptrate = sum(acceptations) / self.Ntheta
        self.acceptratios.append(acceptrate)
        print "acceptance rate: %.3f" % (acceptrate)

    def first_step(self):
        """
        First step: generate Ntheta theta-particles from the prior, and then
        for each theta_i, simulate Nx x-particles from the initial distribution
        p(x_0 | theta_i)
        """
        if self.modeltheta.hasInitDistribution:
            print "init distribution is specified, using it..."
            self.thetaparticles[...] = self.modeltheta.rinit(self.Ntheta)
            self.transformedthetaparticles[...] = self.modeltheta.transform(self.thetaparticles)
            for i in range(self.Ntheta):
                self.thetalogweights[0, i] = self.modeltheta.priorlogdensity(self.transformedthetaparticles[:, i]) - \
                        self.modeltheta.dinit(self.transformedthetaparticles[:, i])
        else:
            print "no init distribution is specified, using prior distribution instead..."
            self.thetaparticles[...] = self.modeltheta.priorgenerator(self.Ntheta)
            self.transformedthetaparticles[...] = self.modeltheta.transform(self.thetaparticles)
        for i in range(self.Ntheta):
            self.xparticles[:, :, i] = self.modelx.firstStateGenerator(self.thetaparticles[:, i], size = self.Nx)

    def next_steps(self):
        """
        Perform all the iterations until time T == number of observations.
        """
        for t in range(self.T):
            excluded = t in self.excludedobservations
            print "time %i" % t
            if excluded:
                print "observations", self.observations[t,:], "set to be excluded"
            TandWresults = self.modelx.transitionAndWeight(self.xparticles, self.observations[t], self.thetaparticles, t)
            self.xparticles[...] = TandWresults["states"]
            if not(excluded):
                self.logxweights[...] = TandWresults["weights"]
                # in case the measure function returns nans or infs, set the weigths very low
                self.logxweights[isnan(self.logxweights)] = -(10**150)
                self.logxweights[isinf(self.logxweights)] = -(10**150)
                self.constants[:] = numpymax(self.logxweights, axis = 0)
                self.logxweights[...] -= self.constants[:]
            else:
                self.logxweights = zeros((self.Nx, self.Ntheta))
                self.constants[:] = numpymax(self.logxweights, axis = 0)
            self.xweights[...] = exp(self.logxweights)
            self.logLike[:] = log(mean(self.xweights, axis = 0)) + self.constants[:]

            if t > 0:
                self.evidences[t] = self.getEvidence(self.thetalogweights[t-1, :], self.logLike)
                self.totalLogLike[:] += self.logLike[:]
                self.thetalogweights[t, :] = self.thetalogweights[t-1, :] + self.logLike[:]
            else:
                self.evidences[t] = self.getEvidence(self.thetalogweights[t, :], self.logLike)
                self.totalLogLike[:] += self.logLike[:]
                self.thetalogweights[t, :] = self.thetalogweights[t, :] + self.logLike[:]
            self.thetalogweights[t, :] -= max(self.thetalogweights[t, :])
            self.xresample()
            self.ESS[t] = self.ESSfunction(exp(self.thetalogweights[t, :]))
            if self.AP["dynamicNx"]:
                print "ESS: %.3f, Nx: %i" % (self.ESS[t], self.Nx)
            else:
                print "ESS: %.3f" % self.ESS[t]
            if self.ESS[t] < (self.AP["ESSthreshold"] * self.Ntheta):
                print "resample step... "
                covdict = self.computeCovarianceAndMean(t)
                if self.AP["proposalkernel"] == "randomwalk":
                    self.proposalcovmatrix = self.AP["rwvariance"] * covdict["cov"]
                    self.proposalmean = None
                elif self.AP["proposalkernel"] == "independent":
                    self.proposalcovmatrix = covdict["cov"]
                    self.proposalmean = covdict["mean"]
                self.thetaresample(t)
                self.resamplingindices.append(t)
                print "move step... "
                for move in range(self.AP["nbmoves"]):
                    self.PMCMCstep(t)
                    if self.acceptratios[-1] < self.AP["dynamicNxThreshold"] and self.Nx < (self.AP["NxLimit"] / 2) \
                            and self.AP["dynamicNx"]:
                        print "increasing the number of x-particles... "
                        self.increaseParticlesNb(t)
            """ filtering and smoothing """
            if self.filteringEnable:
                self.filtering(t)
            if self.smoothingEnable and (t in self.smoothingtimes):
                self.smoothing(t)
            if t in self.savingtimes or t == self.T - 1:
                print "saving particles at time %i" % t
                self.thetahistory[self.alreadystored, ...] = self.thetaparticles.copy()
                self.weighthistory[self.alreadystored, ...] = exp(self.thetalogweights[t, :])
                self.alreadystored += 1
    def filtering(self, t):
        for key in self.filtered.keys():
            tempmatrix = self.modelx.functionals[key](self.xparticles, self.thetaparticles, t)
            tempmatrix = mean(tempmatrix, axis = 0)
            self.filtered[key][t] = average(tempmatrix, weights = exp(self.thetalogweights[t,:]))
    def smoothing(self, t):
        print "smoothing time"
        from singleSIR import SingleSIR
        smoothedx = zeros((t+1, self.statedimension, self.Ntheta))
        randomindices = random.randint(low = 0, high = self.Nx, size = self.Ntheta)
        for indextheta in range(self.Ntheta):
            singleSIR = SingleSIR(self.Nx, self.thetaparticles[:, indextheta], self.observations[0:(t+1),:], self.modelx)
            onepath = singleSIR.retrieveTrajectory(randomindices[indextheta])
            smoothedx[:, :, indextheta] = onepath[1:(t+2), :]
        for key in self.modelx.functionals.keys():
            smoothkey = key + "%i" % (t + 1)
            self.smoothedmeans[smoothkey] = zeros(t + 1)
            self.smoothedvalues[smoothkey] = zeros(self.Ntheta)
            for subt in range(t + 1):
                smoothedxt = smoothedx[subt, ...]
                tempmatrix = self.modelx.functionals[key](smoothedx[newaxis, subt, ...], self.thetaparticles, subt)
                tempmean = mean(tempmatrix, axis = 0)
                self.smoothedmeans[smoothkey][subt] = average(tempmean, weights = exp(self.thetalogweights[t,:]))
                if subt == self.storesmoothingtime:
                    self.smoothedvalues[smoothkey][:] = tempmean
#            print smoothkey
#            print "smoothing", self.smoothedmeans[smoothkey][0:(t+1),...]
#            print "filtering", self.filtered[key][0:(t+1),...]
        print "smoothing done!"

    def computeCovarianceAndMean(self, t):
        X = transpose(self.transformedthetaparticles)
        w = exp(self.thetalogweights[t, :])
        w = w / numpysum(w)
        weightedmean = average(X, weights = w, axis = 0)
        diagw = diag(w)
        part1 = dot(transpose(X), dot(diagw, X))
        Xtw = dot(transpose(X), w[:, newaxis])
        part2 = dot(Xtw, transpose(Xtw))
        numerator = part1 - part2
        denominator = 1 - numpysum(w**2)
        weightedcovariance = numerator / denominator
        # increase a little bit the diagonal to prevent degeneracy effects
        weightedcovariance += diag(zeros(self.modeltheta.parameterdimension) + 10**(-4)/self.modeltheta.parameterdimension)
        return {"mean": weightedmean, "cov": weightedcovariance}

    def ESSfunction(self, weights):
        """
        Compute the ESS, given unnormalized weights.
        """
        norm_weights = weights / sum(weights)
        sqweights = power(norm_weights, 2)
        return 1 / sum(sqweights)

    def getEvidence(self, thetalogweights, loglike):
        """
        Return the evidence at a given time
        """
        return average(exp(loglike), weights = exp(thetalogweights))

    def getResults(self):
        """
        Return a dictionary with vectors of interest.
        """
        resultsDict = {"trueparameters": self.modelx.parameters, \
                "nbparameters" : self.modeltheta.parameterdimension, \
                "Ntheta" : self.Ntheta, "Nx" : self.Nx, "T": self.T, \
                "thetahistory": self.thetahistory, "weighthistory": self.weighthistory, \
                "savingtimes": self.savingtimes, "ESS": self.ESS, "observations": self.observations, \
                "acceptratios": self.acceptratios, "Nxlist": self.Nxlist, \
                "increaseindices": self.increaseindices, "resamplingindices": self.resamplingindices, \
                "evidences": self.evidences, "filtered": self.filtered, "smoothedmeans": self.smoothedmeans, \
                "smoothedvalues": self.smoothedvalues}
        return resultsDict



