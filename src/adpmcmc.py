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
        newaxis, diag, savez
from scipy.stats import norm
from resampling import IndResample, IndResample2D
from parallelSIRs import ParallelSIRs
from various import progressbar

class AdaptivePMCMC:
    def __init__(self, model, algorithmparameters, autoinit = True):
        # get the model:
        self.modelx = model["modelx"]
        self.modeltheta = model["modeltheta"]
        self.observations = model["observations"]
        self.Nx = algorithmparameters["N"]
        self.nbiterations = algorithmparameters["nbiterations"]
        #
        self.theta = self.modeltheta.priorgenerator(1)[:,0]
        self.thetadim = self.theta.shape[0]
        self.transformedtheta = self.modeltheta.transform(self.theta[:, newaxis])[:,0]
        self.chain = zeros((self.thetadim, self.nbiterations))
        self.transformedchain = zeros((self.thetadim, self.nbiterations))
        self.chain[:, 0] = self.theta.copy()
        self.transformedchain[:, 0] = self.transformedtheta.copy()
        #
        self.acceptations = zeros(self.nbiterations)
        self.fixedproposalcovmatrix = diag(repeat(0.01 / self.thetadim, self.thetadim))
        self.adaptiveproposalcovmatrix = diag(repeat(5.6644 / self.thetadim, self.thetadim))
        self.proposalmixtureweights = [0.5, 0.5]
        currentSIRs = ParallelSIRs(self.Nx, self.theta[:, newaxis], self.observations, self.modelx)
        currentSIRs.first_step()
        currentSIRs.next_steps()
        self.totalLogLike = currentSIRs.getTotalLogLike()[0]
        self.loglikechain = zeros(self.nbiterations)
        self.loglikechain[0] = self.totalLogLike
        print "------------------"
        print "launching adaptive PMCMC, with algorithm parameters:"
        for key, element in algorithmparameters.items():
            print key, ":", element
        print "------------------"
        if autoinit:
            self.allsteps()
    def PMCMCstep(self, iteration):
        #progressbar(iteration / (self.nbiterations - 1))
        progressbar(iteration / (self.nbiterations - 1), text = " acceptance rate: %.3f" % \
                (sum(self.acceptations[0:iteration]) / iteration))
        #print "iteration %i, acceptance rate until now: %.3f" % (iteration, sum(self.acceptations[0:iteration]) / iteration)
        if (random.uniform(size = 1) < self.proposalmixtureweights[0]):
            transformedthetastar = self.modeltheta.proposal(self.transformedtheta[:, newaxis], self.fixedproposalcovmatrix)[:, 0]
        else:
            transformedthetastar = self.modeltheta.proposal(self.transformedtheta[:, newaxis], self.adaptiveproposalcovmatrix)[:, 0]
        thetastar = self.modeltheta.untransform(transformedthetastar[:,newaxis])[:,0]
        proposedSIRs = ParallelSIRs(self.Nx, thetastar[:, newaxis], self.observations, self.modelx)
        proposedSIRs.first_step()
        proposedSIRs.next_steps()
        proposedloglike = proposedSIRs.getTotalLogLike()[0]
        proposedlogomega = proposedloglike + self.modeltheta.priorlogdensity(transformedthetastar)
        currentlogomega = self.totalLogLike + self.modeltheta.priorlogdensity(self.transformedtheta)
        acceptation  = (log(random.uniform(size = 1)) < (proposedlogomega - currentlogomega))
#        if not(acceptation):
#            print "\nproposed loglikelihood:", proposedloglike
#            print "\ncurrent loglikelihood:", self.totalLogLike
#            print "\naccept:", acceptation
        if acceptation:
            self.theta = thetastar.copy()
            self.transformedtheta = transformedthetastar.copy()
            self.totalLogLike = proposedloglike.copy()
        self.chain[:, iteration] = self.theta.copy()
        self.transformedchain[:, iteration] = self.transformedtheta.copy()
        if iteration > 10:
            self.adaptiveproposalcovmatrix = (5.6644 / self.thetadim) * cov(self.transformedchain[:, 1:iteration])
            #print "\n", self.adaptiveproposalcovmatrix
        self.loglikechain[iteration] = self.totalLogLike
        self.acceptations[iteration] = acceptation
    def allsteps(self):
        for t in xrange(1, self.nbiterations):
            self.PMCMCstep(t)
        self.globalacceptrate = sum(self.acceptations) / self.nbiterations
        print "...done!"
        print "Final acceptance rate: ",  self.globalacceptrate
        print "Fixed covariance matrix:"
        print self.fixedproposalcovmatrix
        print "Final adaptive covariance matrix:"
        print self.adaptiveproposalcovmatrix
    def getResults(self):
        resultsDict = {"trueparameters": self.modelx.parameters, \
                "nbparameters" : self.modeltheta.parameterdimension, \
                "N" : self.Nx, "T": self.observations.shape[0], \
                "chain": self.chain, \
                "observations": self.observations, \
                "adaptiveproposalcovmatrix": self.adaptiveproposalcovmatrix, \
                "globalacceptrate": self.globalacceptrate, "acceptations": self.acceptations}
        return resultsDict



