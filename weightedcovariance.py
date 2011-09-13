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

def computeCovarianceAndMean(X, unnormalizedw):
    w = unnormalizedw / numpysum(unnormalizedw)
    weightedmean = average(X, weights = w, axis = 0)
    diagw = diag(w)
    part1 = dot(transpose(X), dot(diagw, X))
    Xtw = dot(transpose(X), w[:, newaxis])
    part2 = dot(Xtw, transpose(Xtw))
    numerator = part1 - part2
    denominator = 1 - numpysum(w**2)
    weightedcovariance = numerator / denominator
    # increase a little bit the diagonal to prevent degeneracy effects
    #weightedcovariance += diag(zeros(self.modeltheta.parameterdimension) + 10**(-4)/self.modeltheta.parameterdimension)
    return {"mean": weightedmean, "cov": weightedcovariance}


def computeCovarianceAndMean2(X, unnormalizedw):
    weights = unnormalizedw / numpysum(unnormalizedw)
    Xbar = average(X, weights = w, axis = 0)
    code = \
    """
    int row,col;
    for (row = 0; row < d(0); row++)
    {
        for(col = 0; col < d(0); col++)
        {
          for (int index = 0; index < N(0); index ++){
            covariance(row, col) += weights(index) * (X(index, row) - Xbar(row)) * (X(index, col) - Xbar(col));
          }
        }
    }
    """
    #y = array([y])
    #Nx = states.shape[0]
    #Ntheta = states.shape[2]
    #weights = zeros((Nx, Ntheta))
    #noise = random.normal(size = (Nx, Ntheta), loc = 0, scale = 1)
    d = X.shape[1]
    print "d:", d
    covariance = zeros((d, d))
    d = array([d])
    N = array([X.shape[0]])
    weave.inline(code,['covariance', 'd', 'N', 'Xbar', 'X', 'weights'], \
        type_converters=weave.converters.blitz, libraries = ["m"])
    weightedcovariance = covariance / (1 - numpysum(power(weights, 2)))
    return {"mean": Xbar, "cov": weightedcovariance}
#    diagw = diag(w)
#    part1 = dot(transpose(X), dot(diagw, X))
#    Xtw = dot(transpose(X), w[:, newaxis])
#    part2 = dot(Xtw, transpose(Xtw))
#    numerator = part1 - part2
#    denominator = 1 - numpysum(w**2)
#    weightedcovariance = numerator / denominator
    # increase a little bit the diagonal to prevent degeneracy effects
    #weightedcovariance += diag(zeros(self.modeltheta.parameterdimension) + 10**(-4)/self.modeltheta.parameterdimension)
    #return {"mean": weightedmean, "cov": weightedcovariance}

random.seed(17)
N = 10000
d = 3
covariance = diag(repeat(1, d))
mean = repeat(0, d)
X = random.multivariate_normal(mean = mean, cov = covariance, size = N)
w = random.exponential(size = N)
#print X
#print w
import cProfile
cProfile.run("""
results = computeCovarianceAndMean(X, w)
results2 = computeCovarianceAndMean2(X, w)
""", "prof")
#print results["cov"]
#print results["mean"]
import pstats
p = pstats.Stats('prof')
#p.strip_dirs().sort_stats(-1).print_stats()
#p.sort_stats('name')
#p.print_stats()
p.sort_stats('cumulative').print_stats(10)
p.sort_stats('time').print_stats(10)

#print computeCovarianceAndMean(X, w)["cov"]
#print computeCovarianceAndMean2(X, w)["cov"]

