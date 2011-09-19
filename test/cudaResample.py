from __future__ import division
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32, zeros_like, newaxis, double
from numpy import max as numpymax
from numpy import mean as numpymean
from scipy.stats import norm, truncnorm, gamma
import scipy.weave as weave
import os
import math
import pycuda.autoinit
import pycuda.driver as drv
from pycuda.compiler import SourceModule


def resample2D(array, weights, Nx, xdim, Ntheta):
    """ Indexed of resampled particles
    (deterministic scheme) """
    code = \
    """
    for(int i = 0; i < Ntheta; i++)
    {
        int j = 0;
        float csw = weights(0, i);
        for(int k = 0; k < Nx; k++)
        {
            while(csw < u(i))
            {
                j++;
                csw += weights(j, i);
            }
            for (int l = 0; l < xdim; l++){
                newarray(k, l, i) = array(j, l, i);
            }
            u(i) = u(i) + 1.;
        }
    }
    """
    u = random.uniform(size = Ntheta, low = 0, high = 1).astype(float32)
    print u[0:10],"..."
    weights = double(weights * Nx / sum(weights, axis = 0))
    weights = weights.astype(float32)
    array = array.astype(float32)
    newarray = zeros_like(array).astype(float32)
    weave.inline(code,['u','Nx', 'Ntheta','weights', 'array', 'newarray', 'xdim'], type_converters=weave.converters.blitz)
    return newarray

Nx = 500
Ntheta = 1000
xdim = 2
states = random.normal(size = (Nx, xdim, Ntheta))
uw = abs(random.normal(size = (Nx, Ntheta)))
#print states

resampleKernelTemplate = """
        __global__ void resample(float* u, float* newstates, float* states, float* weights, int* Nx, int* Ntheta, int* xdim)
        {
            int itheta = threadIdx.x + blockDim.x * blockIdx.x;
            if (itheta < *Ntheta){ 
                float csw = weights[itheta];
                int j = 0;
                for(int k = 0; k < (*Nx); k++)
                {
                    while(csw < u[itheta])
                    {
                        j++;
                        csw += weights[itheta + (*Ntheta) * j];
                    }
                    for (int l = 0; l < *xdim; l++){
                        newstates[l* (*Ntheta) + k * ((*Ntheta) * (*xdim)) + itheta] = states[l* (*Ntheta) + j * ((*Ntheta) * (*xdim)) + itheta];
                    }
                    u[itheta] = u[itheta] + 1.;
                }
            }
        }
"""
resampeKernel = resampleKernelTemplate % {"Ntheta": Ntheta}
resampleGF= SourceModule(resampeKernel).get_function("resample")
def resampleCUDA(states, weights, Nx, xdim, Ntheta):
    THREADS_PER_BLOCK_X = 512
    states = states.astype(float32)
    newstates = zeros_like(states).astype(float32)
    u = random.uniform(size = Ntheta, low = 0, high = 1).astype(float32)
    print u[0:10],"..."
    weights = double(weights * Nx / sum(weights, axis = 0))
    num_blocks_x = int(math.ceil((Ntheta + THREADS_PER_BLOCK_X - 1)/ THREADS_PER_BLOCK_X))
    resampleGF(drv.In(u), \
            drv.InOut(newstates), \
            drv.In(states), \
            drv.In(weights.astype(float32)), \
            drv.In(array(Nx, dtype = int32)), \
            drv.In(array(Ntheta, dtype = int32)), \
            drv.In(array(xdim, dtype = int32)), \
            block = (THREADS_PER_BLOCK_X, 1, 1), grid = (num_blocks_x, 1))
    return newstates


import cProfile
cProfile.run("""
random.seed(923)
resCUDA = resampleCUDA(states, uw, Nx, xdim, Ntheta)
random.seed(923)
res = resample2D(states, uw, Nx, xdim, Ntheta)
""", "prof")
import pstats
p = pstats.Stats('prof')
p.sort_stats('cumulative').print_stats(10)
p.sort_stats('time').print_stats(10)
print numpymax(resCUDA - res)
print numpymean(resCUDA - res)
print sum((resCUDA - res) > 0.1)


