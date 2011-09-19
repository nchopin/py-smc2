from __future__ import division
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32, zeros_like, newaxis
from numpy import max as numpymax
from scipy.stats import norm, truncnorm, gamma
import scipy.weave as weave
import os
import math
import pycuda.autoinit
import pycuda.driver as drv
import numpy

from pycuda.compiler import SourceModule
mod = SourceModule("""
__global__ void multiply_them(float *dest, float *a, float *b)
{
  const int i = threadIdx.x;
  dest[i] = a[i] * b[i];
}
""")

multiply_them = mod.get_function("multiply_them")

random.seed(820)

a = numpy.random.randn(400).astype(numpy.float32)
b = numpy.random.randn(400).astype(numpy.float32)

dest = numpy.zeros_like(a)
multiply_them(
        drv.Out(dest), drv.In(a), drv.In(b),
        block=(400,1,1), grid=(1,1))

print dest-a*b

#! /usr/bin/env python
# -*- coding: utf-8 -*-



def transitionAndWeight(states, y, parameters, t):
    code = \
    """
    float tempmeasure1;
    float tempmeasure2;
    float temptransition;
    for (int j = 0; j < Ntheta; j++)
    {
        tempmeasure1 = -0.9189385 - 0.5 * log(parameters(1, j));
        tempmeasure2 = -0.5 / parameters(1, j);
        temptransition = sqrt(parameters(0, j));
        for(int k = 0; k < Nx; k++)
        {
            states(k, 0, j) = states(k, 0, j) + temptransition * noise(k, j);
            weights(k, j) = tempmeasure1 + 
            tempmeasure2 * ((double) y(0) - states(k, 0, j)) * ((double) y(0) - states(k, 0, j));
        }
    }
    """
    y = array([y])
    Nx = states.shape[0]
    Ntheta = states.shape[2]
    weights = zeros((Nx, Ntheta))
    noise = random.normal(size = (Nx, Ntheta), loc = 0, scale = 1)
    weave.inline(code,['Nx', 'Ntheta', 'states', 'y', 'parameters', 'noise', 'weights'], type_converters=weave.converters.blitz, libraries = ["m"])
    return {"states": states , "weights": weights}



Nx = 400
Ntheta = 300
statedim = 2
states = random.normal(size = (Nx, statedim, Ntheta))
#print "original states"
#print states
y = 1
parameters = ones((2, Ntheta))
#print "CPU"
#random.seed(820)
#print transitionAndWeight(states.copy(), y, parameters, 0)
##print "CPU"
#random.seed(820)
#CPUres = transitionAndWeight(states.copy(), y, parameters, 0)
#print CPUres
random.seed(820)
comboKernelTemplate = """
        __global__ void combo(float* y, float* states, float *parameters, float* noise, int* Nx, int* Ntheta, int* statedim)
        {
            int i = threadIdx.x + blockDim.x * blockIdx.x;
            int j = threadIdx.y + blockDim.y * blockIdx.y;
            if (i < *Nx && j < *Ntheta){ 
                int idim = 0;
                int index = idim * (*Ntheta) + i * (*Ntheta * *statedim) + j;
                states[index] = states[index] + sqrt(parameters[0 * (*Ntheta) + j]) * noise[i * (*Ntheta) + j];
                noise[i * (*Ntheta) + j] = -0.9189385 - 0.5 * log(parameters[1 * (*Ntheta) + j]) - 
                0.5 / parameters[1 * (*Ntheta) + j] * (y[0] - states[index]) * (y[0] - states[index]);
            }
        }
"""
#import pycuda.driver as cuda
#import pycuda.autoinit
#from pycuda.compiler import SourceModule
comboKernel = comboKernelTemplate % {"Ntheta": Ntheta}
comboGF = SourceModule(comboKernel).get_function("combo")
def combo_CUDA(states, observation, parameters, comboGF):
    THREADS_PER_BLOCK_X = 16 
    THREADS_PER_BLOCK_Y = 16 
    y = array(observation, dtype = float32)
    states = states.astype(float32)
    parameters = parameters.astype(float32)
    Nx, statedim, Ntheta = states.shape
    #print Nx, Ntheta
    noise = random.normal(size = (Nx, Ntheta), loc = 0, scale = 1).astype(float32)
    #print "CUDA noise"
    #print noise
    num_blocks_x = int(math.ceil(Nx / THREADS_PER_BLOCK_X))
    num_blocks_y = int(math.ceil(Ntheta / THREADS_PER_BLOCK_Y))
    comboGF(drv.In(y), drv.InOut(states), \
            drv.In(parameters), \
            drv.InOut(noise), \
            drv.In(array(Nx, dtype = int32)), \
            drv.In(array(Ntheta, dtype = int32)), \
            drv.In(array(statedim, dtype = int32)), \
            block = (THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y, 1), grid = (num_blocks_x, num_blocks_y))
    return {"states": states, "weights": noise}

print "CUDA:"
import cProfile
cProfile.run("""
random.seed(820)
CPUres = transitionAndWeight(states.copy(), y, parameters, 0)
random.seed(820)
cudares = combo_CUDA(states, y, parameters, comboGF)
""", "prof")
import pstats
p = pstats.Stats('prof')
p.sort_stats('cumulative').print_stats(10)
p.sort_stats('time').print_stats(10)
print "diff"
print numpymax(CPUres["states"] - cudares["states"])
print numpymax(CPUres["weights"] - cudares["weights"])

