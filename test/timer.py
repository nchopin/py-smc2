import time

tic = time.time()
import numpy
a = numpy.random.normal(size = 100000)
print time.time() - tic 
