**NEWS: If you're interested in this package, you should check this out: http://libbi.org/**

Python implementation of the SMC^2 algorithm described in Chopin, Jacob & Papaspiliopoulos (2011) (first version available here: http://arxiv.org/abs/1101.1528)

The package requires recent versions of python, numpy, scipy, a C compiler like GCC in order to use scipy.weave (to give a serious speed boost).

For the plots, R needs to be installed, along with the ggplot2 library. The program can still run without R, the results might be saved in RData format (to be handled with R) or pickle format (to be handled by the python cPickle package).

There's now a start of a manual, in the manual/ subfolder and here:
http://py-smc2.googlecode.com/hg/manual/manual.pdf

The zip source archive is the easiest way to retrieve the code. Otherwise the mercurial repository is generally more up to date.