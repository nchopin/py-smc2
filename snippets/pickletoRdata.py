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
import numpy
from numpy import random, power, sqrt, exp, zeros, \
        ones, mean, average, prod, log, sum, repeat, \
        array, float32, int32, cov, isnan, load, savez
import cPickle

class RFeeder:
    def __init__(self):
        self.txt = "rm(list = ls())" + "\n"
    def r(self, newtext):
        self.txt += newtext + "\n"
#    def savetxt(self, txtfilename):
#        self.txtfilename = txtfilename
#        f = open(self.txtfilename, "w")
#        f.write(self.txt)
#        f.close()
#    def saveRData(self, RDatafilename):
#        txtfilename = "pickle2RData-tmp.R"
#        self.savetxt(txtfilename)
#        import subprocess
#        subprocess.call(["R", "CMD", "BATCH", "--vanilla", self.txtfilename, "pickle2RData-tmp.out"])
#        self.removeTxtFile()
#    def removeTxtFile(self):
#        os.remove(self.txtfilename)
#        os.remove("pickle2RData-tmp.out")


def R_repr(numpyarray):
    return "c(" + ",".join([str(element) for element in numpyarray]) + ")"

def convertarray(element, key, RF):
    s = element.shape
    if len(s) == 1:
        print "converting numpy vector '%s' to R FloatVector..." % key
        RF.r(""" %s <- %s """ % (key, R_repr(element)))
    elif len(s) == 2:
        print "converting numpy matrix '%s' to R matrix... " % key
        RF.r("""%s <- matrix(nrow = %i, ncol = %i)""" % (key, s[0], s[1]))
        for row in range(s[0]):
            RF.r(""" %s[%i+1,] <- %s """ % (key, row, R_repr(element[row,:])))
        #for column in range(s[1]):
        #    RF.r(""" %s[,%i+1] <- %s """ % (key, column, R_repr(element[:,column])))
    elif len(s) == 3:
        print "converting numpy array (dim = 3) '%s' to R array..." % key
        RF.r("""%s <- array(dim = c(%i, %i, %i))""" % (key, s[0], s[1], s[2]))
        for firstcomponent in range(s[0]):
            for secondcomponent in range(s[1]):
                RF.r(""" %s[%i+1, %i+1,] <- %s """ % \
                        (key, firstcomponent, secondcomponent, R_repr(element[firstcomponent, secondcomponent, :])))
    else:
        print "Don't know how to convert %s !!!" % key


def dictionary2RDataWithoutRPY(d, RDatafilename):
    RF = RFeeder()
    for key in d.keys():
        element = d[key]
        if type(element) is numpy.ndarray:
            convertarray(element, key, RF)
        elif isinstance(element, list):
            print "converting python list '%s'to R vector..." % key
            if len(element) == 0:
                RF.r(""" %s <- list()""" % key)
            else:
                if isinstance(element[0], numpy.ndarray):
                    print "it is a list of arrays"
                    RF.r(""" %s <- list()""" % key)
                    for indexsubelement, subelement in enumerate(element):
                        convertarray(subelement, key + "%i" % (indexsubelement + 1), RF)
                        RF.r(""" %s$"%i" <- %s"""  % (key, indexsubelement + 1, key + "%i" % (indexsubelement + 1)))
                        RF.r(""" rm(list = grep("all.+[0-9]", ls(), value = T)) """)
                else:
                    RF.r(""" %s <- %s """ % (key, R_repr(element)))
        elif isinstance(element, float):
            print "converting python float '%s' to R float..." % key
            RF.r(""" %s <- %.10f""" % (key, element))
        elif isinstance(element, int):
            print "converting python int '%s' to R int..." % key
            RF.r(""" %s <- %i""" % (key, element))
        elif isinstance(element, str):
            print "can't handle text '%s' ..." % key
            pass
        elif isinstance(element, dict):
            print "converting python dict '%s' to R list" % key
            for dkey, dvalue in element.items():
                convertarray(dvalue, key + dkey, RF)
        else:
            print "can't handle object '%' ..." % key

    RF.r("""save(list = ls(), file = "%s")""" % RDatafilename)
    print "saving..."
    txtbasename = "/tmp/unique"
    counter = 0
    txtfilename = "%s(%i).R" % (txtbasename, counter)
    while os.path.isfile(txtfilename):
        counter += 1
        txtfilename = "%s(%i).R" % (txtbasename, counter)
    f = open(txtfilename, "w")
    f.write(RF.txt)
    f.close()
    import subprocess
    subprocess.call(["R", "CMD", "BATCH", "--vanilla", txtfilename, "/dev/null"])
    os.remove(txtfilename)
    print "conversion finished"


def dictionary2RDataWithRPY(d, RDatafilename):
    RF = RFeeder()
    for key in d.keys():
        element = d[key]
        if type(element) is numpy.ndarray:
            convertarray(element, key, RF)
        elif isinstance(element, list):
            print "converting python list '%s'to R vector..." % key
            if len(element) == 0:
                RF.r(""" %s <- list()""" % key)
            else:
                if isinstance(element[0], numpy.ndarray):
                    print "it is a list of arrays"
                    RF.r(""" %s <- list()""" % key)
                    for indexsubelement, subelement in enumerate(element):
                        convertarray(subelement, key + "%i" % (indexsubelement + 1), RF)
                        RF.r(""" %s$"%i" <- %s"""  % (key, indexsubelement + 1, key + "%i" % (indexsubelement + 1)))
                        RF.r(""" rm(%s) """ % (key + "%i" % (indexsubelement + 1)))
                else:
                    RF.r(""" %s <- %s """ % (key, R_repr(element)))
        elif isinstance(element, float):
            print "converting python float '%s' to R float..." % key
            RF.r(""" %s <- %.10f""" % (key, element))
        elif isinstance(element, int):
            print "converting python int '%s' to R int..." % key
            RF.r(""" %s <- %i""" % (key, element))
        elif isinstance(element, str):
            print "can't handle text '%s' ..." % key
            pass
        elif isinstance(element, dict):
            print "converting python dict '%s' to R list" % key
            for dkey, dvalue in element.items():
                convertarray(dvalue, key + dkey, RF)
        else:
            print "can't handle object '%' ..." % key

    print "saving..."
    RF.r("""save(list = ls(), file = "%s")""" % RDatafilename)
    import rpy2.robjects as robjects
    robjects.r(RF.txt)
    print "conversion finished"

