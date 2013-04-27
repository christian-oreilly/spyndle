# -*- coding: utf-8 -*-

import os,sys
parentdir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,parentdir + "\\..\\..\\..") 

from spyndle.io import HarmonieReader
from spyndle.io import EDFReader
from spyndle.detector import SpindleDectectorRMS
from spyndle.detector import DetectorEvaluator

print "Reading the .sig file..."
readerSIG =  HarmonieReader(parentdir + "/test.SIG")

#print "Conversion: .sig -> .bdf..."
#readerSIG.saveAsEDFA(parentdir + "/test.BDF", "BDF")

print "Reading the .bdf file..."
readerEDF = EDFReader(parentdir + "/test.BDF")


listChannels = readerEDF.getChannelLabels()[0:2]

detectorEDF = SpindleDectectorRMS(readerEDF, usePickled=False)
detectorSIG = SpindleDectectorRMS(readerSIG, usePickled=True)


evaluator = DetectorEvaluator(detectorSIG, detectorEDF)
evaluator.printEvaluation(listChannels, ["Sleep stage N2", "Sleep stage 2"])

