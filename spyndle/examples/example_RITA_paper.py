# -*- coding: utf-8 -*-

"""
    Code showing a complete working example of how to use the spyndle toolbox.
    It has been developped for the presentation of paper [1].
    
    Copyright (C) 2013  Christian O'Reilly

    For personnal, educationnal, and research purpose, this software is 
    provided under the GNU GPL (V.3) license: you can redistribute it and/or
    modify it under the terms of the version 3 of the GNU General Public 
    License as published by the Free Software Foundation.
          
    To use this software in commercial application, please contact the author.
    
    If this code is used for research purpose, the reference [1] should be
    cited in the derived publication.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


 Author: Christian O'Reilly (christian.oreilly@umontreal.ca)
 Date  : July 1, 2013

 [1] O'Reilly, C. & Nielsen, T. The Spyndle toolbox: an open source software 
     package for the investigation of EEG transient events, Submitted to the 
     Revue internationale des technologies avancées, August 2013.

"""


from spyndle.io import EDFReader
from spyndle.detector import SpindleDectectorRMS
from spyndle.propagation import computeXCST, computeSPF_allFiles, \
    computeAveragePropagation, applyRejectionC3
from spyndle.EEG import getAdjacentsPairs_10_20, plotArrowsBidirect, \
    plotColorMap
from spyndle.EEG.system10_20 import getEEGChannels

from urllib import  urlretrieve,  ContentTooShortError
from scipy import unique, array
import pandas 
from datetime import datetime                                                  
A = datetime.now()                                                             

fileName        = "RITA_example.BDF"
detectionStages = ["Sleep stage 2"]
eventName       = "SpindleRMS"

print "Loading the data file from Internet. This may take some time, " \
      "the file is about 372 MB."
url = "https://bitbucket.org/christian_oreilly/spyndle/"\
      "downloads/RITA_example.BDF"
try:
    urlretrieve(url, fileName)
except ContentTooShortError:
     print "The retreived file is shorter than expected. The download as "\
           "probably been interrupted"
     exit
B = datetime.now()                                                             
print B-A                                                                      

print "Reading the .bdf file..."
readerEDF    = EDFReader(fileName)
listChannels = getEEGChannels(readerEDF.getChannelLabels())

print "Detecting spindles..."
detector = SpindleDectectorRMS(readerEDF)
detector.setDetectionStages(detectionStages)
detector.detectSpindles(listChannels)
detector.saveSpindle(readerEDF, eventName)
C = datetime.now()                                                             
print C-B                                                                      


print "Computing synchrone comparisons..."
computeXCST(EDFReader, fileName, "", eventName, listChannels, 
                 delta =0.5, beforePad=0.5, afterPad=0.5, 
                 offset=0.0, eventProperties = [],
                 resFilesPrefix="spindleDelays", similIndexType="euclidean")
D = datetime.now()                                                             
print D-C                                                                      

print "Computing asynchrone comparisons..."
computeXCST(EDFReader, fileName, "", eventName, listChannels, 
                 delta =0.5, beforePad=0.5, afterPad=0.5, 
                 offset=5.0, eventProperties = [],
                 resFilesPrefix="spindleDelays_offset", similIndexType="euclidean")

E = datetime.now()                                                             
print E-D                                                                      

print "Computing SPFs..."
computeSPF_allFiles("", offsetFilePattern="spindleDelays_offset_*.BDF_*.txt",
                    alpha = 10.0, verbose = False)
F = datetime.now()                                                             
print F-E                                                                      


print "Computing averages..."    
computeAveragePropagation("", aggeragationlevels = ["ref", "test"],
                  observationVariables = ["delay", "duration", "similarity"], 
                  pattern="correctedData_*.csv",  verbose=False)
G = datetime.now()                                                             
print G-F                                                                      

propData = pandas.read_csv("meanData.csv", sep=";")

print "Plotting the results in result_example_a.png..."
medData    = propData.groupby("ref").median()
electrodes = array(medData.index, dtype=str)
plotColorMap(electrodes, medData.duration_mean, "result_example_a.png")
                          
print "Plotting the results in result_example_b.png..."
#propData = applyRejectionC3(propData)
adjPairs = getAdjacentsPairs_10_20(unique(array(propData.ref)))
medData  = propData.groupby(["ref", "test"]).median()["delay_mean"]
print medData  
plotArrowsBidirect(medData, adjPairs, filename="result_example_b.png")
