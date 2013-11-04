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

from spyndle.io import EDFReader, DatabaseMng
from spyndle.detector import SpindleDectectorRMS
from spyndle.propagation import SPFEvaluator
from spyndle.EEG import getAdjacentsPairs_10_20, plotArrowsBidirect, \
    plotColorMap, getActiveElectrodes

from urllib import  urlretrieve,  ContentTooShortError
from scipy import unique
from datetime import datetime                                                  

fileName        = "RITA_example.BDF"
detectionStages = ["Sleep stage 2"]
eventName       = "SpindleRMS"

print "Loading the data file from Internet. This may take some time, " \
      "the file is about 372 MB."
url = "https://bitbucket.org/christian_oreilly/spyndle/"\
      "downloads/RITA_example.BDF"


t1 = datetime.now()                                                             
try:
    urlretrieve(url, fileName)
except ContentTooShortError:
     print "The retreived file is shorter than expected. The download as "\
           "probably been interrupted"
     exit
t2 = datetime.now()                                                             
print "Duration of data downloading : ", t2-t1                                                                      

print "Creating a  database for recording final and intermediate results..."
dbPath = "C:\\Python27\\Lib\\site-packages\\spyndle\\examples\\RITA.db"
dbMng = DatabaseMng("sqlite:///" + dbPath)

print "Reading the .bdf file..."
readerEDF    = EDFReader(fileName)

print "Detecting spindles..."
t1 = datetime.now()  
detector = SpindleDectectorRMS(readerEDF)
detector.setDetectionStages(detectionStages)
detector.detectSpindles()
detector.saveSpindle(readerEDF, eventName, dbSession=dbMng.session)
t2 = datetime.now()                                                             
print "Duration of spindle detection: ", t2-t1                                                                      

print "Plotting the results in result_example_a.png..."
filteringDict = {"psgNight":fileName, "eventName":eventName}
spindles = dbMng.dmm.getTransientEvents(filteringDict, pandasFormat=True)                       
spindles["electrode"] = getActiveElectrodes(spindles["channelName"])
medDat = spindles.groupby("electrode").median()
plotColorMap(unique(spindles["electrode"]), medDat.duration, "example_a.png")

evaluator = SPFEvaluator(fileName, eventName, dbSession = dbMng.session, verbose=True)
         
print "Computing synchrone comparisons..."
t1 = datetime.now()  
evaluator.computeXCST(EDFReader, offset=0.0)
t2 = datetime.now()                                                              
print "Duration of synchronous comparison computation: ", t2-t1                                                                     

print "Computing asynchrone comparisons..."
t1 = datetime.now()  
evaluator.computeXCST(EDFReader, offset=0.5)
t2 = datetime.now()                                                           
print "Duration of asynchronous comparison computation: ", t2-t1                                                                    

t1 = datetime.now()  
print "Computing SPFs..."
evaluator.computeSPF()            

print "Computing propagation averages..." 
evaluator.computeAveragePropagation(minNbValid=20)
t2 = datetime.now()                                                           
print "Duration of spindle propagation field computation: ", t2-t1                                                                     

print "Plotting the results in result_example_b.png..."
propagations = evaluator.getAveragePropagation(applyC3=True, applyC4=True)
print propagations
propagations["ref"] = getActiveElectrodes(propagations["sourceChannelName"])
propagations["test"] = getActiveElectrodes(propagations["sinkChannelName"])
adjPairs = getAdjacentsPairs_10_20(unique(propagations["ref"]))
medDat = propagations.groupby(["ref", "test"]).median()["delay_mean"]
plotArrowsBidirect(medDat, adjPairs, filename="example_b.png")

