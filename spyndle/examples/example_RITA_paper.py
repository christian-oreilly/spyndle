from spyndle.io import EDFReader
from spyndle.detector import SpindleDectectorRMS
from spyndle.propagation import computeXCST, computeSPF_allFiles,computeAveragePropagation
from spyndle.EEG import getAdjacentsPairs_10_20, plotArrowsBidirect, plotColorMap
from spyndle.EEG.system10_20 import getEEGChannels

from urllib import  urlretrieve,  ContentTooShortError
from scipy import unique, array
import pandas 
from datetime import datetime                                                  ##########################################
A = datetime.now()                                                             ##########################################

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
B = datetime.now()                                                             ##########################################
print B-A                                                                      ##########################################

print "Reading the .bdf file..."
readerEDF    = EDFReader(fileName)
listChannels = getEEGChannels(readerEDF.getChannelLabels())

print "Detecting spindles..."
detector = SpindleDectectorRMS(readerEDF)
detector.setDetectionStages(detectionStages)
detector.detectSpindles(listChannels)
detector.saveSpindle(readerEDF, eventName)
C = datetime.now()                                                             ##########################################
print C-B                                                                      ##########################################


print "Computing synchrone comparisons..."
computeXCST(EDFReader, fileName, "", eventName, listChannels, 
                 delta =0.5, beforePad=0.5, afterPad=0.5, 
                 offset=0.0, eventProperties = [],
                 resFilesPrefix="spindleDelays")
D = datetime.now()                                                             ##########################################
print D-C                                                                      ##########################################

print "Computing asynchrone comparisons..."
computeXCST(EDFReader, fileName, "", eventName, listChannels, 
                 delta =0.5, beforePad=0.5, afterPad=0.5, 
                 offset=5.0, eventProperties = [],
                 resFilesPrefix="spindleDelays_offset")

E = datetime.now()                                                             ##########################################
print E-D                                                                      ##########################################

print "Computing SPFs..."
computeSPF_allFiles("", offsetFilePattern="spindleDelays_offset_*.BDF_*.txt",
                    alpha = 10.0, verbose = False)
F = datetime.now()                                                             ##########################################
print F-E                                                                      ##########################################


print "Computing averages..."    
computeAveragePropagation("", aggeragationlevels = ["ref", "test"],
                          observationVariables = ["delay", "duration", "similarity"], 
                          pattern="correctedData_*.csv",  verbose=True)
G = datetime.now()                                                             ##########################################
print G-F                                                                      ##########################################

propData = pandas.read_csv("meanData.csv", sep=";")

print "Plotting the results in result_example_a.png..."
medData = propData.groupby("ref").median()
electrodes = array(medData.index, dtype=str)
plotColorMap(electrodes, medData.duration_mean, "result_example_a.png")
                            
print "Plotting the results in result_example_b.png..."
adjPairs = getAdjacentsPairs_10_20(unique(array(propData.ref)))
medData = propData.groupby(["ref", "test"]).median()["delay_mean"]
plotArrowsBidirect(medData, adjPairs, filename="result_example_b.png")
