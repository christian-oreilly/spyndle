# -*- coding: utf-8 -*-

import pylab
import os,sys
parentdir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,parentdir + "\\..\\..\\..") 

from scipy import zeros, where

from spyndle.io import HarmonieReader
from spyndle.io import EDFReader
from spyndle.detector import SpindleDectectorRMS
from spyndle.detector import DetectorEvaluator
from spyndle import cycleDefinitions, computeDreamCycles

print "Reading the .sig file..."
readerSIG =  HarmonieReader(parentdir + "/test.SIG")

#print "Conversion: .sig -> .bdf..."
#readerSIG.saveAsEDF(parentdir + "/test.BDF", "BDF")

print "Reading the .bdf file..."
readerEDF = EDFReader(parentdir + "/test.BDF")


listChannels = readerEDF.getChannelLabels()[0:2]

aeschbachCycleDef = cycleDefinitions()
aeschbachCycleDef.setAeschbach()

cyclesSIG = computeDreamCycles([e for e in readerSIG.events if e.groupeName == "Stage"], aeschbachCycleDef)
cyclesEDF = computeDreamCycles([e for e in readerEDF.events if e.groupeName == "Stage"], aeschbachCycleDef)


listChannels = [channel for channel in listChannels if channel in readerSIG.getChannelLabels()] 

readerSIG.pickleCompleteRecord(listChannels)   
   
for channel in listChannels:    
    
    signalSIG      = readerSIG.readChannel(channel, True).signal
    TSIG           = readerSIG.getChannelTime(channel)
    
    signalEDF      = readerEDF.readChannel(channel, False).signal
    TEDF           = readerEDF.getChannelTime(channel)
   
    endTime = min(TEDF[-1], TSIG[-1])
   
    print channel
    print len(TSIG), len(TEDF)
    print TSIG
    print TEDF
   
    offset = 900   
   
    indSIG = where( (TSIG >= endTime-30-offset)*(TSIG <= endTime-offset) )[0]
    indEDF = where( (TEDF >= endTime-30-offset)*(TEDF <= endTime-offset) )[0]

    print len(indSIG), len(indEDF)
       
    pylab.plot(TEDF[indEDF], signalEDF[indEDF])
    pylab.plot(TSIG[indSIG], signalSIG[indSIG])

    pylab.show()
       
   
   
   
"""   
    ################################## FILTERING ##########################
    if verbose:   print "Filtering..."   
    
    # Defining EEG filters
    bandPassFilter = Filter(fs)
    bandPassFilter.create(low_crit_freq=self.lowFreq, 
                          high_crit_freq=self.highFreq, order=1001, 
                          btype="bandpass", ftype="FIR", useFiltFilt=True)          

    # filtering can take a lot of memory. By making sure that the 
    # garbage collector as passed just before the filtering, we
    # increase our chances to avoid a MemoryError  
    gc.collect()                      
    signal = bandPassFilter.applyFilter(signal)     

    
    ################################# RMS COMPUTATION #####################
    if verbose:  print "RMS computation..."
    
    # Using teager operator?
    if self.useTeager :
         signal = concatenate(([0], signal[1:-1]**2 - signal[:-2]*signal[2:], [0]))
    
    
    windowNbSample = int(round(self.averagingWindowSize*fs))
    if mod(windowNbSample, 2) == 0 : # We need an odd number.
        windowNbSample += 1

    self.RMS = self.averaging(signal, windowNbSample)
    

    ############################### SPINDLE DETECTION #####################  
    if verbose:  print "Spindle detection..."

    # Thereshold are computed for each sleep cycles and each sleep stage
    # since the average amplitude of EEG signal can vary accoss night
    # and stages.

    channelTime = self.reader.getChannelTime(channel)    
    assert(len(channelTime) == len(signal))
    
    for cycle in cycles :
        stageEvent = self.reader.getEventsByTime(cycle.timeStart(), cycle.timeEnd())               
        stageEvent = filter(lambda e: e.groupeName == "Stage", stageEvent)     
        
        currentStageEvents = []
        for stage in self.detectionStages :
            currentStageEvents.extend(filter(lambda e: e.name.lower() == stage.lower(), stageEvent)  ) 

        if len(currentStageEvents) :
            samplesIndexes = concatenate([where((channelTime >= event.timeStart())*
                                                (channelTime <= event.timeEnd()))[0] 
                                          for event in currentStageEvents])


            treshold = mquantiles(self.RMS[samplesIndexes], self.quantileThreshold)[0]
            
            #0.375179523976 #
            #print "Threshold:", treshold
            
            spinMarkers = concatenate(([0], diff((self.RMS[samplesIndexes] > treshold).astype(int))))
            
            # The signal indexes corresponding to spinMarkers are kept in samplesIndexes                
            
            startSpinInd = samplesIndexes[where(spinMarkers == 1)[0]]
            stopSpinInd  = samplesIndexes[where(spinMarkers == -1)[0]]
            
            #The first marker should be a start marker.                
            if stopSpinInd[0] < startSpinInd[0]:
                stopSpinInd = stopSpinInd[1:]

            indEnd = min(len(startSpinInd), len(stopSpinInd))
            
            startSpinInd = startSpinInd[:indEnd]
            stopSpinInd  = stopSpinInd[:indEnd] 

            assert(all(stopSpinInd - startSpinInd > 0))
            assert(all(startSpinInd[1:] - stopSpinInd[:-1] > 0))

            gapToKeepInd = where(channelTime[startSpinInd[1:]] - channelTime[stopSpinInd[:-1]] > self.maxAllowableGapBellowThreshold)[0]
            
            startSpinInd = startSpinInd[concatenate(([0], gapToKeepInd+1))]
            stopSpinInd  = stopSpinInd[concatenate((gapToKeepInd, [len(stopSpinInd)-1]))]
            
            
            ## Careful : Non contiguous datas not presently supported. 
            duration = channelTime[stopSpinInd] - channelTime[startSpinInd]
                           
            startSpinInd = startSpinInd[duration >= self.minimalSpindleDuration]
            stopSpinInd  = stopSpinInd[duration >= self.minimalSpindleDuration]
            
            newSpindles = [DetectedSpindle(channel, start, end) for start, end in zip(channelTime[startSpinInd], channelTime[stopSpinInd])]   
            self.detectedSpindles.extend(newSpindles) 
        #print [s.startTime() for s in self.detectedSpindles]            
"""