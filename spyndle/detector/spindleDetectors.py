# -*- coding: utf-8 -*-
"""
Created on Tue Jan 08 16:54:05 2013

@author: REVESTECH
"""

"""
    Implementation of spindle detectors. 

    Copyright (C) 2013  Christian O'Reilly

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.


    For personnal, educationnal, and research purpose, this software is 
    provided under the GNU GPL (V.3) license: you can redistribute it and/or
    modify it under the terms of the version 3 of the GNU General Public 
    License as published by the Free Software Foundation.
          
    To use this software in commercial application, please contact the author. 
    If used for research purpose, the reference [1] should be cited in the
    derived publication to refere the reader to the description 
    of the methodology.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


 Author: Christian O'Reilly (christian.oreilly@umontreal.ca)
 Date  : January 8, 2013

 [1] O'Reilly, C. & Nielsen, T. Sleep spindle detection: Automation and 
     performance evaluation using fine temporal resolution, Submitted to
     IEEE Transactions on Biomedical Engineering, april 2013. 

"""


from scipy import concatenate, zeros, mean, sqrt, mod, diff, where, fft
from scipy import array, arange, ones, unique
from scipy.stats.mstats import mquantiles
from scipy.fftpack import fftfreq
from scipy.io import savemat, loadmat
from scipy.integrate import trapz

from spyndle import Filter
from spyndle import cycleDefinitions, computeDreamCycles
from spyndle import computeMST

import os, gc, copy



# Class representing a detected spindle.
class DetectedSpindle:
    def __init__(self, channel, startSample, endSample):

        # Channel on which spindle has been detected
        self.channel = channel
        
        # Spindle starting sample
        self.startSample = startSample
        
        # Spindle ending sample
        self.endSample   = endSample

        self.RMSamp         = 0.0
        self.meanFreq       = 0.0
        self.timeDuration   = 0.0 
        self.sleepStage     = ""


    def computeRMS(self, reader, fmin=10, fmax=16):
        
        # The filters need the signal to be at least 3*nbTaps
        nbTaps = 1001        
        sampleDuration = max(nbTaps*3, self.endSample-self.startSample)        
        data = reader.read(self.channel, self.startSample, sampleDuration=sampleDuration)
        
        
        #print self.channel #, data
        signal      = data[1][0]
        fs          = data[0][0]                       # sampling rate                     
            
            
        # Defining EEG filters
        bandPassFilter = Filter(fs)
        bandPassFilter.create(low_crit_freq=fmin, high_crit_freq=fmax, order=nbTaps, btype="bandpass", ftype="FIR", useFiltFilt=True)          
             
        signal = bandPassFilter.applyFilter(signal)[0:(self.endSample-self.startSample)]                                
    
        self.RMSamp = sqrt(mean(signal**2))
        
    
                
        
    def computeMeanFreq(self, reader, fmin=10, fmax=16):
        data = reader.read(self.channel, self.startSample, sampleDuration=self.endSample-self.startSample)

        samplingRate = data[0][0]
        signal       = data[1][0]

        if signal.size < 512:
            signal = concatenate((signal, zeros(512 - signal.size)))
        
        
        FFT = abs(fft(signal))
        freqs = fftfreq(signal.size, 1.0/samplingRate)        
        
        FFT   = FFT[(freqs >= fmin)*(freqs <= fmax)]
        freqs = freqs[(freqs >= fmin)*(freqs <= fmax)]
        
        self.meanFreq = sum(freqs*FFT)/sum(FFT)
        
        
        
    def computeTimeDuration(self, reader):
        self.timeDuration = float(self.endSample - self.startSample)/reader.getSamplingRate(self.channel)


    def setStage(self, reader):
        # Select the stage where the spindle begin as the sleep stage
        # of the spindle.
        stageEvent23 = filter(lambda e: e.groupeName == "Stage" and
                                            e.sampleStart() <=  self.startSample and 
                                            e.sampleEnd() >= self.startSample, reader.events)                                 
        if len(stageEvent23) > 0 :
            self.sleepStage = stageEvent23[0].name.lower()
                
        else:
            print "No stage", self.startSample
            


# Generic spindle detector.
class SpindleDectector:
    
    
    def __init__(self):
        ###############################################################################
        # Detection patameters
        ###############################################################################
        
        # Sleep stages in which we want to detect spindles. Should be a list of
        # event names used to score the stages we want to extract spindle from.        
        self.detectionStages    = ["Stage2", "Stage3"]
        
        # Low-end of the spindle frequency band
        self.lowFreq  = 11.0
        
        # High-end of the spindle frequency band        
        self.highFreq    = 16.0
        
        # Minimal duration of a valid spindle
        self.minimalSpindleDuration = 0.5 # in seconds
        
        # Sleep cycle definition used to compute the quantile distribution
        # of amplitude per stage/cycles
        self.aeschbachCycleDef = cycleDefinitions()
        self.aeschbachCycleDef.setAeschbach()

        # List of detected spindles
        self.detectedSpindles = []

        ###############################################################################

    def computeRMS(self, reader, fmin=11, fmax=16):
        for spindle in self.detectedSpindles:        
            spindle.computeRMS(reader, fmin, fmax)
        
    # computeRMS2 is less elegant than computeRMS but is 40 times faster at 
    # esecution
    def computeRMS2(self, reader, fmin=11, fmax=16):
        listChannels = unique([s.channel for s in self.detectedSpindles])        
        
        # Pickle data for each channel separatelly so simplify and accelerate
        # the reading of large files.       
        reader.pickleCompleteRecord(listChannels)   
   
        for channel in listChannels:    
            data        = reader.readPickledChannel(channel)

            signal      = data.signal
            fs          = data.samplingRate                      # sampling rate    
            
            # Defining EEG filters   
            bandPassFilter = Filter(fs)
            bandPassFilter.create(low_crit_freq=fmin, high_crit_freq=fmax, order=1001, btype="bandpass", ftype="FIR", useFiltFilt=True)          
       
                   
            # filtering can take a lot of memory. By making sure that the 
            # garbage collector as passed just before the filtering, we
            # increase our chances to avoid a MemoryError  
            gc.collect()                     
            signal = bandPassFilter.applyFilter(signal)             

            channelSpindles = filter(lambda s: s.channel == channel, self.detectedSpindles)     
            for spindle in channelSpindles:    
                spindle.RMSamp = sqrt(mean(signal[spindle.startSample:(spindle.endSample+1)]**2))                
        

        
        
        
    def computeMeanFreq(self, reader, fmin=10, fmax=16):
        for spindle in self.detectedSpindles:        
            spindle.computeMeanFreq(reader, fmin, fmax)
        
    def computeTimeDuration(self, reader):
        for spindle in self.detectedSpindles:        
            spindle.computeTimeDuration(reader)




    def setStagesToSpindles(self, reader):
        for spindle in self.detectedSpindles:        
            spindle.setStage(reader)


    # setStagesToSpindles2 is less elegant thant setStagesToSpindles but is 50 times faster at 
    # esecution
    def setStagesToSpindles2(self, reader):
        stageIndicator, lstStages = self.computeStageIndicator(reader)
        for spindle in self.detectedSpindles:        
            spindle.sleepStage = lstStages[stageIndicator[spindle.startSample]]



    def computeStageIndicator(self, reader):

        stages = filter(lambda e: e.groupeName == "Stage" , reader.events)   
        lstStage = concatenate((unique([stage.name.lower() for stage in stages]), ["No stage"]))    
        stageIndicator = ones(reader.nbSamples)*(len(lstStage)-1)
            
        for i in range(len(lstStage)):
            stageEvent23 = filter(lambda e: e.name.lower() == lstStage[i], stages)    
            if len(stageEvent23) :                
                index = concatenate([range(event.sampleStart(), event.sampleEnd()) for event in stageEvent23])  
                stageIndicator[index] = i
                
        return stageIndicator.astype(int), lstStage
        


    def averaging(self, signal, windowNbSample):
     
        result = copy.copy(signal[:(len(signal)-windowNbSample+1)])
        
        for i in range(1, windowNbSample) : 
            result += signal[i:(i+len(signal)-windowNbSample+1)]

        return concatenate((ones(windowNbSample/2)*result[0], result, ones(windowNbSample/2)* result[-1]))/windowNbSample 



  
    # Detect every spindle in the channels listChannels of the file opened by the reader.
    def detectSpindles(self, reader, listChannels, verbose=True) :
        return

    # TODO: To implement. Used to save detected spindle in EEG data file. 
    def saveSpindleSTS(self, reader, eventName):
        
        for spindle in self.detectedSpindles:            
            reader.addEvent(eventName, spindle.startSample, spindle.endSample-spindle.startSample, spindle.channel)

    # TODO: To implement. Used to save detected spindle in EEG data file. 
    def saveSpindleTxt(self, fileName):

        try:
            f = open(fileName, "w")    
        except IOError:     
            print "Error: The selected file could not be open."
            exit()
        
        for spindle in self.detectedSpindles:            
            f.write(spindle.channel + ";" + str(spindle.startSample) + ";" +  
                    str(spindle.endSample) + ";" + str(spindle.RMSamp) + ";" + str(spindle.meanFreq )
                      + ";" + str(spindle.timeDuration) + ";" + str(spindle.sleepStage) + "\n")             




    # TODO: To implement. Used to save detected spindle in EEG data file. 
    def loadSpindleTxt(self, fileName):
        import csv

        try:
            cr = csv.reader(open(fileName, "rb"), delimiter=';')   
        except IOError:     
            print "Error: The selected file could not be open."
            exit()
        
        self.detectedSpindles = []
        for row in cr:
            spindle = DetectedSpindle(row[0], int(row[1]), int(row[2]))   
            self.detectedSpindles.append(spindle)
  
    # Set the sleep stages in which we want to detect spindles. stages should
    # be a list of valid stages.
    def setDetectionStages(self, stages):
        self.detectionStages = stages
        

 
 
 
 





class SpindleDectectorRMS(SpindleDectector):
    
    
    def __init__(self):
        SpindleDectector.__init__(self)
        
        ###############################################################################
        # Detection patameters
        ###############################################################################
    
        # Quantile of spindle amplitude used to determine spindle activity
        self.quantileThreshold = 0.95

        # We allow for a spindle activity region to contain samples having an
        # amplitude going bellow the quantileThreshold for a 
        # maximum of maxAllowableGapBellowThreshold seconds.
        self.maxAllowableGapBellowThreshold = 0.10 #0.07915697 #0.1 # in seconds
        
        # Width of the window used to computer the RMS amplitude.
        ### Choosing such that it always encompass at least two cycles of the 
        ### smaller frequency. We have a tradeoff between stability and accuracy...
        self.averagingWindowSize = 0.20  #0.508361 #2.0/lowFreq # In seconds

        # If true, teager operator of amplitude is used instead of amplitude
        # itself.
        self.useTeager = False

        ###############################################################################


    def setThreshold(self, value): self.quantileThreshold  = value

  
    # Detect every spindle in the channels listChannels of the file opened by the reader.
    def detectSpindles(self, reader, listChannels, verbose=True, usePickled=False) :

        if isinstance(listChannels, str) :
            listChannels = [listChannels]

        if not type(listChannels) is list and not type(listChannels) is tuple:
            raise TypeError
            
        # List of detected spindles
        self.detectedSpindles = []

        # Computing sleep cycles
        cycles = computeDreamCycles([e for e in reader.events if e.groupeName == "Stage"], self.aeschbachCycleDef)
    
    
        #################################### READING ##########################
        if verbose:   print "Start reading datafile..."
   

        listChannels = [channel for channel in listChannels if channel in reader.getChannelLabels()] 
        print listChannels

        # Pickle data for each channel separatelly so simplify and accelerate
        # the reading of large files.    
        if usePickled :
            reader.pickleCompleteRecord(listChannels)   
   
        for channel in listChannels:    
            if verbose:   print "Channel " + channel + "..."           
            
            data        = reader.readChannel(channel, usePickled=usePickled)

            signal      = data.signal
            fs          = data.samplingRate                      # sampling rate    
   
   
            ################################## FILTERING ##########################
            if verbose:   print "Filtering..."   
            
            # Defining EEG filters
                     
            
            bandPassFilter = Filter(fs)
            bandPassFilter.create(low_crit_freq=self.lowFreq, high_crit_freq=self.highFreq, order=1001, btype="bandpass", ftype="FIR", useFiltFilt=True)          

            # filtering can take a lot of memory. By making sure that the 
            # garbage collector as passed just before the filtering, we
            # increase our chances to avoid a MemoryError  
            gc.collect()                     
            
            #import datetime
            #a = datetime.datetime.now()        
            signal = bandPassFilter.applyFilter(signal, True)     
            #b = datetime.datetime.now()
            #print(b-a) # 0:02:20.937000
            #signal = bandPassFilter.applyFilter(signal, False)     
            #c = datetime.datetime.now()     
            #print(c-b) # 0:00:03.348000                
            
            
            
            ################################# RMS COMPUTATION #####################
            if verbose:  print "RMS computation..."
            
            # Using teager operator?
            if self.useTeager :
                 signal = concatenate(([0], signal[1:-1]**2 - signal[:-2]*signal[2:], [0]))
            
            
            windowNbSample = int(round(self.averagingWindowSize*fs))
            if mod(windowNbSample, 2) == 0 : # We want an odd number.
                windowNbSample += 1
    
            self.RMS = self.averaging(signal, windowNbSample)
            
            # Thereshold are computed for each sleep cycles and each sleep stage
            # since the average amplitude of EEG signal can vary accoss night
            # and stages.
            
            
            ############################### SPINDLE DETECTION #####################  
            if verbose:  print "Spindle detection..."
            
            for cycle in cycles :
                stageEvent = reader.getEventsBySample(cycle.sampleStart(), cycle.sampleEnd()) 
                stageEvent = filter(lambda e: e.groupeName == "Stage", stageEvent)   
                
                currentStageEvents = []
                for stage in self.detectionStages :
                    currentStageEvents.extend(filter(lambda e: e.name.lower() == stage.lower(), stageEvent)  ) 
    
                if len(currentStageEvents) :
                    samplesIndexes = concatenate([range(event.sampleStart(), event.sampleEnd()) for event in currentStageEvents])
    
                    treshold = mquantiles(self.RMS[samplesIndexes], self.quantileThreshold)[0]
                    
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

                    #assert(all(stopSpinInd - startSpinInd > 0))
                    #assert(all(startSpinInd[1:] - stopSpinInd[:-1] > 0))
    
                    gapToKeepInd = where(startSpinInd[1:] - stopSpinInd[:-1] > self.maxAllowableGapBellowThreshold*fs)[0]
                    
                    startSpinInd = startSpinInd[concatenate(([0], gapToKeepInd+1))]
                    stopSpinInd  = stopSpinInd[concatenate((gapToKeepInd, [len(stopSpinInd)-1]))]
                    
                    
                    ## Careful : data are supporte contiguous. 
                    nbSampleDuration = stopSpinInd - startSpinInd
                                   
                    startSpinInd = startSpinInd[nbSampleDuration/fs >= self.minimalSpindleDuration]
                    stopSpinInd  = stopSpinInd[nbSampleDuration/fs >= self.minimalSpindleDuration]
                    
                    newSpindles = [DetectedSpindle(channel, start, end) for start, end in zip(startSpinInd, stopSpinInd)]   
                    self.detectedSpindles.extend(newSpindles)  


 
 
 
 

 
 
class SpindleDectectorTeager(SpindleDectector):
    
    
    def __init__(self):
        SpindleDectector.__init__(self)
        self.detectStages = []
    

        ###############################################################################
        # Detection patameters
        ###############################################################################
    

        self.threshold = 0.6

        # We allow for a spindle activity region to contain samples having an
        # amplitude going bellow the quantileThreshold for a 
        # maximum of maxAllowableGapBellowThreshold seconds.
        self.maxAllowableGapBellowThreshold = 0.0 # in seconds

        ###############################################################################


    def setThreshold(self, value): self.threshold  = value

  
    # Detect every spindle in the channels listChannels of the file opened by the reader.
    def detectSpindles(self, reader, listChannels, verbose=True) :
        # List of detected spindles
        self.detectedSpindles = []
     
        #################################### READING ##########################
        if verbose:   print "Start reading datafile..."


        # Process only available channels...
        listChannels = [channel for channel in listChannels if channel in reader.getAvailableChannels()] 

        # Pickle data for each channel separatelly so simplify and accelerate
        # the reading of large files.       
        reader.pickleCompleteRecord(listChannels)   
   
        for channel in listChannels:    
            if verbose:   print "Channel " + channel + "..."           
            
            data        = reader.readPickledChannel(channel)

            signal      = data.signal
            fs          = data.samplingRate                      # sampling rate                        
                                                            

            ################################## FILTERING ##########################
            if verbose:   print "Filtering..."   
            
            
            # Defining EEG filters
            bandPassFilter = Filter(fs)
            bandPassFilter.create(low_crit_freq=self.lowFreq, high_crit_freq=self.highFreq, order=1001, btype="bandpass", ftype="FIR", useFiltFilt=True)          

            # filtering can take a lot of memory. By making sure that the 
            # garbage collector as passed just before the filtering, we
            # increase our chances to avoid a MemoryError  
            gc.collect()                     
            signal = bandPassFilter.applyFilter(signal)                      


            ########################## Computing Teager operator ######################
            if verbose:  print "TEAGER computation..."
            teager = concatenate(([0], signal[1:-1]**2 - signal[:-2]*signal[2:], [0]))
            
            
            ############################### SPINDLE DETECTION #####################  
            if verbose:  print "Spindle detection..."
            
            

            treshold = mean(teager)*self.threshold

            spinMarkers = concatenate(([0], diff((teager > treshold).astype(int))))         
            
            # The signal indexes corresponding to spinMarkers are kept in samplesIndexes                
            
            startSpinInd = where(spinMarkers == 1)[0]
            stopSpinInd  = where(spinMarkers == -1)[0]
            
            if len(stopSpinInd) and len(startSpinInd):
                #The first marker should be a start marker.                
                if stopSpinInd[0] < startSpinInd[0]:
                    stopSpinInd = stopSpinInd[1:]

                indEnd = min(len(startSpinInd), len(stopSpinInd))
                
                startSpinInd = startSpinInd[:indEnd]
                stopSpinInd  = stopSpinInd[:indEnd] 

                #assert(all(stopSpinInd - startSpinInd > 0))
                #assert(all(startSpinInd[1:] - stopSpinInd[:-1] > 0))

                gapToKeepInd = where(startSpinInd[1:] - stopSpinInd[:-1] > self.maxAllowableGapBellowThreshold*fs)[0]
                
                startSpinInd = startSpinInd[concatenate(([0], gapToKeepInd+1))]
                stopSpinInd  = stopSpinInd[concatenate((gapToKeepInd, [len(stopSpinInd)-1]))]
                
                
                ## Careful : data are supporte contiguous. 
                nbSampleDuration = stopSpinInd - startSpinInd
                               
                startSpinInd = startSpinInd[nbSampleDuration/fs >= self.minimalSpindleDuration]
                stopSpinInd  = stopSpinInd[nbSampleDuration/fs >= self.minimalSpindleDuration]
                
                newSpindles = [DetectedSpindle(channel, start, end) for start, end in zip(startSpinInd, stopSpinInd)]   
                self.detectedSpindles.extend(newSpindles)                   
                print "Nb detected spindles", len(self.detectedSpindles) 

 
 

 
 
 
 
 
 
class SpindleDectectorSigma(SpindleDectector):
    
    
    def __init__(self):
        SpindleDectector.__init__(self)
        self.detectStages = []
    

        ###############################################################################
        # Detection patameters
        ###############################################################################
    

        self.sigmaThreshold = 4.5

        # We allow for a spindle activity region to contain samples having an
        # amplitude going bellow the quantileThreshold for a 
        # maximum of maxAllowableGapBellowThreshold seconds.
        self.maxAllowableGapBellowThreshold = 0.1 # in seconds

        ###############################################################################


    def setThreshold(self, value): self.sigmaThreshold  = value

  
    # Detect every spindle in the channels listChannels of the file opened by the reader.
    def detectSpindles(self, reader, listChannels, verbose=True) :
        # List of detected spindles
        self.detectedSpindles = []
     
        #################################### READING ##########################
        if verbose:   print "Start reading datafile..."


        # Process only available channels...
        listChannels = [channel for channel in listChannels if channel in reader.getAvailableChannels()] 

        # Pickle data for each channel separatelly so simplify and accelerate
        # the reading of large files.       
        reader.pickleCompleteRecord(listChannels)   
   
        for channel in listChannels:    
            if verbose:   print "Channel " + channel + "..."           
            
            data        = reader.readPickledChannel(channel)

            signal      = data.signal
            fs          = data.samplingRate                      # sampling rate                     
                                                            

            ########################## Computing SIGMA INDEX ######################
            if verbose:   print "Computing sigma index..."
            
            fileName = reader.fname + "_sigmaIndex_" + channel + ".mat"
    
            if os.path.exists(fileName):
                print "Using saved sigma index..."    
                self.sigmaIndex = loadmat(fileName)["sigma"]
                self.sigmaIndex = self.sigmaIndex.reshape(self.sigmaIndex.size)
            else:     
                self.sigmaIndex = zeros(reader.getNbSample()) 
                        
                
                nbPad = int(0.1*fs)
                nbWin = int(4.0*fs)
    
                nbIter = int(reader.getNbSample()/nbWin)
                for i in range(nbIter):
                    if mod(i, 1000) == 0: print i, nbIter
                    
                    if i == 0 :            # First iteration
                        indexes = arange(nbPad + nbWin) 
                    elif i == nbIter-1 :   # Last iteration
                        indexes = arange(i*nbWin-nbPad, reader.getNbSample())
                    else:                       # Other iterations
                        indexes = arange(i*nbWin-nbPad, i*nbWin + nbWin+nbPad)
    
                    X, fX = computeMST(signal[indexes], fs, fmin=4.0, fmax=40.0)  
                    
                    if i == 0 :            # First iteration
                        indexesNoPad = arange(nbWin) 
                    elif i == nbIter-1 :   # Last iteration
                        indexesNoPad = arange(nbPad, nbPad + reader.getNbSample()-i*nbWin)
                    else:                       # Other iterations
                        indexesNoPad = arange(nbPad, nbPad + nbWin)
                        
                    X       = abs(X[indexesNoPad])
                    indexes = indexes[indexesNoPad]
        
                    maxalpha  = X[:, (fX >= 7.5)*(fX <= 10.0)].max(axis=1)
                    maxsigma  = X[:, (fX >= self.lowFreq)*(fX <= self.highFreq)].max(axis=1)
                    meanlow   = X[:, (fX >= 4.0)*(fX <= 10.0)].mean(axis=1) 
                    meanhigh  = X[:, (fX >= 20.0)*(fX <= 40.0)].mean(axis=1) 
    
                    sigmaTMP = array(2*maxsigma/(meanlow + meanhigh))
                    self.sigmaIndex[indexes[where(maxalpha <= maxsigma)[0]]] = sigmaTMP[where(maxalpha <= maxsigma)[0]]
 
                savemat(fileName, {"sigma":self.sigmaIndex})


            ############################### SPINDLE DETECTION #####################  
            if verbose:  print "Spindle detection..."
            
            spinMarkers = concatenate(([0], diff((self.sigmaIndex > self.sigmaThreshold).astype(int))))         
            
            # The signal indexes corresponding to spinMarkers are kept in samplesIndexes                
            
            startSpinInd = where(spinMarkers == 1)[0]
            stopSpinInd  = where(spinMarkers == -1)[0]
            
            if len(stopSpinInd) and len(startSpinInd):
                #The first marker should be a start marker.                
                if stopSpinInd[0] < startSpinInd[0]:
                    stopSpinInd = stopSpinInd[1:]

                indEnd = min(len(startSpinInd), len(stopSpinInd))
                
                startSpinInd = startSpinInd[:indEnd]
                stopSpinInd  = stopSpinInd[:indEnd] 

                #assert(all(stopSpinInd - startSpinInd > 0))
                #assert(all(startSpinInd[1:] - stopSpinInd[:-1] > 0))

                gapToKeepInd = where(startSpinInd[1:] - stopSpinInd[:-1] > self.maxAllowableGapBellowThreshold*fs)[0]
                
                startSpinInd = startSpinInd[concatenate(([0], gapToKeepInd+1))]
                stopSpinInd  = stopSpinInd[concatenate((gapToKeepInd, [len(stopSpinInd)-1]))]
                
                
                ## Careful : data are supporte contiguous. 
                nbSampleDuration = stopSpinInd - startSpinInd
                               
                startSpinInd = startSpinInd[nbSampleDuration/fs >= self.minimalSpindleDuration]
                stopSpinInd  = stopSpinInd[nbSampleDuration/fs >= self.minimalSpindleDuration]
                
                newSpindles = [DetectedSpindle(channel, start, end) for start, end in zip(startSpinInd, stopSpinInd)]   
                self.detectedSpindles.extend(newSpindles)                   
                print "Nb detected spindles", len(self.detectedSpindles) 

 
 
 
  
 
class SpindleDectectorRSP(SpindleDectector):
    
    
    def __init__(self):
        SpindleDectector.__init__(self)
        self.detectStages = []
    

        ###############################################################################
        # Detection patameters
        ###############################################################################
    

        self.threshold = 0.22

        # We allow for a spindle activity region to contain samples having an
        # amplitude going bellow the quantileThreshold for a 
        # maximum of maxAllowableGapBellowThreshold seconds.
        self.maxAllowableGapBellowThreshold = 0.0 # in seconds

        ###############################################################################



    def setThreshold(self, value): self.threshold  = value

  
    # Detect every spindle in the channels listChannels of the file opened by the reader.
    def detectSpindles(self, reader, listChannels, verbose=True) :
        # List of detected spindles
        self.detectedSpindles = []
     
        #################################### READING ##########################
        if verbose:   print "Start reading datafile..."


        # Process only available channels...
        listChannels = [channel for channel in listChannels if channel in reader.getAvailableChannels()] 

        # Pickle data for each channel separatelly so simplify and accelerate
        # the reading of large files.       
        reader.pickleCompleteRecord(listChannels)   
   
        for channel in listChannels:    
            if verbose:   print "Channel " + channel + "..."           
            
            data        = reader.readPickledChannel(channel)

            signal      = data.signal
            fs          = data.samplingRate                      # sampling rate                     
                                                            

            ########################## Computing SIGMA INDEX ######################
            if verbose:   print "Computing relative spindle power..."
            
            fileName = reader.fname + "_RSP_" + channel + ".mat"
    
            if os.path.exists(fileName):
                print "Using saved RSP..."    
                self.RSP = loadmat(fileName)["RSP"]
                self.RSP = self.RSP.reshape(self.RSP.size)
            else:     
                self.RSP = zeros(reader.getNbSample()) 
                        
                stageIndicator, lstStages = self.computeStageIndicator(reader)
                indDetectStage = concatenate([where(lstStages[stageIndicator] == stage.lower())[0] for stage in self.detectionStages])
                stageIndicator = zeros(stageIndicator.shape)
                stageIndicator[indDetectStage] = 1

                nbPad = int(0.1*fs)
                nbWin = int(4.0*fs)
    
                nbIter = int(reader.getNbSample()/nbWin)
                for i in range(nbIter):
                    if mod(i, 1000) == 0: print i, nbIter
                    
                    if i == 0 :            # First iteration
                        indexes = arange(nbPad + nbWin) 
                    elif i == nbIter-1 :   # Last iteration
                        indexes = arange(i*nbWin-nbPad, reader.getNbSample())
                    else:                       # Other iterations
                        indexes = arange(i*nbWin-nbPad, i*nbWin + nbWin+nbPad)
    
    
                    if any(stageIndicator[indexes]):    
                        X, fX = computeMST(signal[indexes], fs, fmin=0.5, fmax=40.0)  
                        
                        if i == 0 :            # First iteration
                            indexesNoPad = arange(nbWin) 
                        elif i == nbIter-1 :   # Last iteration
                            indexesNoPad = arange(nbPad, nbPad + reader.getNbSample()-i*nbWin)
                        else:                       # Other iterations
                            indexesNoPad = arange(nbPad, nbPad + nbWin)
                            
                        X       = abs(X[indexesNoPad])
                        indexes = indexes[indexesNoPad]
            
            
                        spindleBand = (fX >= self.lowFreq)*(fX <= self.highFreq)
                        self.RSP[indexes] = trapz(X[:, spindleBand], fX[:, spindleBand], axis=1)/trapz(X, fX, axis=1)
                    else:
                        self.RSP[indexes] = 0.0
 
                savemat(fileName, {"RSP":self.RSP})


            ############################### SPINDLE DETECTION #####################  
            if verbose:  print "Spindle detection..."
            
            spinMarkers = concatenate(([0], diff((self.RSP > self.threshold).astype(int))))         
            
            # The signal indexes corresponding to spinMarkers are kept in samplesIndexes                
            
            startSpinInd = where(spinMarkers == 1)[0]
            stopSpinInd  = where(spinMarkers == -1)[0]
            
            if len(stopSpinInd) and len(startSpinInd):
                #The first marker should be a start marker.                
                if stopSpinInd[0] < startSpinInd[0]:
                    stopSpinInd = stopSpinInd[1:]

                indEnd = min(len(startSpinInd), len(stopSpinInd))
                
                startSpinInd = startSpinInd[:indEnd]
                stopSpinInd  = stopSpinInd[:indEnd] 

                #assert(all(stopSpinInd - startSpinInd > 0))
                #assert(all(startSpinInd[1:] - stopSpinInd[:-1] > 0))

                gapToKeepInd = where(startSpinInd[1:] - stopSpinInd[:-1] > self.maxAllowableGapBellowThreshold*fs)[0]
                
                startSpinInd = startSpinInd[concatenate(([0], gapToKeepInd+1))]
                stopSpinInd  = stopSpinInd[concatenate((gapToKeepInd, [len(stopSpinInd)-1]))]
                
                
                ## Careful : data are supposed contiguous. 
                nbSampleDuration = stopSpinInd - startSpinInd
                               
                startSpinInd = startSpinInd[nbSampleDuration/fs >= self.minimalSpindleDuration]
                stopSpinInd  = stopSpinInd[nbSampleDuration/fs >= self.minimalSpindleDuration]
                
                newSpindles = [DetectedSpindle(channel, start, end) for start, end in zip(startSpinInd, stopSpinInd)]   
                self.detectedSpindles.extend(newSpindles)                
        