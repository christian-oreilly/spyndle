
"""
    Implementation of spindle detectors. 

    Copyright (C) 2013  Christian O'Reilly

    For personnal, educationnal, and research purpose, this software is 
    provided under the GNU GPL (V.3) license: you can redistribute it and/or
    modify it under the terms of the version 3 of the GNU General Public 
    License as published by the Free Software Foundation.
          
    To use this software in commercial application, please contact the author.
    
    If this code is used for research purpose, the reference [1] should be
    cited in the derived publication to refer the reader to the description 
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
     Exert Systems with Applications, Augustl 2013. 

"""

import os, gc, copy, numpy

from abc import ABCMeta, abstractmethod
import numpy as np
import warnings
import bisect

from scipy import concatenate, zeros, mean, sqrt, mod, diff, where, fft
from scipy import array, arange, ones, unique
from scipy.stats.mstats import mquantiles
from scipy.fftpack import fftfreq
from scipy.io import savemat, loadmat
from scipy.integrate import trapz

from datetime import timedelta

from spyndle import Filter
from spyndle import cycleDefinitions, computeDreamCycles
from spyndle import computeST
from spyndle.io import EEGDBReaderBase, Event
from spyndle.errorMng import ErrPureVirtualCall
from spyndle.EEG import getEEGChannels
from spyndle.io.databaseMng import buildTransientEvent

# Class representing a detected spindle.
class DetectedSpindle:
    def __init__(self, channel, startTime, endTime):

        # Channel on which spindle has been detected
        self.channel = channel
        
        # Spindle starting sample
        self.__startTime = startTime
        
        # Spindle ending sample
        self.__endTime   = endTime

        self.timeDuration   = endTime - startTime        

        self.RMSamp         = 0.0
        self.meanFreq       = 0.0
        self.sleepStage     = ""
        self.slopeOrigin    = 0.0
        self.slope          = 0.0
        self.filteredRMSamp = 0.0

    def startTime(self):
        return self.__startTime

    def endTime(self):
        return self.__endTime


    def computeRMS(self, reader, fmin=10, fmax=16):
        raise UserWarning("The code of this function need to be recoded to use time rather than sample.")        
        
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
        raise UserWarning("The code of this function need to be recoded to use time rather than sample.")        
        
        
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
        raise DeprecationWarning("No need anymore to call this function.")


    def setStage(self, reader):

        # Select the stage where the spindle begin as the sleep stage
        # of the spindle.
        stageEvent23 = filter(lambda e: e.groupName == "Stage" and
                                            e.timeStart() <=  self.startTime and 
                                            e.timeEnd() >= self.startTime, reader.events)                                 
        if len(stageEvent23) > 0 :
            self.sleepStage = stageEvent23[0].name.lower()
                
        else:
            print "No stage", self.startTime
            


# Generic spindle detector.
class SpindleDectector:
    __metaclass__ = ABCMeta    
    
    def __init__(self, reader=None, usePickled=False):
        ###############################################################################
        # Detection patameters
        ###############################################################################
        
        # Sleep stages in which we want to detect spindles. Should be a list of
        # event names used to score the stages we want to extract spindle from.        
        self.detectionStages    = ["Sleep stage 2", "Sleep stage N2"]
        
        # Low-end of the spindle frequency band
        self.lowFreq  = 11.0
        
        # High-end of the spindle frequency band        
        self.highFreq    = 16.0
        
        # Minimal duration of a valid spindle
        self.minimalSpindleDuration = 0.5 # in seconds

        # Maximal duration of a valid spindle. Avoid aberrant spindle extraction
        # such spindle with would last for tens of seconds because of a 
        # thresholding problem. It also avoid some problems related to the fact
        # that spindles are taken as being of short duration (e.g. a script
        # stoping crashing because it has exhausted all the memory making
        # a S-transform of a spindle of 30 seconds).
        self.maximalSpindleDuration = 3.0 # in seconds
        
        # Sleep cycle definition used to compute the quantile distribution
        # of amplitude per stage/cycles
        self.aeschbachCycleDef = cycleDefinitions()
        self.aeschbachCycleDef.setAeschbach()

        # List of detected spindles
        self.detectedSpindles = []

        # EEG data reader
        self.reader = reader

        # Pickling data allow faster spindle detection for some reader such
        # as the HarmonieReader but might not be implemented on other readers.
        # If self.usePickled, the spindle detection uses pickled data.
        self.usePickled = usePickled


        self.computeRMS         = False
        self.computeRMSFiltered = False
        self.computeFreq        = False       
        self.computeSlopeFreq   = False
        


    def setReader(self, reader):
        self.reader = reader


    """
     This function could be more elegantly written as :
         
        def simpleComputeRMS(self, reader, fmin=11, fmax=16):
            for spindle in self.detectedSpindles:        
                spindle.computeRMS(reader, fmin, fmax)
    
      but such an implementation is 40 times slower at execution than the 
      proposed version.
    """
    def computeRMS(self, fmin=11, fmax=16, ):
        channelList = unique([s.channel for s in self.detectedSpindles])        
        
        # Pickle data for each channel separatelly to simplify and accelerate
        # the reading of large files.    
        if self.usePickled :
            self.reader.pickleCompleteRecord(channelList)   
   
        for channel in channelList:    
            data        = self.reader.readChannel(channel, usePickled=self.usePickled)

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
        

        
        
        
    def computeMeanFreq(self, fmin=10, fmax=16):
        for spindle in self.detectedSpindles:        
            spindle.computeMeanFreq(self.reader, fmin, fmax)
        
    def computeTimeDuration(self):
        for spindle in self.detectedSpindles:        
            spindle.computeTimeDuration(self.reader)




    """
     This function could be more elegantly written as :
         
        def simpleSetStagesToSpindles(self, reader):
            for spindle in self.detectedSpindles:        
                spindle.setStage(reader)
        
      but such an implementation is 50 times slower at execution than the 
      proposed version.
    """    
    def setStagesToSpindles(self, reader):
        stageIndicator, lstStages = self.computeStageIndicator(reader)
        for spindle in self.detectedSpindles:        
            spindle.sleepStage = lstStages[stageIndicator[spindle.startSample]]



    def computeStageIndicator(self, reader, channel=None):

        stages = filter(lambda e: e.groupName == "Stage" , reader.events)   
        lstStage = concatenate((unique([stage.name.lower() for stage in stages]), ["No stage"]))    
        stageIndicator = ones(reader.getNbSample(channel))*(len(lstStage)-1)
            
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

        return concatenate((ones(windowNbSample/2)*result[0], result, 
                            ones(windowNbSample/2)* result[-1]))/windowNbSample 



    """
     Detect every spindle in the channels channelList of the file opened 
     by the reader.
    """
    @abstractmethod
    def detectSpindles(self, channelList=None, reader=None, verbose=True) : 
        if isinstance(channelList, str) :
            channelList = [channelList]

        if not reader is None:
            if isinstance(reader, EEGDBReaderBase):
                self.reader = reader
            else:
                raise TypeError
        else:
            if self.reader is None:
                raise ValueError("The reader attribute of the object is set to"
                                 " None and not reader object has been passed "
                                 "as argument to the detectSpindles() method.")
            

        if channelList is None:
            self.channelList = getEEGChannels(self.reader.getChannelLabels())            
        else:
            if isinstance(channelList, list) or isinstance(channelList, tuple):
                
                # Process only available channels to the reader...
                self.channelList = [channel for channel in channelList 
                                            if channel in self.reader.getChannelLabels()]                                
            else:
                raise TypeError            
            
        # List of detected spindles
        self.detectedSpindles = []

        
        
        

    # Used to save detected spindle in EEG data file. 
    def saveSpindle(self, reader, eventName, 
                    eventGroupName="Spindle", fileName = None, dbSession=None):
               
               
        for spindle in self.detectedSpindles:    
            event = Event(name          = eventName, 
                          groupName     = eventGroupName, 
                          channel       = spindle.channel, 
                          startTime     = spindle.startTime(),
                          timeLength    = spindle.timeDuration , 
                          dateTime      = reader.getRecordingStartTime() + timedelta(seconds=spindle.startTime()),
                          properties = {"RMSamp"            :spindle.RMSamp,
                                        "filteredRMSamp"    :spindle.filteredRMSamp,
                                        "meanFreq"          :spindle.meanFreq,
                                        "slopeOrigin"       :spindle.slopeOrigin,
                                        "slope"             :spindle.slope,
                                        "stage"             :spindle.sleepStage})     

            reader.addEvent(event)

            if dbSession :
                dbSession.add(buildTransientEvent(event, reader.fileName))

        if dbSession :
            dbSession.commit()
        
        if fileName is None:
            reader.save()
        else:      
            reader.saveAs(fileName)  
        

    # Used to save detected spindle in EEG data file. 
    def saveSpindleTxt(self, fileName):

        try:
            f = open(fileName, "w")    
        except IOError:     
            print "Error: The selected file could not be open."
            exit()
        
        for spindle in self.detectedSpindles:            
            f.write(spindle.channel + ";" + str(spindle.startTime()) + ";" +  
                    str(spindle.endTime()) + ";" + str(spindle.RMSamp) + ";" + str(spindle.meanFreq )
                      + ";" + str(spindle.timeDuration) + ";" + str(spindle.sleepStage) + "\n")             



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
        
        
    def computeFreqSlope_atDetection(self, signal, fs, startSpinInd, stopSpinInd, newSpindles, channelTime):
        
         for startInd, stopInd, spindle in zip(startSpinInd, stopSpinInd, newSpindles):
            sig       = signal[startInd:stopInd]
   
            X, fX = computeST(sig, fs, fmin=self.lowFreq-1, fmax=self.highFreq+1)  
            
            Y = abs(numpy.transpose(X))
                            
            regx = arange(len(sig))/fs
            regy = []
            try:
                for i in range(len(regx)):
                    regy.append(trapz(fX*Y[:, i], fX)/trapz(Y[:, i], fX))  
            except: 
                print fX.shape, X.shape, Y.shape, regx.shape, sig.shape
                print fs, self.lowFreq-1, self.highFreq+1, sig
                print channelTime[startInd:stopInd]
                raise

            z = numpy.polyfit(regx, regy, deg=1)     

            spindle.slopeOrigin = z[1]
            spindle.slope       = z[0]
 




class SpindleDectectorThreshold(SpindleDectector):
    __metaclass__ = ABCMeta    
        
    def __init__(self, reader=None, usePickled=False):
        SpindleDectector.__init__(self, reader, usePickled)
        
        ###############################################################################
        # Detection patameters
        ###############################################################################
    
        # We allow for a spindle activity region to contain samples having an
        # amplitude going bellow the quantileThreshold for a 
        # maximum of maxAllowableGapBellowThreshold seconds.
        self.maxAllowableGapBellowThreshold = 0.10 #0.07915697 #0.1 # in seconds
        

        # Is the threshold computed separetely for every self.perCycle and or
        # self.perStage.
        self.perCycle = True
        self.perStage = True

        ###############################################################################


    # Function performing any processing steps on the data to compute and return
    # the transformed signal that is an index of spindle presence.
    @abstractmethod   
    def preprocessing(self, data):
        raise ErrPureVirtualCall  
        
   
    ###################### THRESHOLDS #########################################
    # There is two kind of thresholds, the base threshold (self.__threshold)
    # and the effective threshold. Base threshold are a constant. Effective 
    # thresholds may depend on computation based on the value of 
    # the base threshold, the samplesIndexes of the
    # cycle and/or stage (depending on the value of self.perCycle and
    # self.perStage) and on the signal itself. The base threshold is defined
    # as a property such that it can be accessed using the syntax
    #    x = self.threshold
    #    self.threshold = x
    # Because effective thresholds need external parameters to get a proper 
    # value, the can only be properly defined using a method. By default,
    # the effective threshold takes the same value as the base threshold.
    # subclasses most redefine the getEffectiveThreshold(...) method to 
    # implement a more complexe behavior.
    def getEffectiveThreshold(self, samplesIndexes, signal):
        return self.threshold
    
    @property
    def threshold(self):
        return self.__threshold    
    
    @threshold.setter
    def threshold(self, value):
        self.__threshold = value        
        
    # For backward compatibility.
    def setThreshold(self, value): 
        warnings.warn("This method has been deprecated. self.threshold = value"\
                      " should be used instead.", DeprecationWarning )        
        self.threshold = value

  
        
  
  
    # Detect every spindle in the channels channelList of the file opened by the reader.
    def detectSpindles(self, channelList=None, reader=None, verbose=True) :
 
        SpindleDectector.detectSpindles(self, channelList, reader, verbose)

        # Computing sleep cycles
        if self.perCycle :
            cycles = computeDreamCycles([e for e in self.reader.events if e.groupName.lower() == "stage"], self.aeschbachCycleDef)

        #################################### READING ##########################
        if verbose:   print "Start reading datafile..."

        # Pickle data for each channel separatelly to simplify and accelerate
        # the reading of large files.    
        if self.usePickled :
            self.reader.pickleCompleteRecord(self.channelList)   
   
        for channel in self.channelList:    
            if verbose:   print "Channel " + channel + "..."           
            
            data            = self.reader.readChannel(channel, usePickled=self.usePickled)
            data.channel    = channel
            
            # We are working with two kind of signals, the raw EEG signal and
            # the transformed signal indexing the presence of sleep spindle.
            EEGsignal       = data.signal
            indexSignal     = self.preprocessing(data)

            ############################### SPINDLE DETECTION #####################  
            if verbose:  print "Spindle detection..."

            # Thereshold are computed for each sleep cycles and each sleep stage
            # since the average amplitude of EEG signal can vary accoss night
            # and stages.

            channelTime = self.reader.getChannelTime(channel)   
            assert(len(channelTime) == len(indexSignal))
            
            
            if self.perCycle :
                for cycle in cycles :
                    stageEvents = self.reader.getEventsByTime(cycle.timeStart(), cycle.timeEnd())               
                    stageEvents = filter(lambda e: e.groupName.lower() == "stage", stageEvents)     
                    
                    if self.perStage :
                        for stage in self.detectionStages :
                            currentStageEvents = filter(lambda e: e.name.lower() == stage.lower(), stageEvents)
                            self.__detectMain__(currentStageEvents, channelTime, data.samplingRate, 
                                                channel, indexSignal, EEGsignal, stage=stage)
                    else:
                        stageEvents = filter(lambda e: np.in1d(e.name.lower(),  [s.lower() for s in self.detectionStages]), stageEvents)
                        self.__detectMain__(stageEvents, channelTime, data.samplingRate, channel, indexSignal, EEGsignal)      
                        
                        
            else:         
                stageEvents = filter(lambda e: e.groupName.lower() == "stage", self.reader.events)                     
                if self.perStage :
                    for stage in self.detectionStages :
                        currentStageEvents = filter(lambda e: e.name.lower() == stage.lower(), stageEvents)
                        self.__detectMain__(currentStageEvents, channelTime, data.samplingRate, 
                                            channel, indexSignal, EEGsignal, stage=stage)
                else:
                    stageEvents = filter(lambda e: np.in1d(e.name.lower(),  [s.lower() for s in self.detectionStages]), stageEvents)
                    self.__detectMain__(stageEvents, channelTime, data.samplingRate, channel, indexSignal, EEGsignal)                       
                
                
                    




    # Performs the main processing steps involved in spindle detection.
    def __detectMain__(self, stageEvents, channelTime, fs, channel, indexSignal, EEGsignal, stage=None):

        if len(stageEvents) == 0:  return 
        
        #######################################################################
        ## This code....
        samplesIndexes = []
        indMin = 0
        for event in stageEvents:
            start  = bisect.bisect_left( channelTime[indMin:], event.timeStart())
            stop   = bisect.bisect_right(channelTime[indMin:], event.timeEnd())
            samplesIndexes.extend(range(start+indMin, stop+indMin))
            indMin += stop
        samplesIndexes = np.array(samplesIndexes)
        
        # ... replaced this one... 
        #samplesIndexes = concatenate([where((channelTime >= event.timeStart())*
        #                                    (channelTime <= event.timeEnd()))[0] 
        #                                                    for event in stageEvents])
        #sort(samplesIndexes) ...
        
        # ... to avoid trigerring this assert when events overlap ....
        assert(all(diff(samplesIndexes) >0))
        # ... but also to take benefits from the fact that channelTime is an
        # odered list. However, it have the down side to reallocating the 
        # samplesIndexes list many times. This shouldn't be too bad according
        # to http://stackoverflow.com/questions/5833907/repeatedly-appending-to-a-large-list-python-2-6-6
        # but a time benchmark between both approaches should be ran. 
        #######################################################################


        if len(samplesIndexes) == 0:  return
        
        spinMarkers = concatenate(([0], diff((indexSignal[samplesIndexes] > self.getEffectiveThreshold(samplesIndexes, indexSignal)).astype(int))))

        # The signal indexes corresponding to spinMarkers are kept in samplesIndexes                
        startSpinInd = samplesIndexes[where(spinMarkers == 1)[0]]
        stopSpinInd  = samplesIndexes[where(spinMarkers == -1)[0]]
              
              
        if len(stopSpinInd) == 0 or len(startSpinInd) == 0:  return
    
        #The first marker should be a start marker.                
        if stopSpinInd[0] < startSpinInd[0]:
            stopSpinInd = stopSpinInd[1:]

        indEnd = min(len(startSpinInd), len(stopSpinInd))
        
        startSpinInd = startSpinInd[:indEnd]
        stopSpinInd  = stopSpinInd[:indEnd] 

        if len(stopSpinInd) == 0 or len(startSpinInd) == 0:  return

        
        try:
            assert(all(stopSpinInd - startSpinInd > 0))
            assert(all(startSpinInd[1:] - stopSpinInd[:-1] > 0))
        except AssertionError:
            print startSpinInd, stopSpinInd, startSpinInd[1:] - stopSpinInd[:-1]
            raise

        gapToKeepInd = where(channelTime[startSpinInd[1:]] - channelTime[stopSpinInd[:-1]] > self.maxAllowableGapBellowThreshold)[0]
        
        startSpinInd = startSpinInd[concatenate(([0], gapToKeepInd+1))]
        stopSpinInd  = stopSpinInd[concatenate((gapToKeepInd, [len(stopSpinInd)-1]))]
        
        

        duration = channelTime[stopSpinInd] - channelTime[startSpinInd]
        ##Accept discontinuity of up to 3 time samples.
        continuous = [(max(diff(channelTime[start:stop])) <= 3.0/fs if stop - start > 10 else False ) 
                                                        for start, stop in zip(startSpinInd, stopSpinInd)]    
        
        valid        = (duration >= self.minimalSpindleDuration)*\
                       (duration <= self.maximalSpindleDuration)*(continuous)      
        startSpinInd = startSpinInd[valid]
        stopSpinInd  = stopSpinInd[valid]
        duration     = duration[valid]               
        
        newSpindles = [DetectedSpindle(channel, start, end) for start, end in 
                            zip(channelTime[startSpinInd], channelTime[stopSpinInd])]   
        
        #######################################################################
        # As the spindles are localized and the signal is extracted and filtred,
        # it is efficient to compute spindles characteristics directly here.
        
        # The RMS amplitude computed here is the raw RMS amplitude, that is the
        # amplitude computed on the unfiltered signal. It may be different than
        # the amplitude in the spindle band.
        if self.computeRMS :
            for startInd, stopInd, spindle in zip(startSpinInd, stopSpinInd, newSpindles):
                spindle.RMSamp = sqrt(mean(EEGsignal[startInd:stopInd]**2))
    
        if self.computeRMSFiltered :
   
            bandPassFilter = Filter(fs)
            bandPassFilter.create(low_crit_freq=self.lowFreq, 
                              high_crit_freq=self.highFreq, order=1001, 
                              btype="bandpass", ftype="FIR", useFiltFilt=True)     
            filteredEEGSignal = bandPassFilter.applyFilter(EEGsignal)

            for startInd, stopInd, spindle in zip(startSpinInd, stopSpinInd, newSpindles):
                spindle.filteredRMSamp= sqrt(mean(filteredEEGSignal[startInd:stopInd]**2))            
    
    
        if self.computeFreq :
            for startInd, stopInd, spindle in zip(startSpinInd, stopSpinInd, newSpindles):
                sig       = EEGsignal[startInd:stopInd]

                if sig.size < 512:
                    sig = concatenate((sig, zeros(512 - sig.size)))
                
                FFT = abs(fft(sig))
                freqs = fftfreq(sig.size, 1.0/fs)        
                
                FFT   = FFT[(freqs >= self.lowFreq)*(freqs <= self.highFreq)]
                freqs = freqs[(freqs >= self.lowFreq)*(freqs <= self.highFreq)]
                
                spindle.meanFreq = sum(freqs*FFT)/sum(FFT)
                                    
        
        if self.computeSlopeFreq :
            self.computeFreqSlope_atDetection(EEGsignal, fs, startSpinInd, 
                                              stopSpinInd, newSpindles, channelTime)

        #######################################################################

        if not stage is None:
            for spindle in newSpindles:
                spindle.sleepStage = stage

        self.detectedSpindles.extend(newSpindles)        


        
        

class SpindleDectectorRMS(SpindleDectectorThreshold):
    
    
    def __init__(self, reader=None, usePickled=False):
        SpindleDectectorThreshold.__init__(self, reader, usePickled)
        
        ###############################################################################
        # Detection patameters
        ###############################################################################
    
        # Quantile of spindle amplitude used to determine spindle activity
        self.threshold = 0.925

        # Width of the window used to computer the RMS amplitude.
        ### Choosing such that it always encompass at least two cycles of the 
        ### smaller frequency. We have a tradeoff between stability and accuracy...
        self.averagingWindowSize = 0.20  #0.508361 #2.0/lowFreq # In seconds

        ###############################################################################


    # Function performing any processing steps on the data to compute and return
    # the transformed signal that is an index of spindle presence.
    def preprocessing(self, data):

        # Defining EEG filters
        bandPassFilter = Filter(data.samplingRate)
        bandPassFilter.create(low_crit_freq=self.lowFreq, 
                              high_crit_freq=self.highFreq, order=1001, 
                              btype="bandpass", ftype="FIR", useFiltFilt=True)          

        # filtering can take a lot of memory. By making sure that the 
        # garbage collector as passed just before the filtering, we
        # increase our chances to avoid a MemoryError  
        gc.collect()                      
        signal = bandPassFilter.applyFilter(data.signal)     

        ################################# RMS COMPUTATION #####################
        windowNbSample = int(round(self.averagingWindowSize*data.samplingRate))
        if mod(windowNbSample, 2) == 0 : # We need an odd number.
            windowNbSample += 1

        # For selecting samples using a quantile-based thresholds, using abs(X)
        # or X**2 to rectify the X signal will give exactly the same result
        # since X**2 eqauls abs(X)**2 (for real numbers) and and the transformation
        # from abs(X) to abs(X)**2 is monotonically increasing, meaning that 
        # rank based statistics (such as quatiles) will give exactly the same
        # result. We use the numpy implementation of abs() because it is the
        # faster alternative.
        return self.averaging(numpy.abs(signal), windowNbSample)
    
    
    
    def getEffectiveThreshold(self, samplesIndexes, signal):
        return mquantiles(signal[samplesIndexes], self.threshold)[0]





# For backward compatibility.
SpindleDectectorAmp = SpindleDectectorRMS

 
 
 
class SpindleDectectorTeager(SpindleDectectorThreshold):
    
    
    def __init__(self, reader=None, usePickled=False):
        SpindleDectectorThreshold.__init__(self, reader, usePickled)
        
        ###############################################################################
        # Detection patameters
        ###############################################################################
        self.threshold = 0.6

        # We allow for a spindle activity region to contain samples having an
        # amplitude going bellow the quantileThreshold for a 
        # maximum of maxAllowableGapBellowThreshold seconds.
        self.maxAllowableGapBellowThreshold = 0.0 # in seconds
        ###############################################################################


    # Function performing any processing steps on the data to compute and return
    # the transformed signal that is an index of spindle presence.
    def preprocessing(self, data):

        # Defining EEG filters
        bandPassFilter = Filter(data.samplingRate)
        bandPassFilter.create(low_crit_freq=self.lowFreq, high_crit_freq=self.highFreq, 
                              order=1001, btype="bandpass", ftype="FIR", useFiltFilt=True)          
    
        # filtering can take a lot of memory. By making sure that the 
        # garbage collector as passed just before the filtering, we
        # increase our chances to avoid a MemoryError  
        gc.collect()                     
        signal = bandPassFilter.applyFilter(data.signal)                      
    
    
        ########################## Computing Teager operator ######################
        return concatenate(([0], signal[1:-1]**2 - signal[:-2]*signal[2:], [0]))
    
      
    def getEffectiveThreshold(self, samplesIndexes, signal):
        return mean(signal[samplesIndexes])*self.threshold


 
 
class SpindleDectectorSigma(SpindleDectectorThreshold):
    
    
    def __init__(self, reader=None, usePickled=False):
        SpindleDectectorThreshold.__init__(self, reader, usePickled)        

        ###############################################################################
        # Detection patameters
        ###############################################################################
        self.threshold = 4.5

        # We allow for a spindle activity region to contain samples having an
        # amplitude going bellow the quantileThreshold for a 
        # maximum of maxAllowableGapBellowThreshold seconds.
        self.maxAllowableGapBellowThreshold = 0.1 # in seconds

        # Duration of the window used to compute S-transform
        # (in seconds)
        self.computationWindow = 4.2   

        # Overlapping of two consecutive windows used to compute the S-transform
        # (in seconds)
        self.windowOverlapping = 0.2

        ###############################################################################



    # Function performing any processing steps on the data to compute and return
    # the transformed signal that is an index of spindle presence.
    def preprocessing(self, data):
    
        fileName = self.reader.fileName + "_sigmaIndex_" + data.channel   + ".mat"
    
        if os.path.exists(fileName):
            print "Using saved sigma index..."    
            self.sigmaIndex = loadmat(fileName)["sigma"]
            self.sigmaIndex = self.sigmaIndex.reshape(self.sigmaIndex.size)
        else:     
            self.sigmaIndex = zeros(self.reader.getNbSample(data.channel)) 
                    
            
            nbPad = int(self.windowOverlapping/2.0*data.samplingRate)
            nbWin = int(self.computationWindow*data.samplingRate) - 2*nbPad
    
            nbIter = int(self.reader.getNbSample(data.channel)/nbWin)
            for i in range(nbIter):
                if mod(i, 1000) == 0: print i, nbIter
                
                if i == 0 :            # First iteration
                    indexes = arange(nbPad + nbWin) 
                elif i == nbIter-1 :   # Last iteration
                    indexes = arange(i*nbWin-nbPad, self.reader.getNbSample(data.channel  ))
                else:                       # Other iterations
                    indexes = arange(i*nbWin-nbPad, i*nbWin + nbWin+nbPad)
    
                X, fX = computeST(data.signal[indexes], data.samplingRate, fmin=4.0, fmax=40.0)  
                
                if i == 0 :            # First iteration
                    indexesNoPad = arange(nbWin) 
                elif i == nbIter-1 :   # Last iteration
                    indexesNoPad = arange(nbPad, nbPad + self.reader.getNbSample(data.channel)-i*nbWin)
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
        return self.sigmaIndex
              

 
 
class SpindleDectectorRSP(SpindleDectectorThreshold):
    
    
    
    def __init__(self, reader=None, usePickled=False):
        SpindleDectectorThreshold.__init__(self, reader, usePickled)        
        
        ###############################################################################
        # Detection patameters
        ###############################################################################
        self.threshold = 0.22

        # We allow for a spindle activity region to contain samples having an
        # amplitude going bellow the quantileThreshold for a 
        # maximum of maxAllowableGapBellowThreshold seconds.
        self.maxAllowableGapBellowThreshold = 0.0 # in seconds

        ###############################################################################


    # Detect every spindle in the channels channelList of the file opened by the reader.
    def preprocessing(self, data):

        fileName = self.reader.getFileName() + "_RSP_" + data.channel + ".mat"

        if os.path.exists(fileName):
            print "Using saved RSP..."    
            self.RSP = loadmat(fileName)["RSP"]
            self.RSP = self.RSP.reshape(self.RSP.size)
        else:     
            self.RSP = zeros(self.reader.getNbSample(data.channel)) 

            nbPad = int(0.1*data.samplingRate)
            nbWin = int(4.0*data.samplingRate)

            nbIter = int(self.reader.getNbSample(data.channel)/nbWin)
            for i in range(nbIter):
                if mod(i, 1000) == 0: print i, nbIter
                
                if i == 0 :            # First iteration
                    indexes = arange(nbPad + nbWin) 
                elif i == nbIter-1 :   # Last iteration
                    indexes = arange(i*nbWin-nbPad, self.reader.getNbSample(data.channel))
                else:                       # Other iterations
                    indexes = arange(i*nbWin-nbPad, i*nbWin + nbWin+nbPad)


                #if any(stageIndicator[indexes]):    
                X, fX = computeST(data.signal[indexes], data.samplingRate, fmin=0.5, fmax=40.0)  
                
                if i == 0 :            # First iteration
                    indexesNoPad = arange(nbWin) 
                elif i == nbIter-1 :   # Last iteration
                    indexesNoPad = arange(nbPad, nbPad + self.reader.getNbSample(data.channel)-i*nbWin)
                else:                       # Other iterations
                    indexesNoPad = arange(nbPad, nbPad + nbWin)
                    
                X       = abs(X[indexesNoPad])
                indexes = indexes[indexesNoPad]
    
    
                spindleBand = (fX >= self.lowFreq)*(fX <= self.highFreq)
                self.RSP[indexes] = trapz(X[:, spindleBand], fX[:, spindleBand], axis=1)/trapz(X, fX, axis=1)
                #else:
                #    self.RSP[indexes] = 0.0
 
            savemat(fileName, {"RSP":self.RSP})
        return self.RSP
