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
 Date  : January 13, 2014

 [1] O'Reilly, C. & Nielsen, T. Sleep spindle detection: Automation and 
     performance evaluation using fine temporal resolution, Submitted to
     Exert Systems with Applications, Augustl 2013. 

"""

import copy, gc

from abc import ABCMeta, abstractmethod
import numpy as np


from scipy import concatenate, ones, unique, diff, where

from datetime import timedelta
import bisect

from spyndle import __version__
from spyndle import cycleDefinitions
from spyndle.EEG import getEEGChannels
from spyndle import computeDreamCycles
from spyndle.io import EEGDBReaderBase, Event, DataManipulationProcess, \
    PSGNight, Channel, TransientEvent, EventClass

from spyndle import Filter





# Class representing a detected event.
class DetectedEvent(object):
    def __init__(self, channel, startTime, endTime):

        # Channel on which event has been detected
        self.channel = channel
        
        # Event starting sample
        self.__startTime = startTime
        
        # Event ending sample
        self.__endTime   = endTime

        self.timeDuration   = endTime - startTime        


    def startTime(self):
        return self.__startTime

    def endTime(self):
        return self.__endTime


    def setStage(self, reader):

        # Select the stage where the event begin as the sleep stage
        # of the event.
        stageEvent23 = filter(lambda e: e.groupName == "Stage" and
                                            e.timeStart() <=  self.startTime and 
                                            e.timeEnd() >= self.startTime, reader.events)                                 
        if len(stageEvent23) > 0 :
            self.sleepStage = stageEvent23[0].name.lower()
                
        else:
            print "No stage", self.startTime
            


# Generic transient detector.
class TransientDetector:
    __metaclass__ = ABCMeta    
    
    def __init__(self, reader=None, usePickled=False):

        self.__spyndle_version = __version__

        ###############################################################################
        # Detection patameters
        ###############################################################################
        
        # Sleep stages in which we want to detect transient events. Should be a list of
        # event names used to score the stages we want to extract transient events from.        
        self.detectionStages    = ["Sleep stage 2", "Sleep stage N2"]
        
        # Sleep cycle definition used to compute the quantile distribution
        # of amplitude per stage/cycles
        self.aeschbachCycleDef = cycleDefinitions()
        self.aeschbachCycleDef.setAeschbach()

        # List of detected events
        self.detectedEvents = []

        # Minimal duration of a valid event
        self.minimalDuration = 0.0 # in seconds

        # Maximal duration of a valid event. 
        self.maximalDuration = np.inf # in seconds


        # EEG data reader
        self.reader = reader

        # Pickling data allow faster detection for some reader such
        # as the HarmonieReader but might not be implemented on other readers.
        # If self.usePickled, the event detection uses pickled data.
        self.usePickled = usePickled

        
        
        
    def __repr__(self):
        # serpent does not manage a whole lot of types such as numpy types. 
        # It look premature to use it for now.
        #return serpent.dumps(self, indent=False, set_literals=False)
        # TODO: Improve on this representation...
    
        # We keep every fields of the object, except for the detectedEvents 
        # which would create a much too large representation. Moreover, 
        # the detected events are stored as an attribute only for 
        # convenience. They do not define the detector per se. 
        return str({key: value for (key, value) in self.__dict__.iteritems() if key != 'detectedEvents'})


       
    def computeTimeDuration(self):
        for event in self.detectedEventss:        
            event.computeTimeDuration(self.reader)




    """
     This function could be more elegantly written as :
         
        def simpleSetStagesToEvents(self, reader):
            for spindle in self.detectedEvents:        
                spindle.setStage(reader)
        
      but such an implementation is 50 times slower at execution than the 
      proposed version.
    """    
    def setStagesToEvents(self, reader=None):
        if reader is None:
            reader = self.reader
        
        stageIndicator, lstStages = self.computeStageIndicator(reader, self.detectedEvents[0].channel)
        t = reader.getChannelTime(self.detectedEvents[0].channel)
        for event in self.detectedEvents:        
            indStart = bisect.bisect_left(t, event.startTime()) 
            event.sleepStage = lstStages[stageIndicator[indStart]]


    def setCyclesToEvents(self, reader = None):
        if reader is None:
            reader = self.reader

        self.cycleDef = cycleDefinitions()
        self.cycleDef.setAeschbach()
        cycles = computeDreamCycles([e for e in reader.events if e.groupName.lower() == "stage"], self.cycleDef)
        
        for event in self.detectedEvents:
            if cycles[0].timeStart() > event.startTime:
                event.cycle = 1
                
        for i, cycle in enumerate(cycles) :
            for event in self.detectedEvents:
                if cycle.timeStart() <= event.startTime and cycle.timeEnd() > event.startTime:
                    event.cycle = i+1

        for event in self.detectedEvents:
            if cycles[-1].timeEnd() <= event.startTime:
                event.cycle = len(cycles)




    def computeStageIndicator(self, reader, channel):

        stages = filter(lambda e: e.groupName.lower() == "stage" , reader.events)   
        lstStage = list(concatenate((unique([str(stage.name) for stage in stages]), ["No stage"])))    
        stageIndicator = ones(reader.getNbSample(channel))*(len(lstStage)-1)
            
        for i in range(len(lstStage)):
            index = where(self.reader.getEventIndicator(lstStage[i], channel, globalEvent=True))[0]
            stageIndicator[index] = i
                
        return stageIndicator.astype(int), lstStage




    """
     This function could be more elegantly written as :
         
        def computeRMS(self, reader, fmin=11, fmax=16):
            for event in self.detectedSpindles:        
                event.computeRMS(reader, fmin, fmax)
    
      but such an implementation is 40 times slower at execution than the 
      proposed version.
    """
    def computeRMS(self, fmin, fmax):
        channelList = unique([s.channel for s in self.detectedEvents])        
        
        # Pickle data for each channel separatelly to simplify and accelerate
        # the reading of large files.    
        if self.usePickled :
            self.reader.pickleCompleteRecord(channelList)   
   
        for channel in channelList:    
            data        = self.reader.readChannel(channel, usePickled=self.usePickled)

            signal      = data.signal
            fs          = data.samplingRate                      # sampling rate    
            t           = self.reader.getChannelTime(self.detectedEvents[0].channel)
                    
            # Defining EEG filters   
            bandPassFilter = Filter(fs)
            bandPassFilter.create(low_crit_freq=fmin, high_crit_freq=fmax, order=1001, btype="bandpass", ftype="FIR", useFiltFilt=True)          
       
                   
            # filtering can take a lot of memory. By making sure that the 
            # garbage collector as passed just before the filtering, we
            # increase our chances to avoid a MemoryError  
            gc.collect()                     
            signal = bandPassFilter.applyFilter(signal)             

            channelEvents = filter(lambda s: s.channel == channel, self.detectedEvents)     
            for event in channelEvents:    
                indStart = bisect.bisect_left(t, event.startTime()) 
                indEnd   = indStart + int(event.timeDuration*fs)
                event.RMSamp = np.sqrt(np.mean(signal[indStart:(indEnd+1)]**2))                
        


        


    def averaging(self, signal, windowNbSample):
     
        result = copy.copy(signal[:(len(signal)-windowNbSample+1)])
        
        for i in range(1, windowNbSample) : 
            result += signal[i:(i+len(signal)-windowNbSample+1)]

        # TODO: This line causes MemoryError in some case. THis need to be fixed.
        return concatenate((ones(windowNbSample/2)*result[0], result, 
                            ones(windowNbSample/2)* result[-1]))/windowNbSample 



    """
     Detect every transient event in the channels channelList of the file opened 
     by the reader.
    """
    @abstractmethod
    def detectEvents(self, channelList=None, reader=None, verbose=True) : 
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
                                 "as argument to the detectEvents() method.")
            

        if channelList is None:
            self.channelList = getEEGChannels(self.reader.getChannelLabels())            
        else:
            if isinstance(channelList, list) or isinstance(channelList, tuple):
                
                # Process only available channels to the reader...
                self.channelList = [channel for channel in channelList 
                                            if channel in self.reader.getChannelLabels()]                                
            else:
                raise TypeError            
            
        # List of detected events
        self.detectedEvents = []

        
        
        

    # Used to save detected events in EEG data file. 
    def saveEvents(self, reader, eventName, 
                    eventGroupName, fileName = None, dbSession=None):
              
        if dbSession :
            
            # Create the data manipulation process record.
            dataManipObj = DataManipulationProcess(reprStr  = repr(self))
            dbSession.add(dataManipObj)                             

            # Create the event class record if none corresponding to this event class exist.                                       
            if dbSession.query(EventClass).filter_by(name=eventName).count() == 0:
                dbSession.add(EventClass(name=eventName))     

            # Create the PSG night record if none corresponding to this night exist. 
            if dbSession.query(PSGNight).filter_by(fileName=reader.fileName).count() == 0:
                dbSession.add(PSGNight(fileName=reader.fileName))                            

            # Create the channel record if none corresponding to this channel exist.                    
            for channel in unique([e.channel for e in self.detectedEvents]):
                if dbSession.query(Channel).filter_by(name=channel).count() == 0:
                    dbSession.add(Channel(name=channel))                    
            dbSession.flush()                          
             


        if dbSession :
            transientEvents = []
          
        for detectedEvent in self.detectedEvents:    
            event = Event(name          = eventName, 
                          groupName     = eventGroupName, 
                          channel       = detectedEvent.channel, 
                          startTime     = detectedEvent.startTime(),
                          timeLength    = detectedEvent.timeDuration , 
                          dateTime      = reader.getRecordingStartTime() + timedelta(seconds=detectedEvent.startTime()),
                          properties = {"stage"             :detectedEvent.sleepStage})     

            reader.addEvent(event)

            if dbSession :
                transientEvents.append(TransientEvent.fromEvent(event, reader.fileName, dataManipObj.no))

        if dbSession :
            dbSession.add_all(transientEvents)

        
        if fileName is None:
            reader.save()
        else:      
            reader.saveAs(fileName)  
            
        
  
    # Set the sleep stages in which we want to detect events. stages should
    # be a list of valid stages.
    def setDetectionStages(self, stages):
        self.detectionStages = stages
        
        
        
        
        
        
        





"""
 All event detectors using some variable bearing information on
 the presence of events to be detected can be thought using a same general 
 framework. This class implement this general framework. Threshold-based 
 detector can be implemented by subclassing this classe.
"""
class ThresholdDetector(TransientDetector):
    __metaclass__ = ABCMeta    
        
    def __init__(self, reader=None, usePickled=False):
        super(ThresholdDetector, self).__init__(reader, usePickled)
        
        ###############################################################################
        # Detection patameters
        ###############################################################################
    
        # We allow for a activity region to contain samples having an
        # amplitude going bellow the threshold for a 
        # maximum of maxAllowableGapBellowThreshold seconds.
        self.maxAllowableGapBellowThreshold = 0.10 # in seconds
        

        # Is the threshold computed separetely for every self.perCycle and or
        # self.perStage.
        self.perCycle = True
        self.perStage = True

        ###############################################################################


    # Function performing any processing steps on the data to compute and return
    # the transformed signal that is an index of transient event presence.
    @abstractmethod   
    def preprocessing(self, data):
        raise NotImplementedError("This method must be implemented in derived classes.")  
        
   
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
        
        
  
  
    # Detect every event in the channels channelList of the file opened by the reader.
    def detectEvents(self, channelList=None, reader=None, verbose=True) :
 
        TransientDetector.detectEvents(self, channelList, reader, verbose)

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
            # the transformed signal indexing the presence of the event to be detected.
            EEGsignal       = data.signal
            indexSignal     = self.preprocessing(data)

            ############################### DETECTION #####################  
            if verbose:  print "Detection..."

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
                            self._detectMain_(currentStageEvents, channelTime, data.samplingRate, 
                                                channel, indexSignal, EEGsignal, stage=stage)
                    else:
                        stageEvents = filter(lambda e: np.in1d(e.name.lower(),  [s.lower() for s in self.detectionStages]), stageEvents)
                        self._detectMain_(stageEvents, channelTime, data.samplingRate, channel, indexSignal, EEGsignal)      
                        
                        
            else:         
                stageEvents = filter(lambda e: e.groupName.lower() == "stage", self.reader.events)                     
                if self.perStage :
                    for stage in self.detectionStages :
                        currentStageEvents = filter(lambda e: e.name.lower() == stage.lower(), stageEvents)
                        self._detectMain_(currentStageEvents, channelTime, data.samplingRate, 
                                            channel, indexSignal, EEGsignal, stage=stage)
                else:
                    stageEvents = filter(lambda e: np.in1d(e.name.lower(),  [s.lower() for s in self.detectionStages]), stageEvents)
                    self._detectMain_(stageEvents, channelTime, data.samplingRate, channel, indexSignal, EEGsignal)                       
                
                
                    

    


    # Performs the main processing steps involved in  detection.
    def _detectMain_(self, stageEvents, channelTime, fs, 
                       channel, indexSignal, EEGsignal, stage=None, DetectedClass=DetectedEvent):

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
        
        eventMarkers = concatenate(([0], diff((indexSignal[samplesIndexes] > self.getEffectiveThreshold(samplesIndexes, indexSignal)).astype(int))))

        # The signal indexes corresponding to spinMarkers are kept in samplesIndexes                
        startInd = samplesIndexes[where(eventMarkers == 1)[0]]
        stopInd  = samplesIndexes[where(eventMarkers == -1)[0]]
              
              
        if len(stopInd) == 0 or len(startInd) == 0:  return
    
        #The first marker should be a start marker.                
        if stopInd[0] < startInd[0]:
            stopInd = stopInd[1:]

        indEnd = min(len(startInd), len(stopInd))
        
        startInd = startInd[:indEnd]
        stopInd  = stopInd[:indEnd] 

        if len(stopInd) == 0 or len(startInd) == 0:  return

        
        try:
            assert(all(stopInd - startInd > 0))
            assert(all(startInd[1:] - stopInd[:-1] > 0))
        except AssertionError:
            print startInd, stopInd, startInd[1:] - stopInd[:-1]
            raise

        gapToKeepInd = where(channelTime[startInd[1:]] - channelTime[stopInd[:-1]] > self.maxAllowableGapBellowThreshold)[0]
        
        startInd = startInd[concatenate(([0], gapToKeepInd+1))]
        stopInd  = stopInd[concatenate((gapToKeepInd, [len(stopInd)-1]))]
        
        

        duration = channelTime[stopInd] - channelTime[startInd]
        ##Accept discontinuity of up to 3 time samples.
        continuous = [(max(diff(channelTime[start:stop])) <= 3.0/fs if stop - start > 10 else False ) 
                                                        for start, stop in zip(startInd, stopInd)]    
        
        valid        = (duration >= self.minimalDuration)*\
                       (duration <= self.maximalDuration)*(continuous)      
        startInd     = startInd[valid]
        stopInd      = stopInd[valid]
        duration     = duration[valid]               
        
        newEvents = [DetectedClass(channel, start, end) for start, end in 
                            zip(channelTime[startInd], channelTime[stopInd])]   
      
        
        #######################################################################
        # As the detected events are localized and the signal is extracted and 
        # filtred, it is efficient to compute events characteristics directly 
        # here. Such computation can be implemented by reimplementing the    
        # __postDetectionComputation__(...) method in subclasses.
        self.postDetectionComputation(EEGsignal, channelTime, startInd, 
                                                      stopInd, newEvents, fs)  
        

        if not stage is None:
            for spindle in newEvents:
                spindle.sleepStage = stage        
        
        self.detectedSpindles.extend(newEvents)        



    def postDetectionComputation(self, startInd, stopInd, newEvents, fs):
        pass

