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
from spyndle import CycleDefinitions
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
        
        # Event starting time
        self.__startTime = startTime
        
        # Event ending time
        self.__endTime   = endTime

    def __str__(self):
        return "{channel:" + str(self.channel)              \
               + ", startTime:" + str(self.__startTime)     \
               + ", endtime:" + str(self.__endTime) + "}"


    def duration(self):
        return self.timeDuration
        
    @property
    def timeDuration(self):
        return self.endTime() - self.startTime() 
        

    def startTime(self):
        return self.__startTime

    def endTime(self):
        return self.__endTime


    def setStage(self, reader):

        # Select the stage where the event begin as the sleep stage
        # of the event.
        stageEvent = [e for e in reader.events if e.groupName == "Stage" and
                                            e.timeStart() <=  self.startTime and 
                                            e.timeEnd() >= self.startTime]                                 
        if len(stageEvent) > 0 :
            self.sleepStage = stageEvent[0].name.lower()
                
        else:
            print(("No stage", self.startTime))
            


# Generic transient detector.
class TransientDetector(metaclass=ABCMeta):
    def __init__(self, reader=None, usePickled=False, verbose=False):

        self.__spyndle_version = __version__
        self.verbose = verbose

        ###############################################################################
        # Detection patameters
        ###############################################################################
        
        # Sleep stages in which we want to detect transient events. Should be a list of
        # event names used to score the stages we want to extract transient events from.        
        self.detectionStages    = ["Sleep stage 2", "Sleep stage N2"]
        
        # Sleep cycle definition used to compute the quantile distribution
        # of amplitude per stage/cycles
        self.aeschbachCycleDef = CycleDefinitions()
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
        return str({key: value for (key, value) in list(self.__dict__.items()) if key != 'detectedEvents'})


       
    def computeTimeDuration(self):
        for event in self.detectedEvents:        
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

        self.cycleDef = CycleDefinitions()
        self.cycleDef.setAeschbach()
        cycles = computeDreamCycles([e for e in reader.events if e.groupName.lower() == "stage"], self.cycleDef)
        
        nbEvents = 0
        for event in self.detectedEvents:
            if cycles[0].timeStart() > event.startTime():
                event.cycle = 1
                nbEvents += 1       
                
                
        for i, cycle in enumerate(cycles) :
            nbEvents = 0
            for event in self.detectedEvents:
                if cycle.timeStart() <= event.startTime() and cycle.timeEnd() > event.startTime():
                    event.cycle = i+1
                    nbEvents += 1

        nbEvents = 0
        for event in self.detectedEvents:
            if cycles[-1].timeEnd() <= event.startTime():
                event.cycle = len(cycles)
                nbEvents += 1


    def computeStageIndicator(self, reader, channel):

        stages = [e for e in reader.events if e.groupName.lower() == "stage"]   
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
    def computeRMS(self, fmin, fmax, reader=None):
        
        if reader is None:
            reader = self.reader        
        
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

            channelEvents = [s for s in self.detectedEvents if s.channel == channel]     
            for event in channelEvents:    
                indStart = bisect.bisect_left(t, event.startTime()) 
                indEnd   = indStart + int(event.timeDuration*fs)
                event.RMSAmp = np.sqrt(np.mean(signal[indStart:(indEnd+1)]**2))                
        


        


    def averaging(self, signal, windowNbSample, weightingFct=None):
     
        N = len(signal)-windowNbSample+1
        result = np.zeros(N)
        
        if weightingFct is None:        
            for i in range(windowNbSample) : 
                result += signal[i:i+N]
        else:
            for i, w in enumerate(weightingFct(windowNbSample)) :
                result += w*signal[i:i+N]

        # TODO: This line causes MemoryError in some case. This need to be fixed.
        return concatenate((ones(int(windowNbSample/2))*result[0], result, 
                            ones(int(windowNbSample/2))* result[-1]))/windowNbSample 



    """
     Detect every transient event in the channels channelList of the file opened 
     by the reader.
    """
    @abstractmethod
    def detectEvents(self, channelList=None, reader=None, verbose=None) : 
        
        if verbose is None:
            verbose = self.verbose
        
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
                    eventGroupName, fileName = None, dbSession=None, removeExisting=True):
              
              
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
          
        if removeExisting:   
            reader.events.removeType(eventName)
          
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
class ThresholdDetector(TransientDetector, metaclass=ABCMeta):
    def __init__(self, *args, **kwargs):

        if "excludedEventNames" in kwargs:
            self.excludedEventNames = kwargs["excludedEventNames"]
            del kwargs["excludedEventNames"]
        else:
            self.excludedEventNames = []
    
        super(ThresholdDetector, self).__init__(*args, **kwargs)
        
        ###############################################################################
        # Detection patameters
        ###############################################################################
    

        self.detectClass = DetectedEvent

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
    def preprocessing(self, signal, time=None, samplingRate=None, channel=None):
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
    def getEffectiveThreshold(self, signal):
        return self.threshold
    
    @property
    def threshold(self):
        return self.__threshold    
    
    @threshold.setter
    def threshold(self, value):
        self.__threshold = value        
        
        

    # Detect every event in the channels channelList of the file opened by the reader.
    def detectEvents(self, channelList=None, reader=None, verbose=None, retDetectStruc=False) :
        
        if verbose is None:
            verbose = self.verbose        
        
        TransientDetector.detectEvents(self, channelList, reader, verbose)
        
        detectStruct = []
        if self.perCycle :
            cycles = computeDreamCycles([e for e in self.reader.events if e.groupName.lower() == "stage"], self.aeschbachCycleDef)
            for i, cycle in enumerate(cycles) :
                if verbose: print(("Sleep cycle ", i+1))
                events = self.reader.getEventsByTime(cycle.timeStart(), cycle.timeEnd())               
                stageEvents = [e for e in events if e.groupName.lower() == "stage"]              
                excludedEvents = [e for e in events if e.name in self.excludedEventNames]       
                
                
                if self.perStage :
                    for stage in self.detectionStages :
                        if verbose: print(("Stage ", stage))
                        currentStageEvents = [e for e in stageEvents if e.name.lower() == stage.lower()]     
                        detecInfo = self._detectMain_(currentStageEvents, stage=stage, verbose=verbose, 
                                                                  excludedEvents=excludedEvents, retDetectStruc=retDetectStruc)
                        detectStruct.append([detecInfo, stage, i])
                else:
                    stageEvents = [e for e in stageEvents if np.in1d(e.name.lower(),  [s.lower() for s in self.detectionStages])]       
                    detecInfo = self._detectMain_(stageEvents, verbose=verbose, 
                                                              excludedEvents=excludedEvents, retDetectStruc=retDetectStruc) 
                    detectStruct.append([detecInfo, None, i])     
                    
        else:         
            stageEvents = [e for e in self.reader.events if e.groupName.lower() == "stage"]        
            excludedEvents = [e for e in self.reader.events if e.name in self.excludedEventNames]      
            if self.perStage :
                for stage in self.detectionStages :
                    if verbose: print(("Stage ", stage))                    
                    currentStageEvents = [e for e in stageEvents if e.name.lower() == stage.lower()]   
                    detecInfo = self._detectMain_(currentStageEvents, stage=stage, verbose=verbose, 
                                                              excludedEvents=excludedEvents, retDetectStruc=retDetectStruc)
                                                              
                    detectStruct.append([detecInfo, stage, None])
            else:
                stageEvents = [e for e in stageEvents if np.in1d(e.name.lower(),  [s.lower() for s in self.detectionStages])]                             
                detecInfo = self._detectMain_(stageEvents, verbose=verbose, 
                                              excludedEvents=excludedEvents, retDetectStruc=retDetectStruc)                       
            
                detectStruct.append([detecInfo, None, None])
          

        if retDetectStruc:
            return detectStruct



    # Performs the main processing steps involved in  detection.
    def _detectMain_(self, stageEvents, stage=None, DetectedClass=None, 
                     verbose=None, excludedEvents=[], retDetectStruc=False):

        if DetectedClass is None:
            DetectedClass = self.detectClass

        if verbose is None:
            verbose = self.verbose

        if len(stageEvents) == 0:  
            if verbose:
                print("No stageEvents. _detectMain_ is returning.")
            return 

        #################################### READING ##########################
        if verbose:   print("Start reading datafile...")

        # Pickle data for each channel separatelly to simplify and accelerate
        # the reading of large files.    
        if self.usePickled :
            self.reader.pickleCompleteRecord(self.channelList)   
   
        detecInfo = {"channel":[], "detectFct":[], "effectThresh":[]}
        for channel in self.channelList:    
            if verbose:   print(("Channel " + channel + "..."))      


            startTime = min([event.timeStart() for event in stageEvents])
            endTime   = max([event.timeEnd() for event in stageEvents])
            timeDuration = endTime-startTime

            #data            = self.reader.readChannel(channel, usePickled=self.usePickled)  
            data            = self.reader.read([channel], startTime, timeDuration)[channel]  
            channelTime     = self.reader.getChannelTime(channel, startTime=startTime, 
                                                         timeDuration=timeDuration)   
            assert(len(channelTime) == len(data.signal))

            samplesIndicators = self.reader.getEventIndicator(stageEvents, time=channelTime)  
            excludedIndicators = self.reader.getEventIndicator([e for e in excludedEvents if e.channel == channel], time=channelTime)  
            samplesIndicators *= np.logical_not(excludedIndicators) 
            
            if np.sum(samplesIndicators) == 0:  
                if verbose:
                    print("sample Indicators is empty. _detectMain_ is returning.")                
                return

            # We are working with two kind of signals, the raw EEG signal and
            # the transformed signal indexing the presence of the event to be detected.
            EEGsignal       = data.signal[samplesIndicators]
            channelTime     = channelTime[samplesIndicators]
            fs              = data.samplingRate
            indexSignal     = self.preprocessing(EEGsignal, channelTime, fs, channel)
            del data, samplesIndicators
            
            assert(len(EEGsignal) == len(indexSignal))

            ############################### DETECTION #####################  
            if verbose:  print("Detection...")

            # Thereshold are computed for each sleep cycles and each sleep stage
            # since the average amplitude of EEG signal can vary accoss night
            # and stages.       
            
            eventMarkers = concatenate(([0], diff((indexSignal > self.getEffectiveThreshold(indexSignal)).astype(int))))
    
            startInd = where(eventMarkers == 1)[0]
            stopInd  = where(eventMarkers == -1)[0]
                  
            if len(stopInd) == 0 or len(startInd) == 0:  
                if verbose:
                    print("stopInd or startId is empty. _detectMain_ is continuing with next channel.")       
                    print("Effective threshold: ", self.getEffectiveThreshold(indexSignal))
                    print("indexSignal, min and max: ", min(indexSignal), max(indexSignal))
                    print("len(indexSignal): ", len(indexSignal))
                    print("Length of startInd and stopInd: ", len(stopInd), len(startInd))
                continue
        
            #The first marker should be a start marker.                
            if stopInd[0] < startInd[0]:
                if verbose:
                    print("stopInd[0] < startInd[0]. Removing first stopInd")
                stopInd = stopInd[1:]

                if len(stopInd) == 0 or len(startInd) == 0:  
                    if verbose:
                        print("stopInd or startId is empty. _detectMain_ is continuing with next channel..")                
                    continue
        
    
            indEnd = min(len(startInd), len(stopInd))
            
            startInd = startInd[:indEnd]
            stopInd  = stopInd[:indEnd] 
    

            
            try:
                assert(all(stopInd - startInd > 0))
                assert(all(startInd[1:] - stopInd[:-1] > 0))
            except AssertionError:
                print((startInd, stopInd, startInd[1:] - stopInd[:-1]))
                raise
    
            # We remove gaps that are smaller than the threshold (self.maxAllowableGapBellowThreshold)
            gapToKeepInd = where(channelTime[startInd[1:]] - channelTime[stopInd[:-1]] > self.maxAllowableGapBellowThreshold)[0]
            startInd = startInd[concatenate(([0], gapToKeepInd+1))]
            stopInd  = stopInd[concatenate((gapToKeepInd, [len(stopInd)-1]))]
            
            
    
            duration = channelTime[stopInd] - channelTime[startInd]
            ##Accept discontinuity of up to 3 time samples.
            continuous = [(max(diff(channelTime[start:stop])) <= 3.0/fs if stop - start > 10 else False) 
                                                            for start, stop in zip(startInd, stopInd)]    
            
            valid        = (duration >= self.minimalDuration)*\
                           (duration <= self.maximalDuration)*(continuous)      
            startInd     = startInd[valid]
            stopInd      = stopInd[valid]
            duration     = duration[valid]               
            
            assert(np.all([end-start<= self.maximalDuration and end-start >= self.minimalDuration 
                           for start, end in zip(channelTime[startInd], channelTime[stopInd])]))
            
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
            
            self.detectedEvents.extend(newEvents)    
            
            if retDetectStruc:
                detecInfo["channel"].append(channel)
                detecInfo["detectFct"].append(indexSignal)
                detecInfo["effectThresh"].append(self.getEffectiveThreshold(indexSignal))

        if retDetectStruc:
            return detecInfo


    def postDetectionComputation(self, EEGsignal, channelTime, startInd, stopInd, newEvents, fs):
        pass

