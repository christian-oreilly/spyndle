
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

import gc

from scipy import concatenate, diff, where

from spyndle import Filter
from spyndle.detector import ThresholdDetector, DetectedEvent

import numpy as np
import bisect





 
class TeagerKComplexDetector(ThresholdDetector):
    """
     The Teager K-Complex detector, inspired by the work of Erdamar, Duman, & 
     Yetkin (2012).
    """
    
    def __init__(self, reader=None, usePickled=False):
        super(TeagerKComplexDetector, self).__init__(reader, usePickled)
        
        ###############################################################################
        # Detection patameters
        ###############################################################################
        self.threshold  = 0.6
        
        self.highFreq   = 5.0 # in Hz

        self.threshold  = 30
        
        self.t1Min      = 0.35
        self.t2Min      = 0.05 




    # Function performing any processing steps on the data to compute and return
    # the transformed signal that is an index of spindle presence.
    def preprocessing(self, data):

        # Defining EEG filters
        lowPassFilter = Filter(data.samplingRate)
        lowPassFilter.create(low_crit_freq=None, high_crit_freq=self.highFreq, 
                              order=1001, btype="lowpass", ftype="FIR", useFiltFilt=True)          
    
        # filtering can take a lot of memory. By making sure that the 
        # garbage collector as passed just before the filtering, we
        # increase our chances to avoid a MemoryError  
        gc.collect()                     
        signal = lowPassFilter.applyFilter(data.signal)                      
    
    
        ########################## Computing Teager operator ######################
        return concatenate(([0], signal[1:-1]**2 - signal[:-2]*signal[2:], [0]))
    
      


       
  
    """
    # Detect every event in the channels channelList of the file opened by the reader.
    def detectEvents(self, channelList=None, reader=None, verbose=True) :
 
        TransientDectector.detectEvents(self, channelList, reader, verbose)

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
                
                
                    

    """


    # Performs the main processing steps involved in  detection.
    def __detectMain__(self, stageEvents, channelTime, fs, 
                       channel, indexSignal, EEGsignal, stage=None):

        print("Detect K-Complexes")

        if len(stageEvents) == 0:  return 
        
        #######################################################################
        ## This code....
        samplesIndexes = []
        indMin = 0
        for event in stageEvents:
            start  = bisect.bisect_left( channelTime[indMin:], event.timeStart())
            stop   = bisect.bisect_right(channelTime[indMin:], event.timeEnd())
            samplesIndexes.extend(list(range(start+indMin, stop+indMin)))
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
        
        eventMarkersT2 = concatenate(([0], diff((indexSignal[samplesIndexes] > self.getEffectiveThreshold(samplesIndexes, indexSignal)).astype(int))))
        eventMarkersT1 = concatenate(([0], diff((indexSignal[samplesIndexes] > self.getEffectiveThreshold(samplesIndexes, indexSignal)).astype(int))))

        # The signal indexes corresponding to spinMarkers are kept in samplesIndexes                
        startIndT2 = samplesIndexes[where(eventMarkersT2 == 1)[0]]
        stopIndT2  = samplesIndexes[where(eventMarkersT2 == -1)[0]]
        
        startIndT1 = samplesIndexes[where(eventMarkersT1 == 1)[0]]
        stopIndT1  = samplesIndexes[where(eventMarkersT1 == -1)[0]]

        inds   = [(ind, 1) for ind in startIndT1] + \
                              [(ind, 2) for ind in startIndT2] + \
                              [(ind, 3) for ind in stopIndT2] + \
                              [(ind, 4) for ind in stopIndT1]
                        
        inds   = sorted(inds)
        
        newEvents = []
        for ind in range(len(inds)-3):
            if inds[ind][1] == 1 and inds[ind+1][1] == 2 and inds[ind+2][1] == 3 and inds[ind+3][1] == 4:
                if channelTime[inds[ind+3][0]] - channelTime[inds[ind][0]] >= self.t1Min:
                    if channelTime[inds[ind+2][0]] - channelTime[inds[ind][1]] >= self.t2Min:
                        newEvents.append(DetectedEvent(channel, 
                                                   channelTime[inds[ind][0]],
                                                   channelTime[inds[ind+3][0]]))
            
      
        
        #######################################################################
        # As the detected events are localized and the signal is extracted and 
        # filtred, it is efficient to compute events characteristics directly 
        # here. Such computation can be implemented by reimplementing the    
        # __postDetectionComputation__(...) method in subclasses.
        #self.postDetectionComputation(EEGsignal, channelTime, startInd, 
        #                                              stopInd, newEvents, fs)  
        

        if not stage is None:
            for spindle in newEvents:
                spindle.sleepStage = stage        
        
        self.detectedEvents.extend(newEvents)        


