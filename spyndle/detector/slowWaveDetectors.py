
"""
    Implementation of slow wave detectors. 

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

import gc, bisect

from scipy import where

from spyndle import Filter
from spyndle.detector import TransientDetector, DetectedEvent
from spyndle.miscellaneous.utils import findMaxima, findMinima, findZeroCrossings
from spyndle.io import Event, DataManipulationProcess, \
    PSGNight, Channel, SlowWaveEvent, EventClass

from datetime import timedelta
import numpy as np
import operator




 
class MassiminiSlowWaveDetector(TransientDetector):
    """
     The Slow-Wave detector, inspired by the work of Massimini, Huber,
     Ferrarelli, Hill and Tononi (2004).
    """
    
    def __init__(self, reader=None, usePickled=False):
        super(MassiminiSlowWaveDetector, self).__init__(reader, usePickled)
        
        ###############################################################################
        # Detection patameters
        ###############################################################################
        
        self.highFreq           = 4.0 # in Hz
        self.lowFreq            = 0.1 # in Hz

        self.negPeakMinDuration = 0.3 # in seconds
        self.negPeakMaxDuration = 1.0 # in seconds
        
        self.posPeakMaxDuration = 1.0 # in seconds
        
        self.minNegativePeakAmp = 80.0 # in uV         
        self.minPeakToPeakAmp   = 140.0 # in uV         
        


    def computeRMS(self, fmin=0.1, fmax=4.0):
        super(MassiminiSlowWaveDetector, self).computeRMS(self.lowFreq, self.highFreq)



    # Function performing any processing steps on the data to compute and return
    # the transformed signal that is an index of spindle presence.
    def preprocessing(self, data):

        # Defining EEG filters
        bandPassPassFilter = Filter(data.samplingRate)
        bandPassPassFilter.create(low_crit_freq=self.lowFreq, high_crit_freq=self.highFreq, 
                              order=1001, btype="bandpass", ftype="FIR", useFiltFilt=True)          
    
        # filtering can take a lot of memory. By making sure that the 
        # garbage collector as passed just before the filtering, we
        # increase our chances to avoid a MemoryError  
        gc.collect()                     
        return bandPassPassFilter.applyFilter(data.signal)                      
    



    # Detect every event in the channels channelList of the file opened by the reader.
    def detectEvents(self, channelList=None, reader=None, verbose=True) :
 
        TransientDetector.detectEvents(self, channelList, reader, verbose)

        #################################### READING ##########################
        if verbose:   print("Start reading datafile...")

        # Pickle data for each channel separatelly to simplify and accelerate
        # the reading of large files.    
        if self.usePickled :
            self.reader.pickleCompleteRecord(self.channelList)   
   
        for channel in self.channelList:    
            if verbose:   print(("Channel " + channel + "..."))           
            
            data            = self.reader.readChannel(channel, usePickled=self.usePickled)
            data.channel    = channel
            
            # We are working with two kind of signals, the raw EEG signal and
            # the transformed signal indexing the presence of the event to be detected.
            signal     = self.preprocessing(data)

            ############################### DETECTION #####################  
            if verbose:  print("Detection...")

            # Thereshold are computed for each sleep cycles and each sleep stage
            # since the average amplitude of EEG signal can vary accoss night
            # and stages.

            channelTime = self.reader.getChannelTime(channel)   
            assert(len(channelTime) == len(signal))
            
            
            stageEvents = [e for e in self.reader.events if e.groupName.lower() == "stage"]                
            stageEvents = [e for e in stageEvents if np.in1d(e.name.lower(),  [s.lower() for s in self.detectionStages])]
            self._detectMain_(stageEvents, channelTime, data.samplingRate, channel, signal)                       
                
        # These are performed for all detected events, regarless of their channels.
        self.setStagesToEvents()
        self.setCyclesToEvents()

           
          
          

    # Performs the main processing steps involved in  detection.
    def _detectMain_(self, stageEvents, channelTime, fs, 
                       channel, signal):


        def validateSeq(seqInd, channelTime, labelInds, signal):
            
            dt = channelTime[labelInds[seqInd+2]] - channelTime[labelInds[seqInd]]  

            if dt < self.negPeakMinDuration:
                return False
                
            if dt > self.negPeakMaxDuration:
                return False
            

            dt = channelTime[labelInds[seqInd+4]] - channelTime[labelInds[seqInd+2]]              
            if dt > self.posPeakMaxDuration:
                return False            
        
            if signal[labelInds[seqInd+3]] - signal[labelInds[seqInd+1]] < self.minPeakToPeakAmp :
                return False
                
            return True
    
    
    
        def setToEvent(seqInd, channelTime, labelInds, channel):    
           return DetectedEvent(channel, channelTime[labelInds[seqInd]], 
                                          channelTime[labelInds[seqInd+4]])



        def removeDuplicateMax(labels, labelInds, signal):
            
            indToRemove = []
            best = -np.inf
            indBest = 0
            
            for i, (label, value) in enumerate(zip(labels, signal[labelInds])):
                if label == 1:
                    if value > best:
                        best = value
                        if indBest != 0:
                            indToRemove.append(indBest)
                        indBest = i
                    else:
                        indToRemove.append(i)
                else:
                    best = -np.inf
                    indBest = 0
                    
            return np.delete(labels, indToRemove),  np.delete(labelInds, indToRemove)
                        


        def removeDuplicateMin(labels, labelInds, signal):
            
            indToRemove = []
            best = np.inf
            indBest = 0
            
            for i, (label, value) in enumerate(zip(labels, signal[labelInds])):            
                if label == -1:
                    if value < best:
                        best = value
                        if indBest != 0:
                            indToRemove.append(indBest)
                        indBest = i
                    else:
                        indToRemove.append(i)
                else:
                    best = np.inf
                    indBest = 0
                    
            return np.delete(labels, indToRemove),  np.delete(labelInds, indToRemove)
                        




        zeroCrossings = findZeroCrossings(signal)
        minPeaks      = findMinima(signal)
        maxPeaks      = findMaxima(signal)

        minPeaks = minPeaks[where(signal[minPeaks] < -self.minNegativePeakAmp)]

        
        labelInds    = np.concatenate((minPeaks, zeroCrossings, maxPeaks))
        labels      = np.concatenate((-np.ones(len(minPeaks)), 
                                      np.zeros(len(zeroCrossings)), 
                                      np.ones(len(maxPeaks))))




        labelInds, labels = list(zip(*sorted(zip(labelInds, labels), key=operator.itemgetter(0))))
        labels      = np.array(labels)
        labelInds   = np.array(labelInds)

        labels, labelInds = removeDuplicateMax(labels, labelInds, signal)
        try:        
            assert(not np.any(labels[:-1]*labels[1:] == 1))
        except AssertionError:
            ValueError(str(labels[:-1]*labels[1:]))

        labels, labelInds = removeDuplicateMin(labels, labelInds, signal)
        assert(not np.any((labels[0:-1] == -1)*(labels[1:] == -1)))


        # We are looking for (0, -1, 0, 1, 0) label sequences.
        seqInds =  np.where(  (labels[ :-4] == 0)*\
                              (labels[1:-3] == -1)*\
                              (labels[2:-2] == 0)*\
                              (labels[3:-1] == 1)*\
                              (labels[4:  ] == 0))[0]
          
        newEvents = [setToEvent(seqInd, channelTime, labelInds, channel) 
                            for seqInd in seqInds 
                            if validateSeq(seqInd, channelTime, labelInds, signal)]
                
     
        self.postDetectionComputation(signal, channelTime, labelInds, seqInds, newEvents, channel)
        
        self.detectedEvents.extend(newEvents)              











    def postDetectionComputation(self, signal, channelTime, labelInds, seqInds, newEvents, channel) :
                                                          

        fs = 1.0/(channelTime[1]- channelTime[0])
        for newEvent, seqInd in zip(newEvents, seqInds):
    
            newEvent.negativeDuration   = channelTime[labelInds[seqInd+2]] - channelTime[labelInds[seqInd]]
            newEvent.positiveDuration   = channelTime[labelInds[seqInd+4]] - channelTime[labelInds[seqInd+2]] 
            newEvent.ZNSlope            = signal[labelInds[seqInd+1]]/(channelTime[labelInds[seqInd+1]] - channelTime[labelInds[seqInd]])
            newEvent.NPSlope            = (signal[labelInds[seqInd+3]] - signal[labelInds[seqInd+1]])/\
                                          (channelTime[labelInds[seqInd+3]] - channelTime[labelInds[seqInd+1]])
            newEvent.PZSlope            = signal[labelInds[seqInd+3]]/(channelTime[labelInds[seqInd+4]] - channelTime[labelInds[seqInd+3]])
            newEvent.timeMin            = channelTime[labelInds[seqInd+1]]
            newEvent.timeMax            = channelTime[labelInds[seqInd+3]]
            
            #### RMS computation
            # We don't use self.computeRMS() because we want to take benifits from
            # the fact that the signal variable holds the filtered EEG such that 
            # it don't need to be filtered again for RMS computation              
            indStart = bisect.bisect_left(channelTime, newEvent.startTime()) 
            indEnd   = indStart + int(newEvent.timeDuration*fs)
            newEvent.RMSAmp = np.sqrt(np.mean(signal[indStart:(indEnd+1)]**2))                
            
              
            
            





        





    # Used to save detected events in EEG data file. 
    def saveEvents(self, reader, eventName="MassiminiSL", 
                    eventGroupName="SlowWave", fileName = None, dbSession=None):
              
        # TODO: This code should not repeat the code of TransientDetector.saveEvents(...)              
              
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
            for channel in np.unique([e.channel for e in self.detectedEvents]):
                if dbSession.query(Channel).filter_by(name=channel).count() == 0:
                    dbSession.add(Channel(name=channel))                    
            dbSession.flush()                          
             


        if dbSession :
            transientEvents = []
            slowWaveEvents   = []
          
        for detectedEvent in self.detectedEvents:    
            event = Event(name          = eventName, 
                          groupName     = eventGroupName, 
                          channel       = detectedEvent.channel, 
                          startTime     = detectedEvent.startTime(),
                          timeLength    = detectedEvent.timeDuration , 
                          dateTime      = reader.getRecordingStartTime() + timedelta(seconds=detectedEvent.startTime()),
                          properties = {"negativeDuration"  :detectedEvent.negativeDuration,
                                        "positiveDuration"  :detectedEvent.positiveDuration,  
                                        "stage"             :detectedEvent.sleepStage,
                                        "RMSAmp"            :detectedEvent.RMSAmp,
                                        "cycle"             :detectedEvent.cycle,
                                        "ZNSlope"           :detectedEvent.ZNSlope,
                                        "NPSlope"           :detectedEvent.NPSlope,
                                        "PZSlope"           :detectedEvent.PZSlope,
                                        "timeMin"           :detectedEvent.timeMin,
                                        "timeMax"           :detectedEvent.timeMax}) 

            reader.addEvent(event)

            if dbSession :
                transient, sl = SlowWaveEvent.fromEvent(event, reader.fileName, dataManipObj.no)
                transientEvents.append(transient)
                slowWaveEvents.append(sl)
                
                if len(transientEvents) > 100:
                    dbSession.add_all(transientEvents)
                    dbSession.flush()
                    dbSession.add_all(slowWaveEvents)
                    dbSession.flush()                
                    transientEvents = []
                    slowWaveEvents   = []                

        if dbSession :
            dbSession.add_all(transientEvents)
            dbSession.flush()
            dbSession.add_all(slowWaveEvents)
            dbSession.flush()
        
        if fileName is None:
            reader.save()
        else:      
            reader.saveAs(fileName)  
            
        



