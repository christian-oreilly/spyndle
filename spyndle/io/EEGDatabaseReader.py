# -*- coding: utf-8 -*-

###############################################################################
# License: For personnal, educationnal, and research purpose, this software is 
#          provided under the Gnu GPL (V.3) license. To use this software in
#          commercial application, please contact the author. 
#
#
# Author: Christian O'Reilly (christian.oreilly@umontreal.ca)
# Date  : June 14, 2012
#
# [1] O'Reilly, C. & Nielsen, T., A novel method for assessing sleep  
#     spindle propagation based on 2D cross-correlation of S-transform spectra
#     (submitted to Sleep).
#
###############################################################################


from scipy import array, arange, sqrt, mean, concatenate, zeros, fft
from scipy.io import loadmat
from scipy.fftpack import fftfreq
from scipy.integrate import trapz

from copy import deepcopy
from spyndle import computeST
import uuid
import bisect
import re
import datetime
import time
import numpy as np

from functools import total_ordering
from lxml import etree  
from abc import ABCMeta, abstractmethod
    

from spyndle import CycleDefinitions, computeDreamCycles    
from spyndle.errorMng import ErrPureVirtualCall
from spyndle import Filter

"""
 Abstract class describing an EEG reader.
"""
class EEGDBReaderBase(object, metaclass=ABCMeta) :
    def __init__(self, recordDuration = 2):
        self.setRecordDuration(recordDuration)         
        self._recordsInfo    = RecordList()
        self._events        = EventList()
        self.__pageSize = None

    @property
    def recordsInfo(self):
        return self._recordsInfo 
        
    #@recordsInfo.setter
    #def recordsInfo(self, recordsInfo):
    #    self._recordsInfo = recordsInfo 

    def getInfoRecords(self, noRecord=None):
        """
         Usge of the recordInfo property allows for a cleaner syntaxe. Be
         careful however self.getInfoRecords(x) == self.recordInfo[x-1]
        """
        if noRecord is None:
            return self._recordsInfo    
        else:
            return self._recordsInfo[noRecord-1]
    
    # Abstract accessor functions   
     
    """
     Implementation of subclasses must return a list of channel labels (names).
    """



    # Should return a DateTime object indicating the time when the recording
    # has begun. 
    @abstractmethod
    def getChannelFreq(self, channel): raise ErrPureVirtualCall
    
    @abstractmethod
    def getFileName(self): raise ErrPureVirtualCall
    

    # Should return a DateTime object indicating the time when the recording
    # has begun. 
    @abstractmethod
    def getRecordingStartTime(self): raise ErrPureVirtualCall

    @abstractmethod
    def getChannelLabels(self): raise ErrPureVirtualCall
        



    ###########################################################################
    # Events
    ##########################################################################


    @property
    def events(self, startTime=None, endTime=None):
        return self.getEvents(startTime, endTime)
    
    @events.setter
    def events(self, events):
        if not isinstance(events, EventList):
            raise TypeError("The argument events must be of EventList type.")
        self._events = events        
        
    # This function should not be called from outside the class anymore. The
    # event property should be used instead.
    def getEvents(self, startTime=None, endTime=None) :
        if startTime is None and endTime is None:
            return self._events
            
        if startTime is None :
            return [e for e in self._events if e.startTime < endTime or 
                                 e.startTime + e.timeLength < endTime]                  
            
        if endTime is None :            
            return [e for e in self._events if e.startTime >= startTime or 
                                    e.startTime + e.timeLength >= startTime]       
                
        return [e for e in self._events if (e.startTime >= startTime and e.startTime < endTime) or 
                     (e.startTime + e.timeLength >= startTime and e.startTime + e.timeLength < endTime)]       
                            
        

            

    def addEvent(self, event):
        self.events.add(event)   

    def removeEvent(self, event):
        self.events.remove(event)
        
    def removeEventType(self, eventType):
        self.events.removeType(eventType)        




        
    @abstractmethod        
    def getNbRecords(self): raise ErrPureVirtualCall
        
    @abstractmethod              
    def readChannel(self, channel, usePickled=False) :  raise ErrPureVirtualCall     
            
    @abstractmethod              
    def getNbSample(self, channel=None) :  raise ErrPureVirtualCall     


    """
     Return an array of the time associated with every sample of the signal
     recorded by the channel.
    """
    @abstractmethod    
    def getChannelTime(self, channel, startTime, timeDuration)  :  raise ErrPureVirtualCall           
        
  
    def setRecordDuration(self, duration):
        self.recordDuration = duration
        
    def getRecordDuration(self):
        return self.recordDuration
        
        
    def getPages(self):
        stageEvents = [e for e in self.events if e.groupName.lower() == "stage"]    
        return [e.startTime for e in stageEvents], \
               [e.dateTime for e in stageEvents]
        
    def getSleepStage_REM(self):
        return [e for e in self.reader.events if e.name == "Sleep stage R"]      
        
    def getSleepStage_NREM(self):
        return [e for e in self.reader.events if e.name == "Sleep stage 1" or 
                                                e.name == "Sleep stage 2" or 
                                                e.name == "Sleep stage 3" or 
                                                e.name == "Sleep stage 4" or 
                                                e.name == "Sleep stage N" or 
                                                e.name == "Sleep stage N1" or 
                                                e.name == "Sleep stage N2" or 
                                                e.name == "Sleep stage N3"]      
        
    def getSleepStage_1(self):
        return [e for e in self.reader.events if e.name == "Sleep stage 1" or 
                                                 e.name == "Sleep stage N1" ]     
    def getSleepStage_2(self):
        return [e for e in self.reader.events if e.name == "Sleep stage 2" or 
                                                 e.name == "Sleep stage N2" ]     
        
    def getSleepStage_3(self):
        return [e for e in self.reader.events if e.name == "Sleep stage 3" or 
                                                 e.name == "Sleep stage N3" ]      
        
    def getSleepStage_4(self):
        return [e for e in self.reader.events if e.name == "Sleep stage 4"]      
        
    def getSleepStage_SWS(self):
        return [e for e in self.reader.events if e.name == "Sleep stage 3" or 
                                                e.name == "Sleep stage 4" or 
                                                e.name == "Sleep stage N3"]         
        
        
    def getEventsByTime(self, startTime, endTime) :
        return [e for e in self.events if (e.startTime >= startTime and e.startTime < endTime) or 
                         (e.startTime + e.timeLength >= startTime and e.startTime + e.timeLength < endTime)]     


  


    def computeRMS(self, event, fmin=11, fmax=16):

        # The filters need the signal to be at least 3*nbTaps
        nbTaps      = 1001.0        
        fs          = self.getChannelFreq(event.channel) 
        duration    = max((nbTaps*3.0+1)/fs, event.duration())   
        data        = self.read([event.channel], event.timeStart(), duration)
        signal      = data[event.channel].signal
        
        # In some cases, the signal will be shorter than nbTaps*3.0/fs+1,
        # for example, when the spindle is too close to the end of the recording
        # epoch. In these time, we have to accept lower number of taps.
        if len(signal) <  nbTaps*3.0+1:
            nbTaps = int((len(signal)-1)/3)
        
        # Defining EEG filters
        bandPassFilter = Filter(fs)
        bandPassFilter.create(low_crit_freq=fmin, high_crit_freq=fmax, 
                              order=nbTaps, btype="bandpass", ftype="FIR", useFiltFilt=True)          
             
        try :
            signal = bandPassFilter.applyFilter(signal)[0:int(event.duration()*fs)]                                
        except ValueError:
            print((len(signal), nbTaps, max(nbTaps*3.0/fs, event.duration()), nbTaps*3.0/fs, event.duration())) 
            raise
        
        event.properties["RMSAmp"] = np.sqrt(np.mean(signal**2))

    
    def computeMeanFreq(self, event, fmin=11, fmax=16):
   
        fs          = self.getChannelFreq(event.channel) 
        data        = self.read([event.channel], event.timeStart(), event.duration())
        signal      = data[event.channel].signal

        if signal.size < 512:
            signal = concatenate((signal, zeros(512 - signal.size)))
        
        
        FFT = abs(fft(signal))
        freqs = fftfreq(signal.size, 1.0/fs)        
        
        FFT   = FFT[(freqs >= fmin)*(freqs <= fmax)]
        freqs = freqs[(freqs >= fmin)*(freqs <= fmax)]
        
        event.properties["meanFreq"] = sum(freqs*FFT)/sum(FFT)




    def computeFreqSlope(self, event, fmin=11, fmax=16):
   
        fs          = self.getChannelFreq(event.channel) 
        data        = self.read([event.channel], event.timeStart(), event.duration())
        signal      = data[event.channel].signal

        nbMinSamples = self.getChannelFreq(event.channel)
        N            = len(signal)
        if len(signal) < nbMinSamples:
            signal = np.concatenate((signal, np.zeros(nbMinSamples - N)))
            
        X, fX = computeST(signal, fs, fmin=fmin-1, fmax=fmax+1)  

        Y = abs(np.transpose(X))
                        
        regx = arange(N)/fs
        regy = []
        for i in range(N):
            try:
                result = trapz(fX*Y[:, i], fX)/trapz(Y[:, i], fX)
                regy.append(result)  
                #if np.isnan(result) or np.isnan(result):
                #    print fX, Y[:, i], len(signal), event.duration()
            except:
                print((fX.shape, Y.shape, fX.shape))
                print((i, self.recordsInfo[-1].startTime, self.recordsInfo[-1].duration, event.timeStart(), event.duration()))
                raise

        assert(not np.any(np.isnan(regy)))
        assert(not np.any(np.isinf(regy)))
        z = np.polyfit(regx, regy, deg=1)     

        event.properties["slopeOrigin"] = z[1]
        event.properties["slope"]       = z[0]







    def getEventIndicator(self, events, channel=None, time=None, checkUniqueness=False):
 
        if time is None :
            if channel is None:
                raise ValueError("The channel must be specified if no time vector is provided.")
            time = self.getChannelTime(channel)  
        eventInd = []
        indMin = 0
        for event in events:
            start  = bisect.bisect_left(time[indMin:], event.timeStart())
            stop   = bisect.bisect_left(time[indMin:], event.timeEnd())
            eventInd.extend(list(range(start+indMin, stop+indMin)))
            indMin += stop
        eventInd = np.array(eventInd)
        
        if checkUniqueness:
            assert(len(eventInd) == len(np.unique(eventInd)))
        
        indicators = np.zeros(len(time), dtype=bool)
        if len(eventInd):
            indicators[eventInd] = True
        return indicators
        




    def getEventNameIndicator(self, eventName, channel, checkUniqueness=False, globalEvent=False):
        """
         Return an array containing the same number of items than the number
         of sample in the channel signal, with ones for samples falling in a 
         window of event eventName and zeros elsewhere.
         
         globalEvent must be false for events that have no specific channel
         such as stage events. However, a channel must nevertheless be passed
         to get the right number of stamples for a given channel.
        """          
        if globalEvent:
            if isinstance(eventName, list):
               events = [e for e in self.events if e.name in eventName]                
            elif isinstance(eventName, str):
               events = [e for e in self.events if e.name == eventName]
            else:
               raise ValueError("The eventName arguement must be a list or a str object.")
        else:     
            if isinstance(eventName, list):
               events = [e for e in self.events if e.name in eventName and e.channel == channel]                
            elif isinstance(eventName, str):
               events = [e for e in self.events if e.name == eventName and e.channel == channel]
            else:
               raise ValueError("The eventName arguement must be a list or a str object.")
         
        return self.getEventIndicator(events, channel=channel, checkUniqueness=checkUniqueness)

        

    @property
    def pageSize(self):
        if self.__pageSize is None:
            stageEvents = [e for e in self.events if e.groupName.lower() == "stage"] 
            self.__pageSize = np.mean([p.duration() for p in stageEvents])
        return self.__pageSize


    def setStagesToEventTypes(self, eventName):
        for event in [e for e in self.events if e.name == eventName] :
            self.setStagesToEvent(event)

    def setStagesToEvent(self, event):

        # Select the stage where the event begin as the sleep stage
        # of the event.
        indMin = self.events.getIndexMinStartTime(event.startTime-self.pageSize)
        indMax = self.events.getIndexMaxStartTime(event.startTime)
        if indMax <= indMin:
            event.properties["stage"] = "No stage"
        else:
            stages = [e for e in self.events[indMin:(indMax+1)] if e.groupName.lower() == "stage"]                  
            if len(stages) == 1 :
                event.properties["stage"] = stages[0].name
            elif len(stages) == 0 :    
                event.properties["stage"] = "No stage"
            else:
                event.properties["stage"] = stages[0].name           
                print(("Waring: " + str(len(stages)) + " staged are including the event starting at time " + str(event.startTime) + ". Keeping only the first one."))



  
    def getIndexMinStartTime(self, timeMin, inclusive=True):
        """
        Return the index of the self.__events list that correspond to the
        first event with a starting time larger (inclusive == False) or
        larger or equal (include == True) than timeMin.
        """   
    def getIndexMaxStartTime(self, timeMax, inclusive=True):
        """
        Return the index of the self.__events list that correspond to the
        last event with a starting time smaller (inclusive == False) or
        smaller or equal (include == True) than timeMin.
        """   




    def setCyclesToEventTypes(self, eventNames):

        if isinstance(eventNames, str):
            eventNames = [eventNames]
        
        self.cycleDef = CycleDefinitions()
        self.cycleDef.setAeschbach()
        cycles = computeDreamCycles([e for e in self.events if e.groupName.lower() == "stage"], self.cycleDef)
        
        nbEvents = len([e for e in self.events if e.name in eventNames])
        gotNbEvents = 0

        #cycles[0].timeStart() > event.startTime():
        try:
            indUp = self.events.getIndexMaxStartTime(cycles[0].timeStart(), inclusive=False)    
            events = [e for e in self.events[:indUp+1] if e.name in eventNames] 
            gotNbEvents += len(events)
            for event in events:
                event.cycle = 1
        except ValueError:
            pass
       
        for i, cycle in enumerate(cycles) :
            indUp  = self.events.getIndexMaxStartTime(cycle.timeEnd(), inclusive=False)    
            indLow = self.events.getIndexMinStartTime(cycle.timeStart(), inclusive=True)    
            #cycle.timeStart() <= event.startTime() and cycle.timeEnd() > event.startTime()
            events = [e for e in self.events[indLow:indUp+1] if e.name in eventNames]                
            gotNbEvents += len(events)
            for event in events:
                event.cycle = i+1

        # cycles[-1].timeEnd() <= event.startTime()
        try:
            indLow = self.events.getIndexMinStartTime(cycles[-1].timeEnd(), inclusive=True)   
            events = [e for e in self.events[indLow:] if e.name in eventNames]           
            gotNbEvents += len(events)
            for event in events:
                    event.cycle = len(cycles)
        except ValueError:
            pass

        assert(gotNbEvents == nbEvents)


    def getDiscontinuityInfo(self, deltaMax = 0.002):
        starts = [self.recordsInfo[0].startTime]
        ends   = []
        noStarts = [0]
        noEnds   = []
        for noRecord, (recInfoLast, recInfo) in enumerate(zip(self.recordsInfo[:-1], self.recordsInfo[1:])):
            finLast = recInfoLast.startTime+recInfoLast.duration
            if recInfo.startTime - finLast > deltaMax:                
                ends.append(finLast)
                starts.append(recInfo.startTime)
                noStarts.append(noRecord+1)
                noEnds.append(noRecord)

        ends.append(self.recordsInfo[-1].startTime+self.recordsInfo[-1].duration)
        noEnds.append(len(self.recordsInfo)-1)
        return list(zip(starts, ends, noStarts, noEnds))


  
"""
  Class implementing a list of events. It must manage its internal list such
  that it is always sorted to allow to find events as function of their timing
  of occurence rapidely. 
"""
class TimeOrderedList(object):
    """
     Inheriting classes should defined setter and getter for a memberList 
     property.
    """    
    
    
    def __init__(self): 
        #self.__events = []
        self.memberList = []


    def add(self, event):
        bisect.insort(self.memberList, event)

    def __getitem__(self, index):
        return self.memberList[index]
        
    def __len__(self):
        return len(self.memberList)

    def __str__(self):
        eventNames = [e.name for e in self.memberList]
        
        retStr = "The list object contains " + str(len(eventNames)) + " memberss:\n"
        for eventName in np.unique(eventNames):
            retStr += eventName + ":" + str(np.sum(np.in1d(eventNames, [eventName]))) + "\n"
        return retStr


    def removeType(self, eventType):
        self.memberList = [e for e in self.memberList if e.name != eventType]

    def remove(self, events):
        if isinstance(events, list):
            for event in events :
                if event in self.memberList:
                    self.memberList.remove(event)
                    
        elif isinstance(events, Event):
            self.memberList.remove(events)
        else:
            raise TypeError
  

    def getIndexMinStartTime(self, timeMin, inclusive=True):
        """
        Return the index of the self.memberList list that correspond to the
        first event with a starting time larger (inclusive == False) or
        larger or equal (include == True) than timeMin.
        """         
         
        def find_ind_gt(a, x):
            # Find the index of the leftmost value greater than x
            i = bisect.bisect_right(a, x)
            if i != len(a):
                return i
            raise ValueError            

        def find_ind_ge(a, x):
            # Find the index of the leftmost item greater than or equal to x
            i = bisect.bisect_left(a, x)
            if i != len(a):
                return i
            raise ValueError
            
            
        dummyMember = deepcopy(self.memberList[0])         
        dummyMember.sortingAttribute = timeMin
        if inclusive:
            return find_ind_ge(self.memberList, dummyMember)
        else:
            return find_ind_gt(self.memberList, dummyMember)            



    def getItemMinStartTime(self, timeMin, inclusive=True):
        try:
            return self.memberList[self.getIndexMinStartTime(timeMin, inclusive)]
        except ValueError:
            return None

              
  

    def getIndexMaxStartTime(self, timeMax, inclusive=True):
        """
        Return the index of the self.__events list that correspond to the
        last event with a starting time smaller (inclusive == False) or
        smaller or equal (include == True) than timeMin.
        """                

        def find_ind_lt(a, x):
            # Find the index of the rightmost value less than x
            i = bisect.bisect_left(a, x)
            if i:
                return i-1
            raise ValueError



        def find_ind_le(a, x):
            # Find the index of the rightmost value less than or equal to x
            i = bisect.bisect_right(a, x)
            if i:
                return i-1
            raise ValueError
        dummyMember = deepcopy(self.memberList[0])     
        dummyMember.sortingAttribute = timeMax
        if inclusive:
            return find_ind_le(self.memberList, dummyMember)
        else:
            return find_ind_lt(self.memberList, dummyMember)     
            
              

    def getItemMaxStartTime(self, timeMin, inclusive=True):
        try:
            return self.memberList[self.getIndexMaxStartTime(timeMin, inclusive)]
        except ValueError:
            return None


  
  

"""
  Class implementing a list of events. It must manage its internal list such
  that it is always sorted to allow to find events as function of their timing
  of occurence rapidely. 
"""
class EventList(TimeOrderedList):
    
    def __init__(self):
        self.__events = []
        
    @property
    def memberList(self):
        return self.__events

    @memberList.setter
    def memberList(self, eventList):
        self.__events = eventList

  
  
class RecordList(TimeOrderedList):
    
    def __init__(self):
        self.__records = []
    
    @property
    def memberList(self):
        return self.__records

    @memberList.setter
    def memberList(self, recordList):
        self.__records = recordList

      
  
  
  
  
# TODO: Manage discontinuous signals.

class EEGRecord:  
    
    def __init__(self, samplingRates, recordStart, recordedSignals, recordingStartTime) :
        
        if not isinstance(recordStart,          datetime.datetime) : raise TypeError
        if not isinstance(recordingStartTime, datetime.datetime) : raise TypeError
        if not isinstance(samplingRates,      dict) : raise TypeError
        if not isinstance(recordedSignals,    dict) : raise TypeError
            

        self.recordStartTime      = recordStart
        self.recordingStartTime = recordingStartTime
        
        # Dictionaries indexed by the channel label.
        self.recordedSignals    = recordedSignals
        self.samplingRates      = samplingRates


        
    """
     Return the datetime object corresponding to the starting of the record.
    """
    def getStartDateTime(self):
        return self.recordingStartTime
  
    """
     Return the time corresponding to the starting of the record in seconds
     since the begining of the recording.
    """
    def getStartTime(self):
        return (self.recordStartTime  - self.recordingStartTime).total_seconds()   
  
  

    def getDuration(self):
        for channel in self.recordedSignals:
            return float(len(self.recordedSignals[channel]))/self.samplingRates[channel]
        
  
    def getSignal(self, channel):
        return self.recordedSignals[channel] 
        
    def getSignalTime(self, channel):
        return arange(len(self.recordedSignals[channel]))/self.samplingRates[channel] + self.getStartTime()
  
  
  
  
@total_ordering
class SortedMember(metaclass=ABCMeta):
    def __eq__(self, other):
        return self.sortingAttribute == other.sortingAttribute
    def __lt__(self, other):
        return self.sortingAttribute < other.sortingAttribute       
        
    @property
    def sortingAttribute(self): raise ErrPureVirtualCall
        
    @sortingAttribute.setter
    def sortingAttribute(self, val): raise ErrPureVirtualCall
        

"""
 Gives some really basic information on a record.
"""
class EEGRecordInfo(SortedMember):
       
    @property
    def sortingAttribute(self):
        return self.startTime  
        
    @sortingAttribute.setter
    def sortingAttribute(self, val): 
        self.startTime = val
        
        
        
    def __init__(self, startSample=None, endSample=None, 
                       startTime=None,   duration=None, isComplete=True):
        
        if isinstance(startSample, (int, np.int32)) or startSample is None:
            self.startSample = startSample
        else:
            raise TypeError("Type should be and int or None. Receive type is " + str(type(startSample)))             
            
        if isinstance(endSample, (int, np.int32)) or endSample is None:
            self.endSample   = endSample
        else:
            raise TypeError("Type should be and int or None. Receive type is " + str(type(endSample)))            
            
        
        if isinstance(startTime, float) or startTime is None:
            self.startTime = startTime
        else:
            raise TypeError            
            
        if isinstance(duration, float) or duration is None:
            self.duration   = duration
        else:
            raise TypeError            
            
        if isinstance(isComplete, bool) or isComplete is None:
            self.isComplete  = isComplete
        else:
            raise TypeError            
            
        
        
        
    def getNbSamples(self):
        return self.endSample - self.startSample+1
  
  
  
  
  
  
  
  
# TODO: Manage discontinuous signals.
class RecordedChannel:
    def __init__(self): #, ISignalFile):
        self.signal         = array([])
        self.samplingRate   = 0.0
        self.type           = ""
        self.startTime      = ""
        
    # This mimic the behavior of Python's dictionary (needed to use the 
    #savemat function).
    def items(self):
        return [("signal", self.signal), ("samplingRate", self.samplingRate), 
                ("type", self.type), ("startTime", self.startTime)]
        
        
    def readPickledChannel(self, fileName):
        data = loadmat(fileName)        
        self.signal         = data["signal"].reshape(data["signal"].size)   
        self.samplingRate   = data["samplingRate"][0][0]    
        self.type           = data["type"][0][0]  
        self.startTime      = time.strptime(data["startTime"][0], "%a, %d %b %Y %H:%M:%S +0000")

    def __str__(self):
        return ("signal:" + str(self.signal) + ", samplingRate:" + 
                 str(self.samplingRate) + ", type:" + str(self.type) +
                 ", startTime:" + str(self.startTime) )
                  
        
        
        
        
        

"""
class EventGroup:
    
    def __init__():
                     
        self.groupName   = groupName
        self.color       = "black"
"""



# Return true if c is a valid character for an XML string
def valid_xml_char_ordinal(c):
    if isinstance(c, str):
        codepoint = ord(c)
    else:
        raise TypeError("Must be a int or a str. Found: " + str(type(c)))
                    
        
    # conditions ordered by presumed frequency
    return (
        0x20 <= codepoint <= 0xD7FF or
        codepoint in (0x9, 0xA, 0xD) or
        0xE000 <= codepoint <= 0xFFFD or
        0x10000 <= codepoint <= 0x10FFFF
        )    

def valid_xml_string(input_string):
    if not isinstance(input_string, str):
        raise TypeError("Must be a str. Found: " + str(type(input_string)))
    
    for c in input_string:
        if not valid_xml_char_ordinal(c):
            return False
    return True


class Event(SortedMember, metaclass=ABCMeta):
    """      
    Abstract class representing an EEG event.    
    
    Note:   
    In order to facilitate the cross-operation of readers for multiple data
    format, some event names and group names are standerdized. It is to the
    specific reader to ensure that data format using differents names transale
    them in correct standerdized designation if these events are to be 
    recognize correctly by other classes, such as detectors. 
        
    Standard group names : 
        Stage : For all events used in sleep stage scoring.
    
    Standard names: (these are compatible with http://www.edfplus.info/specs/edftexts.html#annotations)
        Sleep stage 1
        Sleep stage 2
        Sleep stage 3
        Sleep stage 4
        Sleep stage R
        Sleep stage W
        Sleep stage ?
        Sleep stage N 	
        Sleep stage N1 	
        Sleep stage N2 	
        Sleep stage N3
    """    
    

    
    
    @property
    def sortingAttribute(self): 
        return self.startTime    
        
    @sortingAttribute.setter
    def sortingAttribute(self, val): 
        self.startTime = val
        
    
    
    @property
    def channel(self):
        return self._channel
    
    @channel.setter
    def channel(self, value):
        if valid_xml_string(value):
            self._channel = value
        else:
            raise ValueError("Must be XML compatible: Unicode or ASCII, no NULL bytes or control characters.")
    
    
    @property
    def groupName(self):
        return self._groupName
    
    @groupName.setter
    def groupName(self, value):
        if valid_xml_string(value):
            self._groupName = value
        else:
            raise ValueError("Must be XML compatible: Unicode or ASCII, no NULL bytes or control characters.")
    
    @property
    def name(self):
        return self._name
    
    @name.setter
    def name(self, value):
        if valid_xml_string(value):
            self._name = value
        else:
            raise ValueError("Must be XML compatible: Unicode or ASCII, no NULL bytes or control characters.")
    
    def __init__(self, name = "", groupName = "", channel = "", startTime = -1.0,
                 timeLength = -1.0, dateTime = None, properties = {}):
                     
        self.ID          = str(uuid.uuid1())
        self._groupName   = groupName
        self._channel     = channel
        self._name        = name
        self.startTime   = startTime   # In seconds, since the begining of the recording of the EEG file.
        self.timeLength  = timeLength  # Duration in seconds      
        self.dateTime    = dateTime    # datetime  object giving the begining time of the event.

        self.properties  = properties
    
        
        
    def timeEnd(self):
        return self.startTime + self.timeLength        
    def timeStart(self):
        return self.startTime   
    def duration(self):
        return self.timeLength  

            
    def __str__(self):
        return("groups:" +str(self.groupName) + " channel:" + str(self.channel)
                + " name:" + str(self.name) + " startTime:" + str(self.startTime) 
                + " duration:"  
                + str(self.timeLength) + " properties:" + str(self.properties))


    #def __eq__(self, other): 
    #    return self.__dict__ == other.__dict__


    def getXml(self):
        # create XML 
        try:
            root = etree.Element('Event', name=self.name, groupName=self.groupName, channel=self.channel)
            for propKey in self.properties:
                #propertyElem = etree.Element('Property')
                
                # XML properties cannot contain space characters. Substituting them by "_".                
                #propertyElem.set(propKey.replace(' ', '_'), str(self.properties[propKey]))    
                
                #root.append(propertyElem)
                
                # We replace all non alphanumeric characters by "_" to make sure
                # the property name is valid for XML representation.
                root.set(str(re.sub('[^0-9a-zA-Z]+', '_', propKey)), str(self.properties[propKey])) 
                
        except ValueError :
            print(self.name, self.groupName, self.channel)
            print(self.properties)
            print(self.properties.keys())
            #print(("propKey:", str(propKey),propKey)) 
            #print(("property:", str(self.properties[propKey]), self.properties[propKey]))
            raise
            
            
        return etree.tostring(root, encoding=str) #, pretty_print=True)        

    
            
    def toEDFStr(self):
        
        if self.groupName.lower() == "stage":
            eventStr = self.name
        else:
            eventStr = self.getXml()
        
        return "+" + str(self.startTime) + "\x15" + str(self.timeLength) + "\x14" + eventStr + "\x14\0"    
        