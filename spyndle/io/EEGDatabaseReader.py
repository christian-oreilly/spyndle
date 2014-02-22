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

import uuid
import bisect
import re
import datetime
import time
import numpy as np

from functools import total_ordering
from lxml import etree  
from abc import ABCMeta, abstractmethod
    
from spyndle.errorMng import ErrPureVirtualCall
from spyndle import Filter

"""
 Abstract class describing an EEG reader.
"""
class EEGDBReaderBase(object) :
    __metaclass__ = ABCMeta

    # List of EEGPageInfo of the recording.
    def __init__(self, pageDuration = 30):
        self.setPageDuration(pageDuration)         
        self.__infoPages    = []
        self._events        = EventList()

    def getInfoPages(self, noPage=None):
        if noPage is None:
            return self.__infoPages    
        else:
            return self.__infoPages[noPage-1]
    
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
        self._events = events        
        
    # This function should not be called from outside the class anymore. The
    # event property should be used instead.
    def getEvents(self, startTime=None, endTime=None) :
        if startTime is None and endTime is None:
            return self._events
            
        if startTime is None :
            return filter(lambda e: e.startTime < endTime or 
                                 e.startTime + e.timeLength < endTime , self._events)                  
            
        if endTime is None :            
            return filter(lambda e: e.startTime >= startTime or 
                                    e.startTime + e.timeLength >= startTime, self._events)       
                
        return filter(lambda e: (e.startTime >= startTime and e.startTime < endTime) or 
                     (e.startTime + e.timeLength >= startTime and e.startTime + e.timeLength < endTime) , self._events)       
                            
        

            

    def addEvent(self, event):
        self.events.add(event)   

    def removeEvent(self, event):
        self.events.remove(event)
        
    def removeEventType(self, eventType):
        self.events.removeType(eventType)        




        
    @abstractmethod        
    def getNbPages(self): raise ErrPureVirtualCall
        
    @abstractmethod              
    def readChannel(self, channel, usePickled=False) :  raise ErrPureVirtualCall     
            
    @abstractmethod              
    def getNbSample(self, channel=None) :  raise ErrPureVirtualCall     


    """
     Return an array of the time associated with every sample of the signal
     recorded by the channel.
    """
    @abstractmethod    
    def getChannelTime(self, channel)  :  raise ErrPureVirtualCall           
        
  
    def setPageDuration(self, duration):
        self.pageDuration = duration
        
    def getPageDuration(self):
        return self.pageDuration
        
        
    def getPages(self):
        stageEvents = filter(lambda e: e.groupName.lower() == "stage", self.events)    
        return map(lambda e: e.startTime  , stageEvents), \
               map(lambda e: e.dateTime   , stageEvents)
        
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
        return filter(lambda e: (e.startTime >= startTime and e.startTime < endTime) or 
                         (e.startTime + e.timeLength >= startTime and e.startTime + e.timeLength < endTime) , self.events)     


  


    def computeRMS(self, event, fmin=11, fmax=16):

        # The filters need the signal to be at least 3*nbTaps
        nbTaps      = 1001.0        
        fs          = self.getChannelFreq(event.channel) 
        duration    = max(nbTaps*3.0/fs+1, event.duration())   
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
            print len(signal), nbTaps, max(nbTaps*3.0/fs, event.duration()), nbTaps*3.0/fs, event.duration() 
            raise
        
        event.properties["RMSamp"] = sqrt(mean(signal**2))

    
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



    def getEventIndicator(self, eventName, channel, checkUniqueness=False, globalEvent=False):
        """
         Return an array containing the same number of items than the number
         of sample in the channel signal, with ones for samples falling in a 
         window of event eventName and zeros elsewhere.
         
         globalEvent must be false for events that have no specific channel
         such as stage events. However, a channel must nevertheless be passed
         to get the right number of stamples for a given channel.
        """
        t = self.getChannelTime(channel)  
          
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
         

        """
        The following code is equivalent to :
            eventInd = np.concatenate([np.where((t >= e.timeStart())*(t < e.timeEnd()))[0] for e in events]) 
        but it runs much faster by benifiting from the fact that t is an ordered
        vector.
        """
        eventInd = []
        indMin = 0
        for event in events:
            start  = bisect.bisect_left(t[indMin:], event.timeStart())
            stop   = bisect.bisect_left(t[indMin:], event.timeEnd())
            eventInd.extend(range(start+indMin, stop+indMin))
            indMin += stop
        eventInd = np.array(eventInd)


        
        if checkUniqueness:
            assert(len(eventInd) == len(np.unique(eventInd)))
        
        eventIndicator = zeros(self.getNbSample(channel))
        if len(eventInd):
            eventIndicator[eventInd] = 1            
        return eventIndicator
        




    def setStagesToEventTypes(self, eventName):
        for event in filter(lambda e: e.name == eventName, self.events) :
            self.setStagesToEvent(event)

    def setStagesToEvent(self, event):

        # Select the stage where the event begin as the sleep stage
        # of the event.

        indMin = self.events.getIndexMinStartTime(event.startTime-self.getPageDuration())
        indMax = self.events.getIndexMaxStartTime(event.startTime)
        if indMax <= indMin:
            event.properties["stage"] = "No stage"
        else:
            stages = filter(lambda e: e.groupName.lower() == "stage", self.events[indMin:(indMax+1)])                  
            if len(stages) == 1 :
                event.properties["stage"] = stages[0].name
            elif len(stages) == 0 :    
                event.properties["stage"] = "No stage"
            else:
                event.properties["stage"] = stages[0].name           
                print "Waring: " + str(len(stages)) + " staged are including the event starting at time " + str(event.startTime) + ". Keeping only the first one."



"""
  Class implementing a list of events. It must manage its internal list such
  that it is always sorted to allow to find events as function of their timing
  of occurence rapidely. 
"""
class EventList:
    def __init__(self): 
        self.__events = []


    def add(self, event):
        bisect.insort(self.__events, event)

    def __getitem__(self, index):
        return self.__events[index]
        
    def __len__(self):
        return len(self.__events)

    def __str__(self):
        eventNames = [e.name for e in self.__events]
        
        retStr = "The EventList object contains " + str(len(eventNames)) + " events:\n"
        for eventName in np.unique(eventNames):
            retStr += eventName + ":" + str(np.sum(np.in1d(eventNames, [eventName]))) + "\n"
        return retStr


    def removeType(self, eventType):
        self.__events = [e for e in self.__events if e.name != eventType]

    def remove(self, events):
        if isinstance(events, list):
            for event in events :
                if event in self.__events:
                    self.__events.remove(event)
                    
        elif isinstance(events, Event):
            self.__events.remove(events)
        else:
            raise TypeError
  

    def getIndexMinStartTime(self, timeMin, inclusive=True):
        """
        Return the index of the self.__events list that correspond to the
        first event with a starting time larger (inclusive == False) or
        larger or equal (include == True) than timeMin.
        """         
         
        def find_ind_gt(a, x):
            # Find the index of the leftmost value greater than x
            return bisect.bisect_right(a, x)
            
        def find_ind_ge(a, x):
            # Find the index of the leftmost item greater than or equal to x
            return bisect.bisect_left(a, x)

        dummyEvent = Event()
        dummyEvent.startTime = timeMin
        if inclusive:
            return find_ind_ge(self.__events, dummyEvent)
        else:
            return find_ind_gt(self.__events, dummyEvent)            
              
  

    def getIndexMaxStartTime(self, timeMax, inclusive=True):
        """
        Return the index of the self.__events list that correspond to the
        last event with a starting time smaller (inclusive == False) or
        smaller or equal (include == True) than timeMin.
        """                

        def find_ind_lt(a, x):
            # Find the index of the rightmost value less than x
            return bisect.bisect_left(a, x)-1


        def find_ind_le(a, x):
            # Find the index of the rightmost value less than or equal to x
            return bisect.bisect_right(a, x) -1


        dummyEvent = Event()
        dummyEvent.startTime = timeMax
        if inclusive:
            return find_ind_le(self.__events, dummyEvent)
        else:
            return find_ind_lt(self.__events, dummyEvent)     
            
            
  
# TODO: Manage discontinuous signals.
class EEGPage:
    
    def __init__(self, samplingRates, pageStart, recordedSignals, recordingStartTime) :
        
        if not isinstance(pageStart,          datetime.datetime) : raise TypeError
        if not isinstance(recordingStartTime, datetime.datetime) : raise TypeError
        if not isinstance(samplingRates,      dict) : raise TypeError
        if not isinstance(recordedSignals,    dict) : raise TypeError
            

        self.pageStartTime      = pageStart
        self.recordingStartTime = recordingStartTime
        
        # Dictionaries indexed by the channel label.
        self.recordedSignals    = recordedSignals
        self.samplingRates      = samplingRates


        
    """
     Return the datetime object corresponding to the starting of the page.
    """
    def getStartDateTime(self):
        return self.recordingStartTime
  
    """
     Return the time corresponding to the starting of the page in seconds
     since the begining of the recording.
    """
    def getStartTime(self):
        return (self.pageStartTime  - self.recordingStartTime).total_seconds()   
  
  

    def getDuration(self):
        for channel in self.recordedSignals:
            return float(len(self.recordedSignals[channel]))/self.samplingRates[channel]
        
  
    def getSignal(self, channel):
        return self.recordedSignals[channel] 
        
    def getSignalTime(self, channel):
        return arange(len(self.recordedSignals[channel]))/self.samplingRates[channel] + self.getStartTime()
  
  
  
  

"""
 Gives some really basic information on a page.
"""
class EEGPageInfo:
    def __init__(self, startSample, endSample, isComplete=True):
        
        if isinstance(startSample, int):
            self.startSample = startSample
        else:
            raise TypeError            
            
        if isinstance(endSample, int):
            self.endSample   = endSample
        else:
            raise TypeError            
            
        if isinstance(isComplete, bool):
            self.isComplete  = isComplete
        else:
            raise TypeError            
            
            


        
    def getNbSamples(self):
        return self.endSample - self.startSample
  
  
  
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


@total_ordering
class Event:
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
    
    
    __metaclass__ = ABCMeta
    
    def __init__(self, name = "", groupName = "", channel = "", startTime = -1.0,
                 timeLength = -1.0, dateTime = None, properties = {}):
                     
        self.ID          = str(uuid.uuid1())
        self.groupName   = groupName
        self.channel     = channel
        self.name        = name
        self.startTime   = startTime   # In seconds, since the begining of the recording of the EEG file.
        self.timeLength  = timeLength  # Duration in seconds      
        self.dateTime    = dateTime    # datetime  object giving the begining time of the event.

        self.properties  = properties
        
        
    def __eq__(self, other):
        return self.startTime == other.startTime
    def __lt__(self, other):
        return self.startTime < other.startTime       
        
        
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
                root.set(unicode(re.sub('[^0-9a-zA-Z]+', '_', propKey)), unicode(self.properties[propKey])) 
                
        except ValueError :
            print self.properties
            print self.properties.keys()
            print "propKey:", str(propKey),propKey 
            print "property:", str(self.properties[propKey]), self.properties[propKey]
            raise
            
            
        return etree.tostring(root, encoding=unicode) #, pretty_print=True)        

    
            
    def toEDFStr(self):
        
        if self.groupName.lower() == "stage":
            eventStr = self.name
        else:
            eventStr = self.getXml()
        
        return "+" + str(self.startTime) + "\x15" + str(self.timeLength) + "\x14" + eventStr + "\x14\0"    
        