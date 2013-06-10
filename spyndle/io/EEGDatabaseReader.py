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


from scipy import array, arange
from scipy.io import loadmat

import re
import datetime
import time
from lxml import etree  
from abc import ABCMeta, abstractmethod
    
from spyndle.errorMng import ErrPureVirtualCall

"""
 Abstract class describing an EEG reader.
"""
class EEGDBReaderBase(object) :
    __metaclass__ = ABCMeta

    # List of EEGPageInfo of the recording.
    def __init__(self, pageDuration = 30):
        self.setPageDuration(pageDuration)         
        self.__infoPages   = []

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
    
    

    # Should return a DateTime object indicating the time when the recording
    # has begun. 
    @abstractmethod
    def getRecordingStartTime(self): raise ErrPureVirtualCall

    @abstractmethod
    def getChannelLabels(self): raise ErrPureVirtualCall
        
    @abstractmethod
    def getEvents(self): raise ErrPureVirtualCall
        
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
    In order to facilitate the cross-operation of readers for multiple data
    format, some event names and groupe names are standerdized. It is to the
    specific reader to ensure that data format using differents names transale
    them in correct standerdized designation if these events are to be 
    recognize correctly by other classes, such as detectors. 
    
    Standard groupe names : 
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
class Event:
    def __init__(self, name = "", groupeName = "", channel = "", startTime = -1.0,
                 timeLength = -1.0, dateTime = None, properties = {}):
                     
        self.groupeName  = groupeName
        self.channel     = channel
        self.name        = name
        self.startTime   = startTime   # In seconds, since the begining of the recording of the EEG file.
        self.timeLength  = timeLength  # Duration in seconds      
        self.dateTime    = dateTime    # datetime  object giving the begining time of the event.

        self.properties  = properties

    
    def sampleEnd(self):
        return self.startSample + self.sampleLength        
    def sampleStart(self):
        return self.startSample    
        
    def timeEnd(self):
        return self.startTime + self.timeLength        
    def timeStart(self):
        return self.startTime

            
    def __str__(self):
        return(str(self.groupeName) + " " + str(self.channel)
                + " " + str(self.name) + " " + str(self.startTime) + " " + str(self.timeLength))

    def getXml(self):
        # create XML 
        try:
            root = etree.Element('Event', name=self.name, groupeName=self.groupeName, channel=self.channel)
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
        
        if self.groupeName.lower() == "stage":
            eventStr = self.name
        else:
            eventStr = self.getXml()
        
        return "+" + str(self.startTime) + "\x15" + str(self.timeLength) + "\x14" + eventStr + "\x14\0"    
        