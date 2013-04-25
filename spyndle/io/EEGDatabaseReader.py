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


from scipy import array

from scipy.io import loadmat

import time
from lxml import etree  
from abc import ABCMeta, abstractmethod
    
from spyndle.errorMng import ErrPureVirtualCall

"""
 Abstract class describing an EEG reader.
"""
class EEGDBReaderBase :
    __metaclass__ = ABCMeta
    
    # Abstract accessor functions   
    
    """
     Implementation of subclasses must return a list of channel labels (names).
    """
    @abstractmethod
    def getChannelLabels(self): raise ErrPureVirtualCall
        
    @abstractmethod
    def getEvents(self): raise ErrPureVirtualCall
        
    @abstractmethod        
    def getNbPages(self): raise ErrPureVirtualCall
        
    @abstractmethod              
    def readChannel(self, channel, usePickled=False) :  raise ErrPureVirtualCall     
            
    @abstractmethod              
    def getNbSample(self) :  raise ErrPureVirtualCall     

        
        
        
    def __init__(self, pageDuration = 30): # en secondes
        self.setPageDuration(pageDuration) 
  
    def setPageDuration(self, duration):
        self.pageDuration = duration
        
    def getPageDuration(self):
        return self.pageDuration
        
        
  
  
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
    def __init__(self): #, ISignalFile):
        self.no          = ""
        self.groupeName  = ""
        self.channel     = ""
        self.name        = ""
        self.startTime   = ""
        self.dateTime    = ""      
        self.timeLength  = ""        # Duration in seconds      
        self.startSample = ""
        self.sampleLength= ""
        self.color       = ""

        self.properties  = {}

    
    def sampleEnd(self):
        return self.startSample + self.sampleLength        
    def sampleStart(self):
        return self.startSample    
        
    def timeEnd(self):
        return self.startTime + self.timeLength        
    def timeStart(self):
        return self.startTime

            
    def __str__(self):
        return( str(self.no) + " " + str(self.groupeName) + " " + str(self.channel)
                + " " + str(self.name) + " " + str(self.startTime) + " " + str(self.timeLength))

    def getXml(self):
        # create XML 
        root = etree.Element('Event', name=self.name, groupeName=self.groupeName)
        for propKey in self.properties:
            propertyElem = etree.Element('Property')
            propertyElem.set(propKey, self.properties[propKey])            
            root.append(propertyElem)
            
        # pretty string
        return etree.tostring(root) #, pretty_print=True)        

            
        