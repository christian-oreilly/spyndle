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
    


class ErrPureVirtualCall(Exception): 
    
#    def __init__(self, err):
 #       self.value = value
        
    def __str__(self):
        return "Trying to access a pure virtual function."



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
        
        
        
        
        
    def __init__(self, pageDuration = 30): # en secondes
        self.setPageDuration(pageDuration) 
  
    def setPageDuration(self, duration):
        self.pageDuration = duration
        
    def getPageDuration(self):
        return self.pageDuration
        
        
  
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

            
        