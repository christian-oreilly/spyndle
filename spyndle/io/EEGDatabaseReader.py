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
from abc import ABCMeta, abstractmethod
    

class EEGDBReaderBase :
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def getChannelLabels(self):
        return
        
    @abstractmethod
    def getEvents(self):
        return
    
    def __init__(self):
        self.pageDuration = 30 # en secondes

        
    def setPageDuration(self, duration):
        self.pageDuration = duration
    
    

    
  
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
        self.timeLength  = ""              
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



            
        