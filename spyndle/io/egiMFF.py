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
###############################################################################

# -*- coding: utf-8 -*-

'''
 We consider the record concept of EDF as equivalent to the record concept in
 Harmony and in sleep stage scoring in general. Thus, to use 30 second records, 
 the EDF file must be reformatted if the record duration is different than
 30 seconds.
'''


from .EEGDatabaseReader import EEGDBReaderBase, EEGRecord

import uuid
import os, io
import re, datetime, logging
import numpy as np
#import scipy.weave as weave
import tempfile
import warnings
import struct
#from scipy.weave.build_tools import CompileError
import pickle as pickle


from scipy import array, arange, concatenate, where
from copy import deepcopy, copy
from lxml import etree
from tempfile import gettempdir 
from time import sleep
from warnings import warn
from shutil import copyfile, move

from .EEGDatabaseReader import Event, RecordedChannel, EEGRecordInfo

import numpy

log = logging.getLogger(__name__)

from wonambi.ioeeg.egimff import EgiMff

   
class EgiMffReader(EEGDBReaderBase) :

    def __init__(self, fileName):        
        self.fileName = fileName
        self.wonambiReader = EgiMff(fileName)
        self.header = self.wonambiReader.return_hdr()
        self.preloaded = {}
 
    def getFileName(self): 
        return self.fileName       

    def getRecordingStartTime(self): 
        return self.header[1]      

    def getNbRecords(self):    
        return np.nan
       
    def getChannelLabels(self):
        return self.header[3]       
     

    def getChannelFreq(self, channel=None): 
        return float(self.header[2])   


    def getChannelTime(self, channel=None, startTime=0, timeDuration=np.inf) :        
        endTime = startTime + timeDuration
        time = np.arange(self.getNbSample())/self.getChannelFreq()       
        #return time[np.where((time >= startTime)*(time < endTime))]

        begsam = int(startTime*self.getChannelFreq())
        if endTime is None:
            endsam = int(self.getNbSample())        
        else:
            endsam = int(endTime*self.getChannelFreq())

        return time[begsam:endsam]



    """
     Read the complete signal recorded by a given channel.
    """
    def readChannel(self, signalName, startTime=0, endTime=None):
        #print("readChannel:", signalName, startTime, endTime)

        if not isinstance(signalName, str) :
            raise TypeError("The signalName argument must be a string. Received: " + str(type(signalName)))
        if not signalName in self.getChannelLabels():
            raise ValueError("Unknown channel.")

        data = RecordedChannel()         
        data.samplingRate   = self.getChannelFreq(signalName)
        data.type           = "EEG"
        data.startTime      = startTime

        begsam = int(startTime*data.samplingRate)
        if endTime is None:
            endsam = int(self.getNbSample())        
        else:
            endsam = int(endTime*data.samplingRate)

        if signalName in self.preloaded:
            data.signal = deepcopy(self.preloaded[signalName].signal[begsam:endsam])
        else:
            data.signal = self.wonambiReader.return_dat([self.getChannelLabels().index(signalName)], 
                                                         begsam=begsam, 
                                                         endsam=endsam)[0]
        return data    



    
    def renameChannel(self, oldName, newName):
        self.header[3][self.header[3].index(oldName)] = newName
    

    def preload(self, signalNames):
        self.preloaded = {}  
        for chan in signalNames:
             self.preloaded[chan] = self.readChannel(chan) 

    
    
    def read(self, signalNames, startTime, timeDuration, debug=False):
        returnData = {}  
        for chan in signalNames:
             returnData[chan] = self.readChannel(chan, startTime=startTime, endTime=startTime+timeDuration) 
        return returnData


    def getNbSample(self, channel=None):
        return self.header[4]

      


