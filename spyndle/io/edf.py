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

Reader for EDF+ files.
TODO:
 - add support for log-transformed channels:
   http://www.edfplus.info/specs/edffloat.html and test with 
   data generated with 
   http://www.edfplus.info/downloads/software/NeuroLoopGain.zip.
 - check annotations with Schalk's Physiobank data.

ORIGINAL CODE from Boris Reuderink (Copyright (c) 2012).
Modification by Christian O'Reilly (Copyright (c) 2013)
'''


'''
 We consider the record concept of EDF as equivalent to the page concept in
 Harmony and in sleep stage scoring in general. Thus, to use 30 second pages, 
 the EDF file must be reformated if the record duration is different than
 30 seconds.
'''


from EEGDatabaseReader import EEGDBReaderBase, EEGPage


import os
import re, datetime, logging
import numpy as np

from scipy import array, arange, concatenate

from lxml import etree
from spyndle.io import Event, RecordedChannel



EVENT_CHANNEL = 'EDF Annotations'
log = logging.getLogger(__name__)


"""
 Class representing the Time-stamped Annotations Lists (TALs) used in EDF+ 
 annotation scheme.
"""
#class Tal:

def tal(tal_str):
  '''Return a list with (onset, duration, annotation) tuples for an EDF+ TAL
  stream.
  '''
  exp = '(?P<onset>[+\-]\d+(?:\.\d*)?)' + \
    '(?:\x15(?P<duration>\d+(?:\.\d*)))?' + \
    '(\x14(?P<annotation>[^\x00]*))?' + \
    '(?:\x14\x00)'

  def annotation_to_list(annotation):
    return unicode(annotation, 'utf-8').split('\x14') if annotation else []

  def parse(dic):
    return (
      float(dic['onset']), 
      float(dic['duration']) if dic['duration'] else 0.,
      annotation_to_list(dic['annotation']))

  return [parse(m.groupdict()) for m in re.finditer(exp, tal_str)]


  
class EDFEvent(Event):
    def __init__(self, talEvent): #, ISignalFile):
        Event.__init__(self)
        self.startTime   = talEvent[0]
        self.timeLength  = talEvent[1]        # Duration in seconds  
        #self.dateTime    = ""          


        self.properties  = {}
        # XML event
        if talEvent[2][:6] == "<Event" :
            root = etree.fromstring(talEvent[2])
            for name, value in sorted(root.items()):
                if name == "name":
                    self.name = value
                elif name == "groupeName":
                    self.groupeName = value
                elif name == "channel":
                    self.channel = value
                else :
                    self.properties[name] = value                    
        else:
            self.name = talEvent[2]
            if ( self.name == "Sleep stage 1" or self.name == "Sleep stage 2" or
                 self.name == "Sleep stage 3" or self.name == "Sleep stage 4" or
                 self.name == "Sleep stage R" or self.name == "Sleep stage W" or
                 self.name == "Sleep stage ?" ) :
               
               
               self.groupeName = "Stage"



    def getXml(self):
        # create XML 
        root = etree.Element('Event', name=self.name, groupeName=self.groupeName)
        for propKey in self.properties:
            propertyElem = etree.Element('Property')
            propertyElem.set(propKey, self.properties[propKey])            
            root.append(propertyElem)
            
        # pretty string
        return etree.tostring(root) #, pretty_print=True)        

            
        








class EDFHeader :
    # BDF is a 24-bit version of the 16-bit EDF format, so we read it with the
    # the same re
    def __init__(self, f, fileName):

        assert f.tell() == 0  # check file position

        tampon = f.read(8)
        if tampon == '0       ':
             self.fileType = "EDF"
             self.nbBytes  = 2     # 16 bits
        elif  tampon == '\xFFBIOSEMI':
             self.fileType = "BDF"
             self.nbBytes  = 3     # 24 bits
        else:
             print tampon
             raise IOError

        
        # recording info
        self.subjectID = f.read(80).strip()
        self.recordingIR = f.read(80).strip()
        
        
        # parse timestamp
        (day, month, year) = [int(x) for x in re.findall('(\d+)', f.read(8))]
        (hour, minute, sec)= [int(x) for x in re.findall('(\d+)', f.read(8))]
        
        # Guessing that the data are less than 100 years old....
        year = year + (2000 if year + 2000 < datetime.date.today().year else 1900) 
        self.startDateTime = datetime.datetime(year, month, day, hour, minute, sec)
        
        
        # misc
        self.headerNbBytes      = int(f.read(8))
        self.subtype            = f.read(44)[:5]
        self.contiguous         = self.subtype != 'EDF+D'
        self.nbRecords          = int(f.read(8))
        self.recordDuration      = float(f.read(8))  # in seconds
        self.nbChannels         = int(f.read(4))
        
        # read channel info
        self.channelLabels      = [f.read(16).strip() for n in range(self.nbChannels)]

        self.transducerType     = {}
        self.units              = {}
        self.physicalMin        = {}
        self.physicalMax        = {}          
        self.digitalMin         = {}
        self.digitalMax         = {}
        self.prefiltering       = {}          
        self.nbSamplesPerRecord = {} 
        
        for channel in self.channelLabels : 
            self.transducerType[channel]        = f.read(80).strip()
        for channel in self.channelLabels : 
            self.units[channel]                 = f.read(8).strip()
        for channel in self.channelLabels : 
            self.physicalMin[channel]           = float(f.read(8))
        for channel in self.channelLabels : 
            self.physicalMax[channel]          = float(f.read(8))            
        for channel in self.channelLabels : 
            self.digitalMin[channel]            = int(f.read(8))
        for channel in self.channelLabels : 
            self.digitalMax[channel]            = int(f.read(8))  
        for channel in self.channelLabels : 
            self.prefiltering[channel]          = f.read(80).strip()            
        for channel in self.channelLabels : 
            self.nbSamplesPerRecord[channel]    = int(f.read(8))  

        f.read(32 * self.nbChannels)  # reserved
          
        assert(f.tell() == self.headerNbBytes)


        # Compute variables to locate pages
        #self.headerSize   = self.header['header_nbytes']    
        self.recordSize = sum(self.nbSamplesPerRecord.values())*self.nbBytes

        # Patch pour corriger le fait qu'Harmonie met -1 comme valeur à 
        # self.n_records.
        # Le nombre de bytes restant dans le fichier divisé par ns/sizeof(interger) 
        # va donner la valeur de records
        if self.nbRecords == -1:
            fileSize = os.path.getsize(fileName)
            self.nbRecords =  int((fileSize-self.headerNbBytes)/self.recordSize);      
    

        # calculate ranges for rescaling
        self.gain = {}
        for channel in self.channelLabels :         
            physicalRange = self.physicalMax[channel] - self.physicalMin[channel]
            digitalRange  = self.digitalMax[channel]  - self.digitalMin[channel]
            assert np.all(physicalRange > 0)
            assert np.all(digitalRange > 0)
            self.gain[channel] = physicalRange / digitalRange        
        






   
class EDFReader(EEGDBReaderBase) :
    def __init__(self, fname):
        #EEGDBReaderBase.__init__(self)
        
        self.fileName = fname
        self.file     = open(fname, 'rb')
        
        self.readHeader()    
        self.readEvents()        
        
        
    def getNbPages(self):    
        return self.header.nbRecords          
                
    def getChannelLabels(self):
        return [v for v in self.header.channelLabels if not v == EVENT_CHANNEL]        
        
        
    def readEvents(self):
        self.events          = []
        self.recordStartTime = []
        self.file.seek(self.header.headerNbBytes)        
        #indEventChannel = self.header.label.index(EVENT_CHANNEL)
        for noPage in range(self.getNbPages()):        
            rawRecord = self.readRawRecord()  
            tals      = tal(rawRecord[EVENT_CHANNEL])
            
            # The first index is the mendatory time keeping event. We don't need it.
            self.recordStartTime.append(EDFEvent((tals[0][0], tals[0][1], tals[0][2])))     
            
            for talEvent in tals[1:] : 
                # One TAL can contain many events wit the same startTime/Duration properties
                for noEventStr in talEvent[2]:
                    self.events.append(EDFEvent((talEvent[0], talEvent[1], noEventStr)))            
            




    def getChannelTime(self, channel) :
        if not isinstance(channel, str) :
            raise TypeError        
        
        #if not channel in self.getChannelLabels() :
        #    raise 
        
        nbSamples = float(self.header.nbSamplesPerRecord[channel])
        recordTime = arange(nbSamples)/nbSamples*self.header.recordDuration
        return concatenate([e.startTime + recordTime for e in self.recordStartTime])        
 

    """
     Convert a String of bytes to an array of integers.
    """
    def byteStr2integers(self, samples):
        if self.header.fileType == "EDF":
            # 2-byte little-endian integers
            return np.fromstring(samples, '<h')
        elif self.header.fileType == "BDF":
            #print "reading BDF..."
            # 3-byte little-endian integers
            
            # Inserting a fourth byte to read it as 4-byte floats
            # byteStr = "\0".join([samples[(noByte*3):(noByte*3+3)] for noByte in range(int(len(samples)/3))]) + "\0"
            # ==>> cannot use this because we do not always need to pad with zeros...                     
            
            # Pad witth zeros or ones depending on wheter the number is 
            # positive or negative
            byteStr = ""
            for noByte in range(int(len(samples)/3)):
                byteStr += samples[(noByte*3):(noByte*3+3)] + ('\0' if samples[noByte*3+2] < '\x80' else '\xff')                        
                                             
            return np.fromstring(byteStr, '<i')          
        else:
            raise IOError



    """
     Digital to physical conversion.
    """
    def digital2physical(self, dig, channel):
        return (dig - self.header.digitalMin[channel]) * self.header.gain[channel] + self.header.physicalMin[channel]




    """
     Read the complete signal recorded by a given channel.
    """
    def readChannel(self, signalName, usePickled=False):

        if not isinstance(signalName, str) :
            raise TypeError        
        
        
        if usePickled:
            raise NotImplementedError
        else:
            self.file.seek(self.header.headerNbBytes)        

            data = RecordedChannel()         
            data.samplingRate   = self.header.nbSamplesPerRecord[signalName]/self.header.recordDuration
            data.type           = "EEG"
            data.startTime      = self.header.startDateTime          
            
            signal = []
            for noPage in range(self.getNbPages()):                  
                dig = self.byteStr2integers(self.readRawRecord()[signalName])
                signal.extend(self.digital2physical(dig, signalName))
                
            data.signal = array(signal)
            return data    




    def getEvents(self): 
        return self.events
     
        
        
    def setPageDuration(self, duration):
        self.changeRecordDuration(duration)
        
    def getPageDuration(self):
        return self.header.recordDuration
                
        
    def readPage(self, channelList, pageId):

        # Position the file cursor to the begin of the page no. pageId :
        pagePosition = self.header.headerNbBytes  + (pageId-1)*self.header.recordSize
        self.file.seek(pagePosition)
   
        record = self.readRecord()
        sigStart = record[0]        
        
        recordedSignals = {}
        samplingRates   = {}
        for channel in channelList:
            samplingRates[channel]    = float(self.header.nbSamplesPerRecord[channel])/self.header.recordDuration
            recordedSignals [channel] = record[1][channel]
                
        return EEGPage(samplingRates, sigStart, recordedSignals, self.header.startDateTime)
        
        
        
        
    def getNbSample(self, channel):
        if not isinstance(channel, str):
            raise TypeError
                    
        assert(channel in self.getChannelLabels())
                                
        return self.header.nbSamplesPerRecord[channel]*self.header.nbRecords
               
        

    def __del__(self):
        self.file.close()


    def readHeader(self):
        self.header = EDFHeader(self.file, self.fileName)
      
      
      
    '''
     Read a record and return a dictionnary indexed by the signal labels
     containing arrays with raw bytes.  
    '''        
    def readRawRecord(self):

        result = {}
        
        for channel in self.header.channelLabels:
            nbBytes = self.header.nbSamplesPerRecord[channel]*self.header.nbBytes
            result[channel] = self.file.read(nbBytes)
            if len(result[channel]) != nbBytes :
                raise EOFError   
                
        return result


    '''
     Read a raw record, convert it to a (time, signals, events) tuple based on
     information in the header, and return it.
    '''    
    def readRecord(self):

        rawRecord = self.readRawRecord()    
        
        #dig_min, phys_min, gain = self.dig_min, self.phys_min, self.gain
            
        #offset_seconds = self.currentRecordInd*self.header.record_length      
        #time = self.header.date_time + datetime.timedelta(0,offset_seconds)
        signals = {}
        events = []
        for channel in rawRecord:
            if channel == EVENT_CHANNEL:
                ann = tal(rawRecord[channel])
                time = self.header.startDateTime + datetime.timedelta(0,ann[0][0]) 
                events.extend(ann[1:])
            else:
                dig = self.byteStr2integers(rawRecord[channel])
                signals[channel] = self.digital2physical(dig, channel)
        
        return time, signals, events



    def changeRecordDuration(self, newDuration):
        #TODO: Implement.
        raise NotImplemented




    


    def records(self):
        '''
        Record generator.
        '''
        try:
          while True:
            yield self.readRecord()
        except EOFError:
          pass


        
        