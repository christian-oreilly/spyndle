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


import os, io
import re, datetime, logging
import numpy as np

from scipy import array, arange, concatenate, where
from copy import deepcopy
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
            try:
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
            except:
                print talEvent[2]
                raise
        else:
            self.name = talEvent[2]
            if ( self.name == "Sleep stage 1" or self.name == "Sleep stage 2" or
                 self.name == "Sleep stage 3" or self.name == "Sleep stage 4" or
                 self.name == "Sleep stage R" or self.name == "Sleep stage W" or
                 self.name == "Sleep stage ?" ) :
               
               
               self.groupeName = "Stage"


   



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
        


    def write(self, f):
        
        # Position the cursor at the start of the file
        f.seek(0)

        #8 ascii : version of this data format (0) 
        if self.fileType == "EDF":
            f.write("0       ")
        elif self.fileType == "BDF":
            f.write("\xFFBIOSEMI")
        
        #80 ascii : local patient identification 
        #The 'local patient identification' field must start with the subfields 
        #(subfields do not contain, but are separated by, spaces):
        #  - the code by which the patient is known in the hospital administration.
        #  - sex (English, so F or M).
        #  - birthdate in dd-MMM-yyyy format using the English 3-character abbreviations 
        #    of the month in capitals. 02-AUG-1951 is OK, while 2-AUG-1951 is not.
        #  - the patients name.
        #Any space inside the hospital code or the name of the patient must be replaced by
        #a different character, for instance an underscore. For instance, the 'local patient 
        #identification' field could start with: MCH-0234567 F 02-MAY-1951 Haagse_Harry. 
        #Subfields whose contents are unknown, not applicable or must be made anonymous 
        #are replaced by a single character 'X'. Additional subfields may follow the ones described here.                
        f.write(self.subjectID  + (80-len(self.subjectID ))*" ")
          
          

        # 80 ascii : local recording identification 
        # The 'local recording identification' field must start with the 
        # subfields (subfields do not contain, but are separated by, spaces):
        #    - The text 'Startdate'.
        #    - The startdate itself in dd-MMM-yyyy format using the English 
        #      3-character abbreviations of the month in capitals.
        #    - The hospital administration code of the investigation, 
        #      i.e. EEG number or PSG number.
        #    - A code specifying the responsible investigator or technician.
        #    - A code specifying the used equipment.
        # Any space inside any of these codes must be replaced by a different 
        # character, for instance an underscore. The 'local recording identification' 
        # field could contain: Startdate 02-MAR-2002 PSG-1234/2002 NN Telemetry03. 
        # Subfields whose contents are unknown, not applicable or must be made anonymous
        #  are replaced by a single character 'X'. So, if everything is unknown then the 
        # 'local recording identification' field would start with: Startdate X X X X. 
        # Additional subfields may follow the ones described here.      
        f.write(self.recordingIR + (80-len(self.recordingIR))*" ")              
          
        #8 ascii : startdate of recording (dd.mm.yy)
        f.write(self.startDateTime.strftime("%d.%m.%y"))
             
        #8 ascii : starttime of recording (hh.mm.ss) 
        f.write(self.startDateTime.strftime("%H.%M.%S"))
         

        # 8 ascii : number of bytes in header record
        f.write("%08d" % self.headerNbBytes)   
        
        # 44 ascii : reserved
        if self.fileType == "EDF":
            f.write(" "*44) 
        elif self.fileType == "BDF":
            f.write("24BIT" + " "*39)            
        
        
        # 8 ascii : number of data records (-1 if unknown)
        #  The 'number of data records' can only be -1 during recording. 
        # As soon as the file is closed, the correct number is known and must be entered. 
        f.write("%08d" % self.nbRecords)  
        
        # 8 ascii : duration of a data record, in seconds
        f.write(("%8.6f" % self.recordDuration)[:8] )  
        
        # 4 ascii : number of signals (ns) in data record
        f.write("%04d" % (self.nbChannels))
        
        # ns * 16 ascii : ns * label (e.g. EEG Fpz-Cz or Body temp)
        for channel in self.channelLabels: 
            f.write("%16s" % channel)


  
        # ns * 80 ascii : ns * transducer type (e.g. AgAgCl electrode)
        for channel in self.channelLabels : 
            f.write("%80s" % self.transducerType[channel])
            
        # ns * 8 ascii : ns * physical dimension (e.g. uV or degreeC)
        for channel in self.channelLabels : 
            f.write("%8s" % self.units[channel]  )  
  
        # ns * 8 ascii : ns * physical minimum (e.g. -500 or 34)
        for channel in self.channelLabels : 
            f.write(("%8.6f" %  self.physicalMin[channel])[:8])  
  
        # ns * 8 ascii : ns * physical maximum (e.g. 500 or 40)
        for channel in self.channelLabels : 
            f.write(("%8.6f" %  self.physicalMax[channel])[:8])
 
        # ns * 8 ascii : ns * digital minimum (e.g. -2048)
        for channel in self.channelLabels : 
            f.write("%08d" %  self.digitalMin[channel])

        # ns * 8 ascii : ns * digital maximum (e.g. 2047)
        for channel in self.channelLabels : 
            f.write("%08d" %  self.digitalMax[channel])  
  
        # ns * 80 ascii : ns * prefiltering (e.g. HP:0.1Hz LP:75Hz)
        for channel in self.channelLabels : 
            f.write("%80s" %  self.prefiltering[channel])  
        
        
        # ns * 8 ascii : ns * nr of samples in each data record
        for channel in self.channelLabels : 
            f.write("%08d" % self.nbSamplesPerRecord[channel])


        # ns * 32 ascii : ns * reserved
        f.write(" "*32*len(self.channelLabels))

        assert(f.tell() == self.headerNbBytes)










   
class EDFReader(EEGDBReaderBase) :
    def __init__(self, fname):        
        self.fileName = fname
        self.file     = io.open(fname, 'rb')
        
        self.header = EDFHeader(self.file, self.fileName)  
        
        super(EDFReader, self).__init__(self.getPageDuration())        
        self.readEvents()        
        
        
        
    def getRecordingStartTime(self): 
        return self.header.startDateTime        
        
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
                       
            #print rawRecord[EVENT_CHANNEL].decode('utf-8')            
            
            # The first index is the mendatory time keeping event. We record it separately.
            self.recordStartTime.append(EDFEvent((tals[0][0], tals[0][1], tals[0][2]))) 
         
            
            for talEvent in tals[1:] : 
                # One TAL can contain many events wit the same startTime/Duration properties
                for noEventStr in talEvent[2]:
                    self.events.append(EDFEvent((talEvent[0], talEvent[1], noEventStr)))            
                    

            
            

    def addEvent(self, event):
        self.events.append(event)   



    """
     Save the header information and the events to the file. The data themself
     are unchanged for two reasons:
         1 - We should always keep the original data as is an apply modifications
             (e.g., filtering) "on-line" to avoid loosing important information
             in the recorded phenomena.
         2 - The data are not loaded in kept within the object, as opposed to 
             the header and the event informations.
    """
    def save(self):

 
        # Because the annotation field may change in size, records can be shifted 
        # so we cannot only alter the information in the annotation fields. We need
        # to create a completely new file and swap it with the original. 
        self.fileWrite  = io.open("temp", 'wb')
        
        # Computing the strings representing the EDF annotations
        eventStings = self.computeEDFAnnotations()       
                       
        # Any modifications to the header for the writing must be made in a copied
        # header as it must not alter the reading of the file for the duration
        # of the saving operation.
        writeHeader = deepcopy(self.header)
        
        # Verifying if the annotation field is large enough to record the annotations.
        # If not, enlarge it on the writing copy.
        if max([len(s) for s in eventStings]) >= self.header.nbSamplesPerRecord[EVENT_CHANNEL]*self.header.nbBytes:       
            writeHeader.nbSamplesPerRecord[EVENT_CHANNEL] = int(max([len(s) for s in eventStings])*1.2)

       
        # Write the header
        writeHeader.write(self.fileWrite)

        # Write the body. Actually, we leave the data intact, updating only
        # the annotation field.
        for eventStr in eventStings:
            rawRecord = self.readRawRecord()    
            for channel in self.header.channelLabels:
                if channel == EVENT_CHANNEL:
                    self.fileWrite.write(eventStr + "\0"*(writeHeader.nbSamplesPerRecord[channel]*self.header.nbBytes - len(eventStr))  )  
                else: 
                    self.fileWrite.write(rawRecord[channel])

        
        assert(self.fileWrite.tell() == writeHeader.headerNbBytes + sum(writeHeader.nbSamplesPerRecord.values())*writeHeader.nbBytes*self.getNbPages())

        self.file.close()  
        self.fileWrite.close()   
        
        # Delete the original file
        os.remove(self.fileName)
 
        # rename the new file       
        os.rename("temp", self.fileName)     
        
        # Reinit the object with the new file.
        self.__init__(self.fileName)



    """ 
     Compute the EDFAnnotations strings for each pages (records) from the 
     events contrained in the reader. Used internally by self.save(...)
    """
    def computeEDFAnnotations(self):

        def prepareEventStr(timeDiff, events):
            if not isinstance(timeDiff, float):
                raise TypeError
            if not isinstance(events, list):
                raise TypeError                
            
            
            # Annotation channel            
            timeKeepingStr = "+" + str(timeDiff) + "\x14\x14\0"       
                           
            for event in events:
                if not isinstance(event, Event):
                    raise TypeError        

                timeKeepingStr += event.toEDFStr()


            return timeKeepingStr


        # Position the file cursor just after the header.
        self.file.seek(self.header.headerNbBytes)

        # Computing the strings representing the EDF annotations
        eventStings = []
        nbEvents = 0
        for nopage, startTimeEvent in enumerate(self.recordStartTime):
            startTime = startTimeEvent.startTime
            for channel in self.header.channelLabels:
                if channel == EVENT_CHANNEL:
                    
                    if nopage == 0:
                        filteredEvents = filter(lambda e: e.startTime < self.recordStartTime[1].startTime, self.events)     
                        
                    elif nopage == len(self.recordStartTime) - 1:
                        filteredEvents = filter(lambda e: e.startTime >= startTime, self.events)  
                        
                    else:                    
                        filteredEvents = filter(lambda e: e.startTime <  self.recordStartTime[nopage+1].startTime and 
                                                          e.startTime >= startTime, self.events)     
                        
                    eventStings.extend([prepareEventStr(startTime, filteredEvents)])         

                    nbEvents += len(filteredEvents)  

        # Verify that the number of recorded events correspond to the number of events within the object.                          
        assert(nbEvents == len(self.events))    

        return eventStings



    def getChannelFreq(self, channel): 
        return float(self.header.nbSamplesPerRecord[channel])/self.header.recordDuration        


    def getChannelTime(self, channel) :
        if not (isinstance(channel, str) or isinstance(channel, unicode)):
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


    def read(self, signalNames, startTime, timeDuration):
                
        pages = []
        for noPage, startTimeEvent in enumerate(self.recordStartTime):
            pageStartTime = startTimeEvent.startTime
            pageDuration  = startTimeEvent.timeLength
            
            if pageStartTime + pageDuration >= startTime:
                pages.append(self.readPage(signalNames, noPage))
            
            if pageStartTime + pageDuration >= startTime + timeDuration:    
                break

            # Caution: The pages need to be in order.

        assert(len(pages)>0)

        info = pages[0]
        if len(pages) > 1:
            for page in pages[1:]:
                for key in page.recordedSignals:
                    info.recordedSignals[key] = concatenate((info.recordedSignals[key], page.recordedSignals[key]))

        returnData = {}  
        for channel in info.recordedSignals:
            time = info.getStartTime() + arange(len(info.recordedSignals[channel]))/info.samplingRates[channel]
            ind  = where((time >= startTime)*(time <= startTime + timeDuration))[0]
            assert(len(ind)>0)


            returnData[channel]                = RecordedChannel()               
            returnData[channel].samplingRate   = info.samplingRates[channel]
            returnData[channel].type           = None
            returnData[channel].signal         = info.recordedSignals[channel][ind]
            returnData[channel].startTime      = time[ind[0]]      
    
        return returnData








    def getEvents(self): 
        return self.events
     
        
        
    def setPageDuration(self, duration):
        if duration != self.getPageDuration():
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


        
        