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
 We consider the record concept of EDF as equivalent to the page concept in
 Harmony and in sleep stage scoring in general. Thus, to use 30 second pages, 
 the EDF file must be reformatted if the record duration is different than
 30 seconds.
'''


from EEGDatabaseReader import EEGDBReaderBase, EEGPage

import uuid
import os, io
import re, datetime, logging
import numpy as np
import scipy.weave as weave

from scipy import array, arange, concatenate, where
from copy import deepcopy, copy
from lxml import etree
from tempfile import gettempdir 
from time import sleep
from warnings import warn
from shutil import copyfile

from spyndle.io import Event, RecordedChannel, EventList


EVENT_CHANNEL    = 'EDF Annotations'
REFORMAT_CHANNEL = 'EDF Reformat'
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
                    if not isinstance(name, unicode):
                        if isinstance(name, str):  
                            name = unicode(name)
                        else:
                            print type(name)
                            raise TypeError

                    if not isinstance(value, unicode):
                        if isinstance(value, str):  
                            value = unicode(value)
                        else:
                            print type(value)
                            raise TypeError
                    
                    
                    
                    if name == u"name":
                        self.name = value
                    elif name == u"groupName":
                        self.groupName = value
                    elif name == u"channel":
                        self.channel = value
                    else :
                        self.properties[name] = value                    
            except:
                print talEvent[2]
                raise
        else:
            self.name = talEvent[2]
            if ( self.name == u"Sleep stage 1" or self.name == u"Sleep stage 2" or
                 self.name == u"Sleep stage 3" or self.name == u"Sleep stage 4" or
                 self.name == u"Sleep stage R" or self.name == u"Sleep stage W" or
                 self.name == u"Sleep stage ?" ) :
               
               
               self.groupName = u"Stage"



"""
 Class used to reformat the EEG montage.

 Channel labes must be enclosed by []. Accepted operators so far are
 +, -, /, *. Pathenthesis can be used.
 E.g. "(2*[C3-A1]-[A2-A1])/2" can be used to obtain a reformatted ear-linked reference.
"""
class ReformatExpressionMng:
    def __init__(self, channelExpr, originalHeader):
        self.channelExpr    = channelExpr
        self.originalHeader = originalHeader
        
        for channel in self.channelExpr :
            exprChannels = self.getExprChannels(self.channelExpr[channel])
            if len(np.unique([originalHeader.nbSamplesPerRecord[chan] for chan in exprChannels])) != 1:
                print channel, self.channelExpr[channel], exprChannels, [originalHeader.nbSamplesPerRecord(chan) for chan in exprChannels]
                raise "All channels used in the expression of a reformatted "\
                       "channel must have the same sampling frequency"
    
    
    def getUsedChannels(self):
        channelLst = []
        for channel in self.channelExpr :
            channelLst.extend(self.getExprChannels(self.channelExpr[channel]) )
        return np.unique(channelLst)
    
    def getExprChannels(self, expression):
        return re.findall("\[(.*?)\]", expression) #"(2*[C3-A1]-[A2-A1])/2") 
        
    #def getExprNonChannels(self, expression):
     #   tokens = re.split("\[(.*?)\]", "(2*[C3-A1]-[A2-A1])/2") 
      #  return tokens[0:len(tokens):2]
                
                
    
    def physicalMinMax(self, channel):    

        expr         = self.channelExpr[channel]
        exprChannels = self.getExprChannels(expr)
        minMaxValues = np.array(np.concatenate([[self.originalHeader.physicalMin[exprChannel],
                                                  self.originalHeader.physicalMax[exprChannel]] 
                                                  for exprChannel in exprChannels]))
          
               
        candidates   = eval(re.sub("\[(.*?)\]", np.array_repr(minMaxValues, precision=32),  expr))
        return min(candidates), max(candidates)
      
      
    def physicalMin(self, channel):  
        return self.physicalMinMax(channel)[0]
        
    def physicalMax(self, channel):
        return self.physicalMinMax(channel)[1]
        
        
        
        
    
    def digitalMinMax(self, channel):    

        expr         = self.channelExpr[channel]
        exprChannels = self.getExprChannels(expr)
        minMaxValues = np.array(np.concatenate([[self.originalHeader.digitalMin[exprChannel],
                                                  self.originalHeader.digitalMax[exprChannel]] 
                                                  for exprChannel in exprChannels]))
          

        candidates   = eval(re.sub("\[(.*?)\]",np.array_repr(minMaxValues, precision=32),  expr))
        return min(candidates), max(candidates)
              
         
    def digitalMin(self, channel):
        return self.digitalMinMax(channel)[0]
        
    def digitalMax(self, channel):
        return self.digitalMinMax(channel)[1]
        
        
        
        
        
        
    def nbSamplesPerRecord(self, channel):
        exprChannels = self.getExprChannels(self.channelExpr[channel])
        return self.originalHeader.nbSamplesPerRecord[exprChannels[0]]  
        
        
        
        
    def getRawRecord(self, channel, page):
        
        expr         = self.channelExpr[channel]
        exprChannels = self.getExprChannels(expr)
        
        np.set_printoptions(threshold=np.nan)
        string = ""
        globals = {}
        novar = 0
        for token in re.split("\[(.*?)\]", expr) :  
            if token in exprChannels:
                globals["X" + str(novar)] = page.recordedSignals[token]
                string += "X" + str(novar)
                novar +=1
            else:
                string += token
                
        return eval(string, globals)
        





class EDFHeader :
    # BDF is a 24-bit version of the 16-bit EDF format, so we read it with the
    # the same reader
    def __init__(self, f=None, fileName=None):
        if not f is None and not fileName is None:
            self.read(f, fileName)


    def read(self, f, fileName):
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
        
        if self.fileType == "EDF":
            if self.subtype == 'EDF+C':             
                self.fileType   = "EDF+"
                self.contiguous = True
            elif self.subtype == 'EDF+D': 
                self.fileType   = "EDF+"
                self.contiguous = False           
            elif self.subtype == '     ': 
                self.fileType   = "EDF"
                self.contiguous = True   
                
        elif self.fileType == "BDF":
            self.contiguous = False           
 

             
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
        elif self.fileType == "EDF+" and self.contiguous:
            f.write("EDF+C" + " "*39)     
        elif self.fileType == "EDF+" and not self.contiguous:
            f.write("EDF+D" + " "*39)              
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




    def reformatHeader(self, reformattedChannels):
        originalHeader = deepcopy(self)        
        
        self.recordingIR += " (reformatted)"
        if len(self.recordingIR) > 80:
            self.recordingIR = self.recordingIR[0:80]
        
        
        self.reformattedChannels = reformattedChannels
           
 
        self.nbChannels         = len(reformattedChannels)+2
        
        # read channel info
        self.channelLabels      = reformattedChannels.keys()
        for i, channel in enumerate(self.channelLabels):
            if len(channel) > 16:
                self.channelLabels[i] = channel[:16]
        
        self.transducerType     = {}
        self.units              = {}
        self.physicalMin        = {}
        self.physicalMax        = {}          
        self.digitalMin         = {}
        self.digitalMax         = {}
        self.prefiltering       = {}          
        self.nbSamplesPerRecord = {} 
        
        for channel in self.channelLabels : 
            expressionMng                       = ReformatExpressionMng(reformattedChannels, originalHeader)
            self.transducerType[channel]        = "Reformatted"
            self.prefiltering[channel]          = "Reformatted"   
            self.units[channel]                 = "Reformat" 
            self.physicalMin[channel]           = expressionMng.physicalMin(channel)
            self.physicalMax[channel]           = expressionMng.physicalMax(channel)   
            self.digitalMin[channel]            = expressionMng.digitalMin(channel)
            self.digitalMax[channel]            = expressionMng.digitalMax(channel)
            self.nbSamplesPerRecord[channel]    = expressionMng.nbSamplesPerRecord(channel)     


        self.channelLabels.extend([EVENT_CHANNEL, REFORMAT_CHANNEL])

        self.transducerType[EVENT_CHANNEL]        = ""
        self.prefiltering[EVENT_CHANNEL]          = ""   
        self.units[EVENT_CHANNEL]                 = "" 
        self.physicalMin[EVENT_CHANNEL]           = -1
        self.physicalMax[EVENT_CHANNEL]           =  1
        self.digitalMin[EVENT_CHANNEL]            = -32768
        self.digitalMax[EVENT_CHANNEL]            =  32767
        self.nbSamplesPerRecord[EVENT_CHANNEL]    = 400

        self.transducerType[REFORMAT_CHANNEL]        = ""
        self.prefiltering[REFORMAT_CHANNEL]          = ""   
        self.units[REFORMAT_CHANNEL]                 = "" 
        self.physicalMin[REFORMAT_CHANNEL]           = -1
        self.physicalMax[REFORMAT_CHANNEL]           =  1
        self.digitalMin[REFORMAT_CHANNEL]            = -32768
        self.digitalMax[REFORMAT_CHANNEL]            =  32767
        self.nbSamplesPerRecord[REFORMAT_CHANNEL]    = len(str(self.reformattedChannels))/self.nbRecords/self.nbBytes + 1


        # calculate ranges for rescaling
        self.gain = {}
        for channel in self.channelLabels :         
            physicalRange = self.physicalMax[channel] - self.physicalMin[channel]
            digitalRange  = self.digitalMax[channel]  - self.digitalMin[channel]
            assert np.all(physicalRange > 0)
            assert np.all(digitalRange > 0)
            self.gain[channel] = physicalRange / digitalRange        
        
        self.headerNbBytes = 8 + 80 + 80 + 8 + 8 + 8 + 44 + 8 + 8 + 4 + \
                              self.nbChannels*(16 + 80+ 8 + 8 + 8 + 8 + 8 + 80 + 8 + 32)







   
class EDFBaseReader(EEGDBReaderBase) :

    def __init__(self, fname, isSplitted=None, annotationFileName=None):        
        self.fileName = fname
        
        with io.open(fname, 'rb') as fileObj:
            self.header = EDFHeader(fileObj, self.fileName)  
        
            super(EDFBaseReader, self).__init__(self.getPageDuration()) 
            
            if self.header.fileType == "EDF+" or self.header.fileType == "BDF" :
                self.readEvents(fileObj)        
        

    def getFileName(self): 
        return self.fileName       
        
    def getRecordingStartTime(self): 
        return self.header.startDateTime        
        
    def getNbPages(self):    
        return self.header.nbRecords          
                
    def getChannelLabels(self):
        return [v for v in self.header.channelLabels if not v == EVENT_CHANNEL and not v == REFORMAT_CHANNEL]        
        
        
    def readEvents(self, fileObj):
        #self.events          = []
        self.recordStartTime = []        
        fileObj.seek(self.header.headerNbBytes)        

        #indEventChannel = self.header.label.index(EVENT_CHANNEL)
        for noPage in range(self.getNbPages()):            
            rawRecord = self.readRawRecord(fileObj)  
            tals      = tal(rawRecord[EVENT_CHANNEL])        

            # The first index is the mendatory time keeping event. We record it separately.
            # The time duration of this time keeping event is left blank but we know that records
            # are of a duration given by self.header.recordDuration
            self.recordStartTime.append(EDFEvent((tals[0][0], self.header.recordDuration, ""))) 

            
            for talEvent in tals[1:] : 
                # One TAL can contain many events wit the same startTime/Duration properties
                for noEventStr in talEvent[2]:
                    self.events.add(EDFEvent((talEvent[0], talEvent[1], noEventStr)))  
                    

            
            

    def addEvent(self, event):
        self.events.add(event)   






    def reformatMontage(self, channelExpressions, saveFileName = ""):
        if saveFileName == "":
            saveFileName = self.fileName[:-4] + "_reformatted" + self.fileName[-4:]
            
            
        reformattedReader = deepcopy(self)
        reformattedReader.fileName = saveFileName
        reformattedReader.header.reformatHeader(channelExpressions)


        with io.open(saveFileName, 'wb') as fileWrite:
            with io.open(self.fileName, 'rb') as fileObj:
            
                # Computing the strings representing the EDF annotations
                eventStings     = reformattedReader.computeEDFAnnotations()   
                reformatStrings = reformattedReader.computeReformatStrings()

                
                # Verifying if the annotation field is large enough to record the annotations.
                # If not, enlarge it on the writing copy.
                nbBytesEvent    = reformattedReader.header.nbSamplesPerRecord[EVENT_CHANNEL]*reformattedReader.header.nbBytes
                nbBytesReformat = reformattedReader.header.nbSamplesPerRecord[REFORMAT_CHANNEL]*reformattedReader.header.nbBytes
                if max([len(s) for s in eventStings]) >= nbBytesEvent:       
                    reformattedReader.header.nbSamplesPerRecord[EVENT_CHANNEL] = int(max([len(s) for s in eventStings])*1.2/reformattedReader.header.nbBytes)
                    nbBytesEvent    = reformattedReader.header.nbSamplesPerRecord[EVENT_CHANNEL]*reformattedReader.header.nbBytes

               
                # Write the header
                reformattedReader.header.write(fileWrite)
                fileObj.seek(self.header.headerNbBytes)
        
                expressionMng = ReformatExpressionMng(reformattedReader.header.reformattedChannels, self.header) 
                for nopage in range(self.getNbPages()):   
                    
                    page = self.readPage(expressionMng.getUsedChannels(), nopage+1) 
                        
                    for channel in reformattedReader.header.channelLabels:
                        if channel == EVENT_CHANNEL:    
                            encodedString = eventStings[nopage].encode("utf8")  
                            fileWrite.write(encodedString + "\0"*(nbBytesEvent - len(encodedString)) )  
                        elif channel == REFORMAT_CHANNEL:
                            encodedString = reformatStrings[nopage].encode("utf8")                     

                            fileWrite.write(encodedString + "\0"*(nbBytesReformat - len(encodedString)) )  
                        else: 
                                                        
                            formatedSignal = expressionMng.getRawRecord(channel, page)
                            phys_range = reformattedReader.header.physicalMax[channel]  - reformattedReader.header.physicalMin[channel]
                            dig_range = reformattedReader.header.digitalMax[channel] - reformattedReader.header.digitalMin[channel]
                            assert np.all(phys_range > 0)
                            assert np.all(dig_range > 0)
                            gain = dig_range/phys_range                          
                            
                            formatedSignal = (formatedSignal -  reformattedReader.header.physicalMin[channel])*gain + reformattedReader.header.digitalMin[channel]       
        
                            if reformattedReader.header.nbBytes == 2: # EDF
                                fileWrite.write(formatedSignal.astype('<h').tostring())
                            elif reformattedReader.header.nbBytes == 3: #BDF
                                # Writing in a string of 32-bit integers and removing every fourth byte
                                # to obtain a string of 24-bit integers
                                formatedSignal = formatedSignal.astype('<i').tostring()
                                formatedSignal = "".join([formatedSignal[noBit] for noBit in range(len(formatedSignal)) if (noBit+1)%4])
                                
                                fileWrite.write(formatedSignal)


        return reformattedReader






    def crop(self, timeStart, timeEnd, saveFileName=None):
        
        
        """
         Note: timeStart and timeEnd are not the exact starting and ending time.
         The starting and ending time will be those of the page containing the
         samples corresponding to timeStart and timeEnd. We avoid to make new
         record divisions to obtain the exact starting and ending time since the
         records can be discontinuous in EDF+ format.
        """
        
        # Valid pages
        noPages = [noPage for noPage, rec in enumerate(self.recordStartTime) if \
                     rec.startTime <= timeEnd and rec.startTime + rec.timeLength >= timeStart] 

        if saveFileName is None:
            self.save(noPages=noPages)
        else:
            self.saveAs(saveFileName, noPages=noPages)            



    """
     Save the header information and the events to the file. The data themself
     are unchanged for two reasons:
         1 - We should always keep the original data as is an apply modifications
             (e.g., filtering) "on-line" to avoid loosing important information
             in the recorded phenomena.
         2 - The data are not loaded in kept within the object, as opposed to 
             the header and the event informations.
             
     noPage can either be left to known or to be set to a list of page number. 
     In the latter case, only the page with numbers included in that list will
     be saved.
     
    """
    def save(self, tempPath=None, noPages = None):

        if tempPath is None:
            tempPath = gettempdir() 
 
        # Because the annotation field may change in size, records can be shifted 
        # so we cannot only alter the information in the annotation fields. We need
        # to create a completely new file and swap it with the original.
        tempFileName = tempPath + "temp-" + str(uuid.uuid1()) + "-" + os.path.basename(self.fileName)
                
        self.saveAs(tempFileName, noPages=noPages)

        # Try to delete the original file. This may failled because the file
        # is used by another process (e.g., when running )        
        for i in range(13):
            try:
                os.remove(self.fileName)
            except WindowsError:
                warn("Failed to remove " + self.fileName + ". Retrying in " + str(2**i) + " seconds.", UserWarning)
                sleep(2**i)
                continue
            break            
             
        # rename the new file       
        for i in range(13):
            try:
                os.rename(tempFileName, self.fileName)  
            except WindowsError:
                warn("Failed to remove " + self.fileName + ". Retrying in " + str(2**i) + " seconds.", UserWarning)
                sleep(2**i)
                continue
            break          
        
        # Reinit the object with the new file.
        self.__init__(self.fileName)


    def saveAs(self, saveFileName, noPages=None):

        with io.open(saveFileName, 'wb') as fileWrite:
            with io.open(self.fileName, 'rb') as fileObj:
            
                # Computing the strings representing the EDF annotations
                eventStings = self.computeEDFAnnotations()       
                               
                # Any modifications to the header for the writing must be made in a copied
                # header as it must not alter the reading of the file for the duration
                # of the saving operation.
                writeHeader = deepcopy(self.header)
                
                if not noPages is None:
                    writeHeader.nbRecords = len(noPages)                 
                
                
                # Verifying if the annotation field is large enough to record the annotations.
                # If not, enlarge it on the writing copy.
                if max([len(s) for s in eventStings]) >= self.header.nbSamplesPerRecord[EVENT_CHANNEL]*self.header.nbBytes:       
                    writeHeader.nbSamplesPerRecord[EVENT_CHANNEL] = int(max([len(s) for s in eventStings])/self.header.nbBytes*1.2)
        
               
                # Write the header
                writeHeader.write(fileWrite)
                fileObj.seek(fileWrite.tell())

                # If it has not been specified that only certain pages are to be
                # save, then specify that all pages are to be saved.
                if noPages is None: 
                    noPages = range(len(eventStings))
        
                # Write the body. Actually, we leave the data intact, updating only
                # the annotation field.
                for noPage, eventStr in enumerate(eventStings):
                    rawRecord = self.readRawRecord(fileObj)    
                    if noPage in noPages:                    
                        for channel in self.header.channelLabels:
                            if channel == EVENT_CHANNEL:    
                                encodedString = eventStr.encode("utf8")                     
                                fileWrite.write(encodedString + "\0"*(writeHeader.nbSamplesPerRecord[channel]*self.header.nbBytes - len(encodedString)) )  
                            else: 
                                fileWrite.write(rawRecord[channel])
            
                
                #assert(fileWrite.tell() == writeHeader.headerNbBytes + sum(writeHeader.nbSamplesPerRecord.values())*writeHeader.nbBytes*self.getNbPages())








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


        # Computing the strings representing the EDF annotations
        eventStings = []
        nbEvents = 0
        for nopage, startTimeEvent in enumerate(self.recordStartTime):
            startTime = startTimeEvent.startTime
                    
            if nopage == 0:                                      # first page 
                filteredEvents = filter(lambda e: e.startTime < self.recordStartTime[1].startTime, self.events)     
                
            elif nopage == len(self.recordStartTime) - 1:       # last page
                filteredEvents = filter(lambda e: e.startTime >= startTime, self.events)  
                
            else:                                               # other pages
                filteredEvents = filter(lambda e: e.startTime <  self.recordStartTime[nopage+1].startTime and 
                                                  e.startTime >= startTime, self.events)     
                
            eventStings.extend([prepareEventStr(startTime, filteredEvents)])         

            nbEvents += len(filteredEvents)  

        # Verify that the number of recorded events correspond to the number of events within the object.                          
        assert(nbEvents == len(self.events))    

        return eventStings





    def computeReformatStrings(self):
        
        # Computing the strings representing the reformating operations
        formatStings = []
        if hasattr(self.header, "reformattedChannels"):
            computeString = str(self.header.reformattedChannels)
        else:
            computeString = ""            
        
        nbSamplesRecord = self.header.nbSamplesPerRecord[REFORMAT_CHANNEL]*self.header.nbBytes       
        
        indStart = 0
        for nopage in range(self.getNbPages()):
            indEnd = min(indStart + nbSamplesRecord, len(computeString))
            formatStings.append(computeString[indStart:indEnd]) 
            indStart = indEnd

        return formatStings









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
        if self.header.fileType == "EDF" or  self.header.fileType == "EDF+":
            # 2-byte little-endian integers
            return np.fromstring(samples, '<i2').astype(np.int32)
        elif self.header.fileType == "BDF":
            #print "reading BDF..."
            # 3-byte little-endian integers
            
            # Inserting a fourth byte to read it as 4-byte floats
            # byteStr = "\0".join([samples[(noByte*3):(noByte*3+3)] for noByte in range(int(len(samples)/3))]) + "\0"
            # ==>> cannot use this because we do not always need to pad with zeros...                     
            # Pad witth zeros or ones depending on wheter the number is 
            # positive or negative

            ##### This implementation is about 75% slower than            
            #byteStr = ""                                            
            #for noByte in range(0, len(samples), 3):
            #    byteStr += samples[noByte:(noByte+3)] + ('\0' if samples[noByte+2] < '\x80' else '\xff')                        
            #return np.fromstring(byteStr, '<i')  
            # A ticket has been placed for the implementation of a '<3i' dtype in numpy (https://github.com/numpy/numpy/issues/664)

            ##### this one... 
            #import struct
            #return array([struct.unpack('<i', samples[noByte:(noByte+3)] + ('\0' if samples[noByte+2] < '\x80' else '\xff'))[0] for noByte in range(0, len(samples), 3)])
            
            ##### which is about ten time slower than this one
            samples = np.fromstring(samples, '<u1')
            
            N   = len(samples)
            ret = np.zeros(N/3, dtype='int32')
            
            code = """
                        int i, i2;
                        for (i=0, i2=0; i<N; i+=3, i2++)
                            ret(i2) = ((samples(i+2) < 128 ? 0 : 255) << 24) | (samples(i+2) << 16) | (samples(i+1) << 8) | samples(i);
                        return_val = 1;
                 """  
            
            weave.inline(code, [ 'N', 'ret', 'samples'], type_converters=weave.converters.blitz, compiler = 'gcc')
            return ret            
            
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

        if not (isinstance(signalName, str) or isinstance(signalName, unicode)) :
            print "Wrong type:", type(signalName), signalName
            raise TypeError        
        
        
        if usePickled:
            raise NotImplementedError
        else:
            with io.open(self.fileName, 'rb') as fileObj:
                fileObj.seek(self.header.headerNbBytes)        

                data = RecordedChannel()         
                data.samplingRate   = self.header.nbSamplesPerRecord[signalName]/self.header.recordDuration
                data.type           = "EEG"
                data.startTime      = self.header.startDateTime          
                
                signal = []
                for noPage in range(self.getNbPages()):                  
                    dig = self.byteStr2integers(self.readRawRecord(fileObj)[signalName])
                    signal.extend(self.digital2physical(dig, signalName))
                    
                data.signal = array(signal)
                return data    



    """
     Return the time associated with the next sample following startTime. Return
     None if there is no sample next to startTime.
    """
    def getNextSampleStartTime(self, signalName, startTime):
             
        for noPage, startTimeEvent in enumerate(self.recordStartTime):
            pageStartTime = startTimeEvent.startTime
            pageDuration  = startTimeEvent.timeLength
            
            if pageStartTime + pageDuration >= startTime:
                
                pageNbRecords = self.header.nbSamplesPerRecord[signalName]
                samplingRate  = float(self.header.nbSamplesPerRecord[signalName])/self.header.recordDuration
                
                time = pageStartTime + arange(pageNbRecords)/samplingRate
                #time = page.getStartTime() + arange(len(page.recordedSignals[signalName]))/page.samplingRates[signalName]
                IND = where(time >= startTime)[0]
                if len(IND) :
                    return time[IND[0]]                      

        return None 


    def read(self, signalNames, startTime, timeDuration, debug=False):
                
        pages = []
        for noPage, startTimeEvent in enumerate(self.recordStartTime):
            pageStartTime = startTimeEvent.startTime
            pageDuration  = startTimeEvent.timeLength
            
            if pageStartTime + pageDuration >= startTime:
                pages.append(self.readPage(signalNames, noPage+1))
            
            if pageStartTime + pageDuration >= startTime + timeDuration:    
                break

            # Caution: The pages need to be in order.

        assert(len(pages)>0)

        info = pages[0]
        if len(pages) > 1:
            for page in pages[1:]:
                for key in page.recordedSignals:
                    info.recordedSignals[key] = concatenate((info.recordedSignals[key], page.recordedSignals[key]))


        if debug:
            channel = info.recordedSignals.keys()[0]
            time = info.getStartTime() + arange(len(info.recordedSignals[channel]))/info.samplingRates[channel]
            print "Reader : ", len(pages), startTime, timeDuration, time, info.samplingRates[channel]




        returnData = {}  
        for channel in info.recordedSignals:
            time = info.getStartTime() + arange(len(info.recordedSignals[channel]))/info.samplingRates[channel]
            ind  = where((time >= startTime)*(time <= startTime + timeDuration))[0]
            
            try:
                assert(len(ind)>0)
            except :
                print "time: ", time
                print info.getStartTime(), len(info.recordedSignals[channel]), info.samplingRates[channel]
                print startTime, timeDuration
                raise


            returnData[channel]                = RecordedChannel()               
            returnData[channel].samplingRate   = info.samplingRates[channel]
            returnData[channel].type           = None
            returnData[channel].signal         = info.recordedSignals[channel][ind]
            returnData[channel].startTime      = time[ind[0]]      
    
        return returnData

        
    def setPageDuration(self, duration):
        nbSamples = copy(self.header.nbSamplesPerRecord)
        del nbSamples[EVENT_CHANNEL]
        
        # len(nbSamples) == 0 for edfa files.
        if len(nbSamples):
            dtMax = self.header.recordDuration/max(array(nbSamples.values()))
            # We change the duration only if its greater than half the size of the 
            # smaller sampling period. Else, because of the digitization resolution,
            # the change has no effet. We remove the EVENT_CHANNEL from the computation
            # of the sampling period because the sampling period of this channel has
            # no meaning.
            if abs(duration - self.getPageDuration()) > dtMax:
                #print dtMax, max(array(self.header.nbSamplesPerRecord.values())), self.header.recordDuration
                #print nbSamples.values()
                #print nbSamples        
                #print abs(duration - self.getPageDuration()), duration, self.getPageDuration()
                self.changeRecordDuration(duration)
        
    def getPageDuration(self):
        return self.header.recordDuration
                
        
    def readPage(self, channelList, pageId):
        # PageId are numbered starting from 1, not from 0.

        # Position the file cursor to the begin of the page no. pageId :
        pagePosition = self.header.headerNbBytes  + (pageId-1)*self.header.recordSize
        with io.open(self.fileName, 'rb') as fileObj:
            fileObj.seek(pagePosition)
            record = self.readRecord(fileObj)
            
        sigStart = record[0]        
        
        recordedSignals = {}
        samplingRates   = {}
        for channel in channelList:
            samplingRates[channel]    = float(self.header.nbSamplesPerRecord[channel])/self.header.recordDuration
            recordedSignals [channel] = record[1][channel]
                
        return EEGPage(samplingRates, sigStart, recordedSignals, self.header.startDateTime)
        
        
        
        
    def getNbSample(self, channel=None):

        if channel is None:
            return max(self.header.nbSamplesPerRecord.values())*self.header.nbRecords
        
        if not isinstance(channel, str):
            raise TypeError
                    
        assert(channel in self.getChannelLabels())
                                
        return self.header.nbSamplesPerRecord[channel]*self.header.nbRecords
               
        


      
      
      
    '''
     Read a record and return a dictionnary indexed by the signal labels
     containing arrays with raw bytes.  
    '''        
    def readRawRecord(self, fileObj):

        result = {}
        
        for channel in self.header.channelLabels:
            nbBytes = self.header.nbSamplesPerRecord[channel]*self.header.nbBytes
            result[channel] = fileObj.read(nbBytes)
            if len(result[channel]) != nbBytes :
                raise EOFError  
                
        return result


    '''
     Read a raw record, convert it to a (time, signals, events) tuple based on
     information in the header, and return it.
    '''    
    def readRecord(self, fileObj):

        rawRecord = self.readRawRecord(fileObj)    
        
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
            elif channel == REFORMAT_CHANNEL:
                signals[channel] = rawRecord[channel]
            else:
                dig = self.byteStr2integers(rawRecord[channel])                                
                signals[channel] = self.digital2physical(dig, channel)
        
        return time, signals, events



    def changeRecordDuration(self, newDuration):
        #TODO: Implement.
        raise NotImplemented






   
   
   
   
   
class EDFReader(EEGDBReaderBase) :
    
    """
     When working with large datafiles such as whole night polysomnographic 
     recording, recording the annotations in the same file as the recording
     data can be inneficient since adding new annotations might requires 
     reformatting the whole data file. In these case, it can be more interesting
     to split the original data and the annotations in two separate files. 
     
     We suggest to always use the same name for both files, except for a 
     supplementary "a" (for annotation) at the end of the extension of the 
     annotation file (e.g. aFile.bdf/aFile.bdfa or someFile.edf/someFile.edfa).
     Nevertheless, this class provide the option of using an arbitrary name 
     for the annotation file.

     The parameter isSplitted (boolean) can be used to specify wheter the 
     data file specified by fname is a splitted data/annotation set of files or
     a single file. If isSplitted == True, the annotationFileName
     parameter will specify the name of the annotation file.
     If not specified, fname + "a" will be used.
     If isSplitted == False, no supplementary annotation file is considered
     regardless of the value of annotationFileName. If isSplitted == None, then
     if an annotationFileName is specified, this annotation file will be used; if
     no annotationFileName is specified but there exist a  fname + "a" file, this
     annotation file will be used; if no annotationFileName is specified and
     there is no file fname + "a", no anotation file will be used.
     
     Note that if an annotation file is used, only the annotations in the 
     annotation file be available. These of the data file (if any) will be 
     masked.
    """
    def __init__(self, fname, isSplitted=None, annotationFileName=None):        
        
        if isSplitted is False:
            self.isSplitted = False
            self.annotationFileName = None
            
        elif isSplitted is True:
            self.isSplitted = True            
            if not annotationFileName is None:
                self.annotationFileName = annotationFileName
            elif os.path.isfile(fname + "a") :
                self.annotationFileName = fname + "a"
            else:
                raise "In EDFReader(...), the parameter isSplitted is set to True but"\
                      " no annotation file name is specified using the annotationFileName"\
                      " parameter and no file " + fname + "a is existing."
                      
        elif isSplitted is None:
            if not annotationFileName is None:
                self.annotationFileName = annotationFileName
                self.isSplitted         = True 
            elif os.path.isfile(fname + "a") :
                self.annotationFileName = fname + "a"
                self.isSplitted         = True 
            else:
                self.isSplitted         = False  
                self.annotationFileName = None
    
        else:
            raise "In EDFReader(...), the parameter isSplitted can only take the "\
                  "values True, False or None. Value " + str(isSplitted) + " used."                   


        self.dataReader        = EDFBaseReader(fname)
        super(EDFReader, self).__init__(self.getPageDuration()) 

        if self.isSplitted :
            self.annotationsReader  = EDFBaseReader(self.annotationFileName)
            self.events             = self.annotationsReader.events
        else:
            self.annotationsReader = None
            self.events             = self.dataReader.events
        

    @property
    def header(self):
        return self.dataReader.header

    @header.setter
    def header(self, value):
        self.dataReader.header = value


    @property
    def fileName(self):
        return self.dataReader.fileName

    @fileName.setter
    def fileName(self, value):
        self.dataReader.fileName = value



    def getFileName(self): 
        #warn("EDFReader.getFileName(...) is deprecated. Please use EDFReader.fileName()", DeprecationWarning)
        return self.fileName       
        
    def getRecordingStartTime(self): 
        return self.header.startDateTime        
        
    def getNbPages(self):    
        return self.header.nbRecords          
                
    def getChannelLabels(self):
        return self.dataReader.getChannelLabels()

    def addEvent(self, event):
        self.events.add(event)   


    def reformatMontage(self, channelExpressions, saveFileName = "", anotationSaveFileName=""):

        if saveFileName == "":
            saveFileName = self.fileName[:-4] + "_reformatted" + self.fileName[-4:]


        # If the original recording was not a splited file and no file name
        # has been given for the annotation file, do not split the file.
        if anotationSaveFileName == "" and self.annotationFileName is None:
            pass

        # If the original recording was a splited file but no file name
        # has been given for the annotation file, use the standard naming 
        # for annotation file (exact same name as the data file but with
        # a "a" suffix to the file extension).        
        elif anotationSaveFileName == "" and not self.annotationFileName is None:
            anotationSaveFileName = saveFileName + "a"
            copyfile(self.annotationFileName, anotationSaveFileName)


        # If the original recording was not a splited file but a file name
        # has been given for the annotation file, split the reformated record.   
        elif anotationSaveFileName != "" and self.annotationFileName is None:
            self.splitRecord(anotationSaveFileName)
            self.save()

        # If the original recording was  a splited file and a file name
        # has been given for the annotation file, save the annotation file using
        # the requested name.           
        else: # anotationSaveFileName != "" and not self.annotationFileName is None
            copyfile(self.annotationFileName, anotationSaveFileName)
        
          

        self.dataReader.reformatMontage(channelExpressions, saveFileName)        
                 
         
         
    def splitRecord(self, annotationFileName):

        self.isSplitted = True
        self.annotationFileName = annotationFileName          

        self.annotationsReader  = deepcopy(self.dataReader)
        self.events             = self.annotationsReader.events

        # TODO: Events in the self.dataReader could be removed (except for the time
        # stamps of the records which are mendatory for EDF+).
        #self.dataReader ...

        self.annotationsReader.header.headerNbBytes  = 8 + 80 + 80 + 8 + 8 + 8 + 44 + 8 + 8 + 4 + (16 + 80+ 8 + 8 + 8 + 8 + 8 + 80 + 8 + 32) 
        self.annotationsReader.header.nbChannels     = 1
        self.annotationsReader.header.channelLabels  = ["EDF Annotations"] 

        self.annotationsReader.header.transducerType     = {"EDF Annotations":""}
        self.annotationsReader.header.prefiltering       = {"EDF Annotations":""} 
        self.annotationsReader.header.units              = {"EDF Annotations":""}
        self.annotationsReader.header.physicalMin        = {"EDF Annotations":-1}
        self.annotationsReader.header.physicalMax        = {"EDF Annotations":1}          
        self.annotationsReader.header.digitalMin         = {"EDF Annotations":-32768}
        self.annotationsReader.header.digitalMax         = {"EDF Annotations":32767}     
        
        # TODO: The nbSamplesPerRecord could be adjusted.
        #self.annotationsReader.header.nbSamplesPerRecord = {"EDF Annotations":annotationFieldLength} 
            
  
  
  



    def crop(self, timeStart, timeEnd, saveFileName=None, saveAnnotationFileName=None):
        
        if saveAnnotationFileName is None and not saveFileName is None:
            saveAnnotationFileName = saveFileName + "a"

        if self.isSplitted:
            self.annotationsReader.crop(timeStart, timeEnd, saveAnnotationFileName)
        self.dataReader.crop(timeStart, timeEnd, saveFileName)


            


    """
     Save the header information and the events to the file. The data themself
     are unchanged for two reasons:
         1 - We should always keep the original data as is and apply modifications
             (e.g., filtering) "on-line" to avoid loosing important information
             in the recorded phenomena.
         2 - The data are not loaded and kept within the object, as opposed to 
             the header and the event informations.
     However, if the record is not split edf/edfa files, changes in events may
     require increasing the space of the event records, forcing to offset 
     appropriately each data record.
    """
    def save(self, tempPath=None):

        """
         If the EEG recording is a splitted data file set, we just need to 
         save the annotation file. Else, we need to save the complete data file.
        """
        if self.isSplitted:
            self.annotationsReader.save(tempPath)
        else:
            self.dataReader.save(tempPath)




    def saveAs(self, saveFileName, saveAnnotationFileName=""):

        if saveAnnotationFileName == "":
            saveAnnotationFileName = saveFileName + "a"

        if self.isSplitted:
            copyfile(self.fileName, saveFileName)
            self.annotationsReader.saveAs(saveAnnotationFileName)
        else:
            self.dataReader.saveAs(saveFileName)
            
            


    def getChannelFreq(self, channel): 
        return self.dataReader.getChannelFreq(channel)


    def getChannelTime(self, channel) :
        return self.dataReader.getChannelTime(channel)   
 


    """
     Read the complete signal recorded by a given channel.
    """
    def readChannel(self, signalName, usePickled=False):
        return self.dataReader.readChannel(signalName, usePickled)           



    """
     Return the time associated with the next sample following startTime. Return
     None if there is no sample next to startTime.
    """
    def getNextSampleStartTime(self, signalName, startTime):
        return self.dataReader.getNextSampleStartTime(signalName, startTime)                 
             

    def read(self, signalNames, startTime, timeDuration, debug=False):
        return self.dataReader.read(signalNames, startTime, timeDuration, debug)                   
          
        
    def setPageDuration(self, duration):
        self.dataReader.setPageDuration(duration)        

        
    def getPageDuration(self):
        return self.dataReader.getPageDuration()
                
        
    def readPage(self, channelList, pageId):
        return self.dataReader.readPage(channelList, pageId)        
        
        
    def getNbSample(self, channel=None):
        return self.dataReader.getNbSample(channel)     




