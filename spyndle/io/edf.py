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
#import pyximport
#pyximport.install(setup_args={"include_dirs": numpy.get_include()},
#                  reload_support=False)
from .xEdf import encodeBDF



"""
class EDFConf:
    __useWeave = True
    
    @property
    def useWeave(self):
        return self.__useWeave

    @useWeave.setter        
    def useWeave(self, useWeave):
        self.__useWeave = useWeave
        


try:
    with open('edfConf.pkl', 'rb') as pkl_file:
        edfConf = pickle.load(pkl_file)
        if not hasattr(edfConf, "useWeave"):
            raise ValueError
except:
    edfConf = EDFConf()    
"""


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
  exp = b'(?P<onset>[+\-]\d+(?:\.\d*)?)'       + \
        b'(?:\x15(?P<duration>\d+(?:\.\d*)?))?' + \
        b'(\x14(?P<annotation>[^\x00]*))?'     + \
        b'(?:\x14\x00)'

  def annotation_to_list(annotation):
    return annotation.split(b'\x14') if annotation else []

  def parse(dic):
    return (
      round(float(dic['onset']), 6), 
      float(dic['duration']) if dic['duration'] else 0.,
      annotation_to_list(dic['annotation']))

  return [parse(m.groupdict()) for m in re.finditer(exp, tal_str)]


  
class EDFEvent(Event):
    def __init__(self, talEvent=None): #, ISignalFile):
        Event.__init__(self)
        if talEvent is None:
            return
            
        self.startTime   = talEvent[0]
        self.timeLength  = talEvent[1]        # Duration in seconds  
        #self.dateTime    = ""          


        self.properties  = {}
        # XML event
        self.label = talEvent[2].decode('latin-1')
        if self.label[:6] == "<Event" :
            try:
                root = etree.fromstring(self.label)
                for name, value in sorted(root.items()):
                    if not isinstance(name, str):
                        if isinstance(name, str):  
                            name = str(name)
                        else:
                            print((type(name)))
                            raise TypeError

                    if not isinstance(value, str):
                        if isinstance(value, str):  
                            value = str(value)
                        else:
                            print((type(value)))
                            raise TypeError
                    
                    
                    
                    if name == "name":
                        self.name = value
                    elif name == "groupName":
                        self.groupName = value
                    elif name == "channel":
                        self.channel = value
                    elif name == "cycle":
                        self.cycle = int(value)
                    else :
                        self.properties[name] = value                    
            except:
                print(( self.label))
                raise
        else:
            self.name =  self.label
            if ( self.name == "Sleep stage 1" or self.name == "Sleep stage 2" or
                 self.name == "Sleep stage 3" or self.name == "Sleep stage 4" or
                 self.name == "Sleep stage R" or self.name == "Sleep stage W" or
                 self.name == "Sleep stage ?" ) :
               
               
               self.groupName = "Stage"
            else:
                raise(ValueError("Unrecognized event format (" + self.label + ")."))


    def initFromData_sample(self, reader, startSample, 
                              durationSample, name, groupName,
                              channel = "", properties={}):

        self.startTime   = reader.timeFromSample(channel, startSample)
        self.timeLength  = durationSample/float(reader.getChannelFreq(channel))
        
        self.properties  = properties
        self.name = name
        self.groupName = groupName
        self.channel = channel
                  









"""
 Class used to reformat the EEG montage.

 Channel labes must be enclosed by []. Accepted operators so far are
 +, -, /, *. Pathenthesis can be used.
 E.g. "(2*[C3-A1]-[A2-A1])/2" can be used to obtain a reformatted ear-linked reference.
"""
class ReformatExpressionMng:
    def __init__(self, channelExpr, originalHeader, transducerTypes = None,
                       prefilterings = None, units = None):
        
        self.__transducerTypes = transducerTypes
        self.__prefilterings = prefilterings
        self.__units = units

        self.channelExpr    = channelExpr
        self.originalHeader = originalHeader
        
        for channel in self.channelExpr :
            exprChannels = self.getExprChannels(self.channelExpr[channel])
            if len(np.unique([originalHeader.nbSamplesPerRecord[chan] for chan in exprChannels])) != 1:
                print((channel, self.channelExpr[channel], exprChannels, [originalHeader.nbSamplesPerRecord(chan) for chan in exprChannels]))
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
          
        """               
        candidates   = eval(re.sub("\[(.*?)\]", np.array_repr(minMaxValues, precision=32),  expr))
        return min(candidates), max(candidates)
        """
        return min(minMaxValues), max(minMaxValues)
      
      
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
          
        """  
        candidates   = eval(re.sub("\[(.*?)\]",np.array_repr(minMaxValues, precision=32),  expr))
        return min(candidates), max(candidates)
        """      
        return min(minMaxValues), max(minMaxValues)
         
    def digitalMin(self, channel):
        return self.digitalMinMax(channel)[0]
        
    def digitalMax(self, channel):
        return self.digitalMinMax(channel)[1]
        
        
    def transducerType(self, channel):
        if self.__transducerTypes is None:
            return "Reformatted"
        elif isinstance(self.__transducerTypes, str) :
            
            if self.__transducerTypes == "keep":
                """
                 We keep the tranducerType of the first channel in the expression.
                """
                expr         = self.channelExpr[channel]
                exprChannels = self.getExprChannels(expr)  
                return self.originalHeader.transducerType[exprChannels[0]] 
                        
                
            else:
                ValueError("Invalid value for self.transducerTypes.")
                
        else:
            ValueError("Invalid value for self.transducerTypes.")

        
    def prefiltering(self, channel):
        if self.__prefilterings is None:
            return "Reformatted"
        elif isinstance(self.__prefilterings, str) :
            
            if self.__prefilterings == "keep":                   
                """
                 We keep the prefiltering of the first channel in the expression.
                """
                expr         = self.channelExpr[channel]
                exprChannels = self.getExprChannels(expr)  
                return self.originalHeader.prefiltering[exprChannels[0]]
                
            else:
                ValueError("Invalid value for self.prefilterings.")
                
        else:
            ValueError("Invalid value for self.prefilterings.")
                    
        
    def units(self, channel):
        if self.__units is None:
            return "Reformat"
        
        elif isinstance(self.__units, str) :
            if self.__units == "keep":                   
                """
                 We keep the units of the first channel in the expression.
                """
                expr         = self.channelExpr[channel]
                exprChannels = self.getExprChannels(expr)  
                return self.originalHeader.units[exprChannels[0]]                
                
            else:
                ValueError("Invalid value for self.units.")
                
        else:
            ValueError("Invalid value for self.units.")
        
        
        
    def nbSamplesPerRecord(self, channel):
        exprChannels = self.getExprChannels(self.channelExpr[channel])
        assert(isinstance(self.originalHeader.nbSamplesPerRecord[exprChannels[0]], int))
        return self.originalHeader.nbSamplesPerRecord[exprChannels[0]]  
        
        
        
        
    def getRawRecord(self, channel, record):
        
        expr         = self.channelExpr[channel]
        exprChannels = self.getExprChannels(expr)
        
        np.set_printoptions(threshold=np.nan)
        string = ""
        globals = {}
        novar = 0
        for token in re.split("\[(.*?)\]", expr) :  
            if token in exprChannels:
                globals["X" + str(novar)] = record.recordedSignals[token]
                string += "X" + str(novar)
                novar +=1
            else:
                string += token
                
        return eval(string, globals)
        

class IntDict(dict):
    def __setitem__(self,key,value):
        if isinstance(value, int):
            dict.__setitem__(self,key,value)
        else:
            raise TypeError("This dictionnary takes only interger values. Value "\
                            + "passed is of type " + str(type(value)))



def rename_duplicated_channels(channels):
    """
     We don't allow for duplicated channel names. First because it poorly 
     identifies the signals, second because it crash thre reader because of 
     the use of dictionnaries indexed by channel names.
    """
    for item in channels:
        if channels.count(item) > 1:
            assert(channels.count(item) < 10)
            no =1
            for i in range(len(channels)):
                if channels[i] == item:
                    channels[i] = channels[i][:14] + "_" + str(no)
                    no += 1        
    return channels

class EDFHeader :
    # BDF is a 24-bit version of the 16-bit EDF format, so we read it with the
    # the same reader
    def __init__(self, f=None, fileName=None, eventChannel=EVENT_CHANNEL):
        if not f is None and not fileName is None:
            self.read(f, fileName)
            self.fileName = fileName
            self.eventChannel = eventChannel

    def anonymize(self):         
        self.subjectID = "X X X X"
        self.recordingIR = "Startdate X X X X"
        
        # Guessing that the data are less than 100 years old....
        self.startDateTime = datetime.datetime(2000, 1, 1, self.startDateTime.hour, 
                                               self.startDateTime.minute, self.startDateTime.second)
        
        
        
        
        
        
    def read(self, f, fileName):
        assert f.tell() == 0  # check file position

        tampon = f.read(8)
        if tampon == b'0       ':
             self.fileType = "EDF"
             self.nbBytes  = 2     # 16 bits
        elif  tampon == b'\xffBIOSEMI':
             self.fileType = "BDF"
             self.nbBytes  = 3     # 24 bits
        else:
             raise IOError("Invalid file header. First character found: " +  
                           str(tampon) + ". Expected values are '0       ' or '\xFFBIOSEMI'.")

        
        # recording info
        self.subjectID = f.read(80).decode('latin-1').strip()
        self.recordingIR = f.read(80).decode('latin-1').strip()
        
        
        # parse timestamp
        (day, month, year) = [int(x) for x in re.findall('(\d+)', f.read(8).decode('latin-1'))]
        (hour, minute, sec)= [int(x) for x in re.findall('(\d+)', f.read(8).decode('latin-1'))]

        try:        
            # Guessing that the data are less than 100 years old....
            year = year + (2000 if year + 2000 < datetime.date.today().year else 1900) 
            self.startDateTime = datetime.datetime(year, month, day, hour, minute, sec)
        except ValueError:
            if month > 12:
                warnings.warn("Found a month larger than 12 for the recording date. " + 
                              "Considering that the month and the day in the date "   +
                              "format have been inverted. Correcting.")
                self.startDateTime = datetime.datetime(year, day, month, hour, minute, sec)
            else:
                raise

        # misc
        self.headerNbBytes      = int(f.read(8))
        self.subtype            = f.read(44).decode('latin-1')[:5]
        
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
        self.channelLabels      = [f.read(16).decode('latin-1').strip() for n in range(self.nbChannels)]
        self.channelLabels      = rename_duplicated_channels(self.channelLabels)
        

        self.transducerType     = {}
        self.units              = {}
        self.physicalMin        = {}
        self.physicalMax        = {}          
        self.digitalMin         = {}
        self.digitalMax         = {}
        self.prefiltering       = {}          
        self.nbSamplesPerRecord = IntDict()
        self.samplingRate       = {}
        
        for channel in self.channelLabels : 
            self.transducerType[channel]        = f.read(80).decode('latin-1').strip()
        for channel in self.channelLabels : 
            self.units[channel]                 = f.read(8).decode('latin-1').strip()
        for channel in self.channelLabels : 
            self.physicalMin[channel]           = float(f.read(8))
        for channel in self.channelLabels : 
            self.physicalMax[channel]          = float(f.read(8))            
        for channel in self.channelLabels : 
            self.digitalMin[channel]            = int(f.read(8))
        for channel in self.channelLabels : 
            self.digitalMax[channel]            = int(f.read(8))  
        for channel in self.channelLabels : 
            self.prefiltering[channel]          = f.read(80).decode('latin-1').strip()           
        for channel in self.channelLabels : 
            self.nbSamplesPerRecord[channel]    = int(f.read(8))  
        for channel in self.channelLabels : 
            self.samplingRate[channel]          = float(self.nbSamplesPerRecord[channel])/float(self.recordDuration)


        f.read(32 * self.nbChannels)  # reserved
          
        assert(f.tell() == self.headerNbBytes)

        f.seek(0)
        self.headerString = f.read(self.headerNbBytes).decode('latin-1')


        # Compute variables to locate records
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
            #try:
            #    assert np.all(physicalRange > 0)
            #    assert np.all(digitalRange > 0)
            #except:
            #    print(physicalRange, digitalRange)
            #    print(self.physicalMax)
            #    print(self.physicalMin)
            #    print(self.digitalMax)
            #    print(self.digitalMin)
            self.gain[channel] = physicalRange / float(digitalRange)        
        


    def write(self, f):
        
        # Position the cursor at the start of the file
        f.seek(0)

        #8 ascii : version of this data format (0) 
        if self.fileType[:3] == "EDF":
            f.write(b"0       ")
        elif self.fileType == "BDF":
            f.write(b"\xFFBIOSEMI")

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
        f.write((self.subjectID  + (80-len(self.subjectID ))*" ").encode('latin-1'))
          
          

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
        f.write((self.recordingIR + (80-len(self.recordingIR))*" ").encode('latin-1'))
          
        #8 ascii : startdate of recording (dd.mm.yy)
        f.write(self.startDateTime.strftime("%d.%m.%y").encode('latin-1'))
             
        #8 ascii : starttime of recording (hh.mm.ss) 
        f.write(self.startDateTime.strftime("%H.%M.%S").encode('latin-1'))
         

        # 8 ascii : number of bytes in header record
        f.write(("%-8d" % self.headerNbBytes).encode('latin-1'))
        
        # 44 ascii : reserved
        if self.fileType == "EDF":
            fileTypeStr = " "*44 
        elif self.fileType == "EDF+" and self.contiguous:
            fileTypeStr = "EDF+C" + " "*39   
        elif self.fileType == "EDF+" and not self.contiguous:
            fileTypeStr = "EDF+D" + " "*39             
        elif self.fileType == "BDF":
            fileTypeStr = "24BIT" + " "*39   
        else:
            raise ValueError("Error in writing EDF header: self.fileType == '" + self.fileType + "' and self.contiguous == " + str(self.contiguous))
        
        f.write(fileTypeStr.encode('latin-1'))
        
        
        # 8 ascii : number of data records (-1 if unknown)
        #  The 'number of data records' can only be -1 during recording. 
        # As soon as the file is closed, the correct number is known and must be entered. 
        f.write(("%-8d" % self.nbRecords).encode('latin-1'))  
        
        # 8 ascii : duration of a data record, in seconds
        f.write((("%8.6f" % self.recordDuration)[:8]).encode('latin-1'))  
        
        # 4 ascii : number of signals (ns) in data record
        f.write(("%-4d" % (self.nbChannels)).encode('latin-1'))
        
        # ns * 16 ascii : ns * label (e.g. EEG Fpz-Cz or Body temp)
        for channel in self.channelLabels: 
            if len(channel) > 16:
                warningStr = "ChannelLabels larger than 16 characters. "\
                             "Truncating " + channel + " to " + channel[:16] +  "."
                warnings.warn(warningStr)
            f.write(("%-16s" % channel[:16]).encode('latin-1'))


  
        # ns * 80 ascii : ns * transducer type (e.g. AgAgCl electrode)
        for channel in self.channelLabels : 
            f.write(("%-80s" % self.transducerType[channel]).encode('latin-1'))
            
        # ns * 8 ascii : ns * physical dimension (e.g. uV or degreeC)
        for channel in self.channelLabels : 
            f.write(("%-8s" % self.units[channel]  ).encode('latin-1'))
  
        # ns * 8 ascii : ns * physical minimum (e.g. -500 or 34)
        for channel in self.channelLabels : 
            f.write((("%8.6f" %  self.physicalMin[channel])[:8]).encode('latin-1'))
  
        # ns * 8 ascii : ns * physical maximum (e.g. 500 or 40)
        for channel in self.channelLabels : 
            f.write((("%8.6f" %  self.physicalMax[channel])[:8]).encode('latin-1'))
 
        # ns * 8 ascii : ns * digital minimum (e.g. -2048)
        for channel in self.channelLabels : 
            f.write(("%-8d" %  self.digitalMin[channel]).encode('latin-1'))

        # ns * 8 ascii : ns * digital maximum (e.g. 2047)
        for channel in self.channelLabels : 
            f.write(("%-8d" %  self.digitalMax[channel]).encode('latin-1'))
  
        # ns * 80 ascii : ns * prefiltering (e.g. HP:0.1Hz LP:75Hz)
        for channel in self.channelLabels : 
            f.write(("%-80s" %  self.prefiltering[channel]).encode('latin-1'))
        
        
        # ns * 8 ascii : ns * nr of samples in each data record
        for channel in self.channelLabels : 
            f.write(("%-8d" % self.nbSamplesPerRecord[channel]).encode('latin-1'))


        # ns * 32 ascii : ns * reserved
        f.write((" "*32*len(self.channelLabels)).encode('latin-1'))

        assert(f.tell() == self.headerNbBytes)
        
        
        

    def reformatHeader(self, reformattedChannels, transducerTypes = None,
                       prefilterings = None, units = None):
        originalHeader = deepcopy(self)        
        
        self.recordingIR += " (reformatted)"
        if len(self.recordingIR) > 80:
            self.recordingIR = self.recordingIR[0:80]
        
        
        self.reformattedChannels = reformattedChannels
           
 
        self.nbChannels         = len(reformattedChannels)+2
        
        # read channel info
        self.channelLabels      = list(reformattedChannels.keys())
        self.channelLabels      = rename_duplicated_channels(self.channelLabels)
                
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
        self.nbSamplesPerRecord = IntDict()
        
        for channel in self.channelLabels : 
            expressionMng   = ReformatExpressionMng(reformattedChannels, originalHeader, 
                                                    transducerTypes, prefilterings, units)
                                                    
            self.transducerType[channel]        = expressionMng.transducerType(channel)
            self.prefiltering[channel]          = expressionMng.prefiltering(channel) 
            self.units[channel]                 = expressionMng.units(channel)                             
            self.physicalMin[channel]           = expressionMng.physicalMin(channel)
            self.physicalMax[channel]           = expressionMng.physicalMax(channel)   
            self.digitalMin[channel]            = expressionMng.digitalMin(channel)
            self.digitalMax[channel]            = expressionMng.digitalMax(channel)
            self.nbSamplesPerRecord[channel]    = expressionMng.nbSamplesPerRecord(channel)     


        self.channelLabels.extend([self.eventChannel, REFORMAT_CHANNEL])

        self.transducerType[self.eventChannel]        = ""
        self.prefiltering[self.eventChannel]          = ""   
        self.units[self.eventChannel]                 = "" 
        self.physicalMin[self.eventChannel]           = -1
        self.physicalMax[self.eventChannel]           =  1
        self.digitalMin[self.eventChannel]            = -32768
        self.digitalMax[self.eventChannel]            =  32767
        self.nbSamplesPerRecord[self.eventChannel]    = 400


        # We store in a new signal the reformattedChannels dict for reference.
        self.transducerType[REFORMAT_CHANNEL]        = ""
        self.prefiltering[REFORMAT_CHANNEL]          = ""   
        self.units[REFORMAT_CHANNEL]                 = "" 
        self.physicalMin[REFORMAT_CHANNEL]           = -1
        self.physicalMax[REFORMAT_CHANNEL]           =  1
        self.digitalMin[REFORMAT_CHANNEL]            = -32768
        self.digitalMax[REFORMAT_CHANNEL]            =  32767
        self.nbSamplesPerRecord[REFORMAT_CHANNEL]    = int(len(str(self.reformattedChannels))/self.nbRecords/self.nbBytes + 1)


        # calculate ranges for rescaling
        self.gain = {}
        for channel in self.channelLabels :         
            physicalRange = self.physicalMax[channel] - self.physicalMin[channel]
            digitalRange  = self.digitalMax[channel]  - self.digitalMin[channel]
            assert np.all(physicalRange > 0)
            assert np.all(digitalRange > 0)
            self.gain[channel] = physicalRange / float(digitalRange)       
        
        self.headerNbBytes = 8 + 80 + 80 + 8 + 8 + 8 + 44 + 8 + 8 + 4 + \
                              self.nbChannels*(16 + 80+ 8 + 8 + 8 + 8 + 8 + 80 + 8 + 32)


    def changeChannelSet(self, chanDict):
        """
         This method can be used to rename channels and to keep only a 
         subset of the available channels. The argument is dictionnary of 
         channels where the key is the current channel name and the value
         is the new channel name. Channel
         labels that are not in the chanDict will be removed.
        """
        originalHeader = deepcopy(self)        

        self.channelLabels      = [new for old, new in list(chanDict.items()) if old in originalHeader.channelLabels]       

        self.transducerType     = {}
        self.units              = {}
        self.physicalMin        = {}
        self.physicalMax        = {}          
        self.digitalMin         = {}
        self.digitalMax         = {}
        self.prefiltering       = {}          
        self.nbSamplesPerRecord = IntDict()
        
        for originalChannelLabel, newChannelKey in list(chanDict.items()) : 
            if originalChannelLabel in originalHeader.channelLabels:
                self.transducerType[newChannelKey]        = originalHeader.transducerType[originalChannelLabel]  
                self.prefiltering[newChannelKey]          = originalHeader.prefiltering[originalChannelLabel]   
                self.units[newChannelKey]                 = originalHeader.units[originalChannelLabel]  
                self.physicalMin[newChannelKey]           = originalHeader.physicalMin[originalChannelLabel]  
                self.physicalMax[newChannelKey]           = originalHeader.physicalMax[originalChannelLabel]    
                self.digitalMin[newChannelKey]            = originalHeader.digitalMin[originalChannelLabel]  
                self.digitalMax[newChannelKey]            = originalHeader.digitalMax[originalChannelLabel]  
                self.nbSamplesPerRecord[newChannelKey]    = originalHeader.nbSamplesPerRecord[originalChannelLabel]     
                self.gain[newChannelKey]                  = originalHeader.gain[originalChannelLabel]  

        self.channelLabels.append(self.eventChannel)

        self.transducerType[self.eventChannel]        = ""
        self.prefiltering[self.eventChannel]          = ""   
        self.units[self.eventChannel]                 = "" 
        self.physicalMin[self.eventChannel]           = -1
        self.physicalMax[self.eventChannel]           =  1
        self.digitalMin[self.eventChannel]            = -32768
        self.digitalMax[self.eventChannel]            =  32767
        self.nbSamplesPerRecord[self.eventChannel]    = originalHeader.nbSamplesPerRecord[self.eventChannel]

        self.nbChannels         = len(self.channelLabels)

        self.headerNbBytes = 8 + 80 + 80 + 8 + 8 + 8 + 44 + 8 + 8 + 4 + \
                              self.nbChannels*(16 + 80+ 8 + 8 + 8 + 8 + 8 + 80 + 8 + 32)



   
class EDFBaseReader(EEGDBReaderBase) :

    def __init__(self, fileName, isSplitted=None, annotationFileName=None, 
                 readEventsOnInit=True, eventChannel=EVENT_CHANNEL):        
        self.fileName = fileName
        self.eventChannel = eventChannel
        
        with io.open(fileName, 'rb') as fileObj:
            self.header = EDFHeader(fileObj, self.fileName)  
        
            super(EDFBaseReader, self).__init__(self.getRecordDuration()) 
            
            if readEventsOnInit:
                if self.header.fileType == "EDF+" or self.header.fileType == "BDF" :
                    self.readEvents(fileObj)
                else:
                    for noRecord in range(self.getNbRecords()):
                        self.recordsInfo.add(EEGRecordInfo(startTime=noRecord*self.getRecordDuration(), 
                                                 duration=self.getRecordDuration()))                    
        

    def getFileName(self): 
        return self.fileName       
        
    def getRecordingStartTime(self): 
        return self.header.startDateTime        
        
    def getNbRecords(self):    
        return self.header.nbRecords          
                
    def getChannelLabels(self):
        return [v for v in self.header.channelLabels if not v == self.eventChannel and not v == REFORMAT_CHANNEL]        
        
        
    def readEvents(self, fileObj):    
        fileObj.seek(self.header.headerNbBytes)        

        #indEventChannel = self.header.label.index(self.eventChannel)
        for noRecord in range(self.getNbRecords()):            
            rawRecord = self.readRawRecord(fileObj)  
            tals      = tal(rawRecord[self.eventChannel])     
           
            # The first index is the mendatory time keeping event. We record it separately.
            # The time duration of this time keeping event is left blank but we know that records
            # are of a duration given by self.header.recordDuration
            self.recordsInfo.add(EEGRecordInfo(startTime=tals[0][0], 
                                                 duration=self.header.recordDuration))

            for talEvent in tals[1:] : 
                # One TAL can contain many events wit the same startTime/Duration properties
                for noEventStr in talEvent[2]:
                    self.events.add(EDFEvent((talEvent[0], talEvent[1], noEventStr)))  
                    

        #print "readEvents", len(self.recordsInfo)
         
            




    def getEventStrings(self):
        """
        This methods is mostly for debugging. It reads and return the plain
        string representation of the event description.
        """

        if not (self.header.fileType == "EDF+" or self.header.fileType == "BDF") :
            return []

        with io.open(self.fileName, 'rb') as fileObj:
            fileObj.seek(self.header.headerNbBytes)        
    
            for noRecord in range(self.getNbRecords()):            
                rawRecord = self.readRawRecord(fileObj)  
                tals      = tal(rawRecord[self.eventChannel])        
                
                for talEvent in tals[1:] : 
                    # One TAL can contain many events wit the same startTime/Duration properties
                    for noEventStr in talEvent[2]:
                       print((talEvent[0], talEvent[1], noEventStr)) 





    def reformatMontage(self, channelExpressions, saveFileName = "", 
                        transducerTypes = None, prefilterings = None, units = None):
        if saveFileName == "":
            saveFileName = self.fileName[:-4] + "_reformatted" + self.fileName[-4:]
            
            
        reformattedReader = deepcopy(self)
        reformattedReader.fileName = saveFileName
        reformattedReader.header.reformatHeader(channelExpressions, transducerTypes, 
                                                prefilterings, units)

        if saveFileName == self.fileName:
            replace = True
            saveFileName += "_tmp"
        else:
            replace = False


        with io.open(saveFileName, 'wb') as fileWrite:
            with io.open(self.fileName, 'rb') as fileObj:
            
                # Computing the strings representing the EDF annotations
                eventStings     = reformattedReader.computeEDFAnnotations()   
                reformatStrings = reformattedReader.computeReformatStrings()

                
                # Verifying if the annotation field is large enough to record the annotations.
                # If not, enlarge it on the writing copy.
                nbBytesEvent    = reformattedReader.header.nbSamplesPerRecord[self.eventChannel]*reformattedReader.header.nbBytes
                nbBytesReformat = reformattedReader.header.nbSamplesPerRecord[REFORMAT_CHANNEL]*reformattedReader.header.nbBytes
                if max([len(s) for s in eventStings]) >= nbBytesEvent:       
                    reformattedReader.header.nbSamplesPerRecord[self.eventChannel] = int(max([len(s) for s in eventStings])*1.2/float(reformattedReader.header.nbBytes))
                    nbBytesEvent    = reformattedReader.header.nbSamplesPerRecord[self.eventChannel]*reformattedReader.header.nbBytes

               
                # Write the header
                reformattedReader.header.write(fileWrite)
                fileObj.seek(self.header.headerNbBytes)
        
                expressionMng = ReformatExpressionMng(reformattedReader.header.reformattedChannels, self.header) 
                for noRecord in range(self.getNbRecords()):   
                    
                    record = self.readRecord(expressionMng.getUsedChannels(), noRecord+1) 
                        
                    for channel in reformattedReader.header.channelLabels:
                        if channel == self.eventChannel:    
                            encodedString = eventStings[noRecord].encode("latin-1")  
                            fileWrite.write(encodedString + b'\x00'*(nbBytesEvent - len(encodedString)) )  
                        elif channel == REFORMAT_CHANNEL:
                            encodedString = reformatStrings[noRecord].encode("latin-1")                     

                            fileWrite.write(encodedString + b'\x00'*(nbBytesReformat - len(encodedString)) )  
                        else: 
                                                        
                            formatedSignal = expressionMng.getRawRecord(channel, record)
                            phys_range = reformattedReader.header.physicalMax[channel]  - reformattedReader.header.physicalMin[channel]
                            dig_range = reformattedReader.header.digitalMax[channel] - reformattedReader.header.digitalMin[channel]
                            assert np.all(phys_range > 0)
                            assert np.all(dig_range > 0)
                            gain = dig_range/float(phys_range)                        
                            
                            formatedSignal = (formatedSignal -  reformattedReader.header.physicalMin[channel])*gain + reformattedReader.header.digitalMin[channel]       
        
                            if reformattedReader.header.nbBytes == 2: # EDF
                                fileWrite.write(formatedSignal.astype('<h').tostring())
                            elif reformattedReader.header.nbBytes == 3: #BDF
                                # Writing in a string of 32-bit integers and removing every fourth byte
                                # to obtain a string of 24-bit integers
                                formatedSignal = formatedSignal.astype('<i').tostring()
                                formatedSignal = "".join([formatedSignal[noBit] for noBit in range(len(formatedSignal)) if (noBit+1)%4])
                                
                                fileWrite.write(formatedSignal)

        if replace:
            os.remove(self.fileName)
            move(saveFileName, self.fileName)

        return reformattedReader







    def changeChannelSet(self, chanDict, saveFileName = ""):
        """
         This method can be used to rename channels and to keep only a 
         subset of the available channels. The argument is dictionnary of 
         channels where the key is the current channel name and the value
         is the new channel name. Channel
         labels that are not in the chanDict will be removed.
        """
        if saveFileName == "":
            isTemp = True
            saveFileName = os.path.join(tempfile.gettempdir(), "renamed.edf")
        else:
            isTemp = False
            
        reversedDict = {value:key for (key, value) in list(chanDict.items())}            
            
        changedReader = deepcopy(self)
        changedReader.fileName = saveFileName
        changedReader.header.changeChannelSet(chanDict)


        with io.open(saveFileName, 'wb') as fileWrite:
            with io.open(self.fileName, 'rb') as fileRead:
            
                ######
                # This might be set as an optionnal operation. It cannot be
                # performed systematically since we do not want for example
                # that an annotation file not containing the EEG signals 
                # sees its events be removed!
                ######
                # Change channels name of events and remove events from 
                # droped channels.
                #eventToRemove = []
                #for event in changedReader.events:
                #    if hasattr(event, "channel"):
                #        if event.channel in chanDict:
                #            event.channel = chanDict[event.channel]
                #        elif event.channel != "":
                #            eventToRemove.append(event)
                #            
                #for event in eventToRemove:
                #    changedReader.removeEvent(event)

                eventStings     = changedReader.computeEDFAnnotations()
                
                # Verifying if the annotation field is large enough to record the annotations.
                # If not, enlarge it on the writing copy.
                nbBytesEvent    = changedReader.header.nbSamplesPerRecord[self.eventChannel]*changedReader.header.nbBytes
                if max([len(s) for s in eventStings]) >= nbBytesEvent:       
                    changedReader.header.nbSamplesPerRecord[self.eventChannel] = int(max([len(s) for s in eventStings])*1.2/float(changedReader.header.nbBytes))
                    nbBytesEvent    = changedReader.header.nbSamplesPerRecord[self.eventChannel]*changedReader.header.nbBytes
               
                # Write the header
                changedReader.header.write(fileWrite)
                fileRead.seek(self.header.headerNbBytes)
        
                for noRecord in range(self.getNbRecords()):   
                    
                    record = self.readRawRecord(fileRead)                   
                        
                    for channel in changedReader.header.channelLabels: 
                        if channel == self.eventChannel:    
                            encodedString = eventStings[noRecord].encode("latin-1")  
                            writeStr = encodedString + b"\x00"*(nbBytesEvent - len(encodedString))
                            assert(len(writeStr) == nbBytesEvent)
                            fileWrite.write(writeStr)  
                        else: 
                            fileWrite.write(record[reversedDict[channel]])

        if isTemp:
            os.remove(self.fileName) 
            move(saveFileName, self.fileName)
            changedReader.fileName = self.fileName
            self = changedReader
            return self
        else:                  
            return changedReader






    def crop(self, startTime, endTime, saveFileName=None):
        
        
        """
         Note: startTime and endTime are not the exact starting and ending time.
         The starting and ending time will be those of the record containing the
         samples corresponding to startTime and endTime. We avoid to make new
         record divisions to obtain the exact starting and ending time since the
         records can be discontinuous in EDF+ format.
        """
        
        # Valid records
        noRecords = [noRecord for noRecord, rec in enumerate(self.recordsInfo) if \
                     rec.startTime <= endTime and rec.startTime + rec.duration >= startTime] 

        newReader = deepcopy(self)
        startDelta = self.recordsInfo[noRecords[0]].startTime - self.recordsInfo[0].startTime
        newReader.header.startDateTime -= datetime.timedelta(0, startDelta)
        for e in newReader.events:
            e.startTime -= startDelta
        for r in newReader.recordsInfo:
            r.startTime -= startDelta
        
        if saveFileName is None:
            newReader.save(noRecords=noRecords)
            self = newReader
        else:
            newReader.saveAs(saveFileName, noRecords=noRecords)            



    """
     Save the header information and the events to the file. The data themself
     are unchanged for two reasons:
         1 - We should always keep the original data as is an apply modifications
             (e.g., filtering) "on-line" to avoid loosing important information
             in the recorded phenomena.
         2 - The data are not loaded in kept within the object, as opposed to 
             the header and the event informations.
             
     noRecord can either be left to known or to be set to a list of record number. 
     In the latter case, only the record with numbers included in that list will
     be saved.
     
    """
    def save(self, tempPath=None, noRecords = None):

        if tempPath is None:
            tempPath = gettempdir() 
 
        # Because the annotation field may change in size, records can be shifted 
        # so we cannot only alter the information in the annotation fields. We need
        # to create a completely new file and swap it with the original.
        tempFileName = os.path.join(tempPath, "temp-" + str(uuid.uuid1()) + "-" + os.path.basename(self.fileName))

        self.saveAs(tempFileName, noRecords=noRecords)

        # New with Python 3.3, there is an os.rename() operation which is
        # atomic in fix these errors renaming or deleting files in Windows.
        #os.replace(tempFileName, self.fileName)    
        ### os.replace does not work if paths are on different drives, in Windows

        # Try to delete the original file. This may fail because the file
        # is used by another process (e.g., when running )        
        for i in range(13):
            try:
                os.remove(self.fileName)
            except OSError as e:
                message = "Failed to remove " + self.fileName + "."
                if i == 12:
                    raise WindowsError(message + str(e))
                warn(message +  "Retrying in " + str(2**i) + " seconds.", UserWarning)
                sleep(2**i)
                continue
            break            
             
        # rename the new file       
        for i in range(13):
            try:
                # Using shutil.move instead of os.rename to account for the fact
                # src and dst can be on different drives.
                move(tempFileName, self.fileName)  
            except OSError as e:
                message = "Failed to rename " + tempFileName + " to " + self.fileName  + "."
                if i == 12:
                    raise WindowsError(message + str(e))
                warn(message + " Retrying in " + str(2**i) + " seconds.", UserWarning)
                sleep(2**i)
                continue
            break          
    
        
        # Reinit the object with the new file.
        # TODO: Reinitializing the event is inefficient. The object should be
        # updated by the saving procedure such that it don't need full reinitialization.
        self.__init__(self.fileName)





    def checkAnnotationSize(self):
        # Computing the strings representing the EDF annotations
        eventStings = self.computeEDFAnnotations()       
        return [len(s) for s in eventStings]     

               


    def saveAs(self, saveFileName, noRecords=None):

        with io.open(saveFileName, 'wb') as fileWrite:
            with io.open(self.fileName, 'rb') as fileObj:
            
                # Computing the strings representing the EDF annotations
                eventStings = self.computeEDFAnnotations()       
                # Any modifications to the header for the writing must be made in a copied
                # header as it must not alter the reading of the file for the duration
                # of the saving operation.
                writeHeader = deepcopy(self.header)
                
                if not noRecords is None:
                    writeHeader.nbRecords = len(noRecords)                 
                
                
                # Verifying if the annotation field is large enough to record the annotations.
                # If not, enlarge it on the writing copy.
                if max([len(s) for s in eventStings]) >= self.header.nbSamplesPerRecord[self.eventChannel]*self.header.nbBytes:       
                    writeHeader.nbSamplesPerRecord[self.eventChannel] = int(max([len(s) for s in eventStings])/float(self.header.nbBytes)*1.2)
               
                # Write the header
                writeHeader.write(fileWrite)
                fileObj.seek(fileWrite.tell())
                assert(self.header.headerNbBytes == fileWrite.tell())

                # If it has not been specified that only certain records are to be
                # save, then specify that all records are to be saved.
                if noRecords is None: 
                    noRecords = list(range(len(eventStings)))
        
                # Write the body. Actually, we leave the data intact, updating only
                # the annotation field.
                for noRecord, eventStr in enumerate(eventStings):
                    #print noRecord, len(eventStings)
                    rawRecord = self.readRawRecord(fileObj)    
                    if noRecord in noRecords:                    
                        for channel in self.header.channelLabels:
                            if channel == self.eventChannel:    
                                encodedString = eventStr.encode("latin-1")                    
                                fileWrite.write(encodedString + b"\x00"*(writeHeader.nbSamplesPerRecord[channel]*self.header.nbBytes - len(encodedString)) )     
                            else: 
                                fileWrite.write(rawRecord[channel])
            
                
                #assert(fileWrite.tell() == writeHeader.headerNbBytes + sum(writeHeader.nbSamplesPerRecord.values())*writeHeader.nbBytes*self.getNbRecords())



 




    """ 
     Compute the EDFAnnotations strings for each records (records) from the 
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
                if not isinstance(event, Event) and str(type(event)).split('\'')[1] != 'spyndle.io.EEGDatabaseReader.Event':
                    # The second part of the condition is used because 
                    # Event.__class__.__name__ == EEGDatabaseReader.Event
                    # and is not recognised as the same class of 
                    # spyndle.io.EEGDatabaseReader.Event
                    # This could probably be arrange by making a more consistent
                    # use of relative imports.
                    raise TypeError("Type " + str(type(event)) +
                                    " received instead of type " + 
                                    str(Event) +".")        

                timeKeepingStr += event.toEDFStr()


            return timeKeepingStr


        # Computing the strings representing the EDF annotations
        eventStings = []
        nbEvents = 0
        for noRecord, recInfo in enumerate(self.recordsInfo):
            startTime = recInfo.startTime
            if len(self.events) == 0:
                eventStings.extend([prepareEventStr(startTime, [])]) 
                continue
                    
            # TODO: Optimize using the bisect module (self.events is ordered
            # accordinf to startTime).
            if noRecord == 0:                                      # first record 
                try:
                    indLim = self.events.getIndexMinStartTime(self.recordsInfo[1].startTime, inclusive=True)           
                    filteredEvents = self.events[:indLim]
                except ValueError:
                    # No events had stating time smaller than self.recordsInfo[1].startTime
                    # Thus, there is no event in this first record.
                    filteredEvents = []
                # Equivalent to but faster than
                # filteredEvents = filter(lambda e: e.startTime < self.recordsInfo[1].startTime, self.events)     
                
            elif noRecord == len(self.recordsInfo) - 1:       # last record
                try:
                    indLim = self.events.getIndexMaxStartTime(startTime, inclusive=False)           
                    filteredEvents = self.events[(indLim+1):]      
                except ValueError:
                    # No events had stating time larger than startTime
                    # Thus, there is no event in this last record.
                    filteredEvents = []                    
                # Equivalent to but faster than            
                #filteredEvents = filter(lambda e: e.startTime >= startTime, self.events)  
                
            else:                                               # other records
                try:
                    indLimMin = self.events.getIndexMinStartTime(startTime, inclusive=True)           
                    indLimMax = self.events.getIndexMaxStartTime(self.recordsInfo[noRecord+1].startTime, inclusive=False)           
                    filteredEvents = self.events[indLimMin:(indLimMax+1)]       
                except ValueError:
                    # No events had stating time larger than startTime or smaller than self.recordsInfo[noRecord+1].startTime
                    # Thus, there is no event in this last record.
                    filteredEvents = []                             
                    
                # Equivalent to but faster than                  
                #filteredEvents = filter(lambda e: e.startTime <  self.recordsInfo[noRecord+1].startTime and 
                #                                  e.startTime >= startTime, self.events)     
                
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
        for noRecord in range(self.getNbRecords()):
            indEnd = min(indStart + nbSamplesRecord, len(computeString))
            formatStings.append(computeString[indStart:indEnd]) 
            indStart = indEnd

        return formatStings









    def getChannelFreq(self, channel): 
        try:
            return self.header.samplingRate[channel]      
        except KeyError:
            raise KeyError("Channel '" + channel + "' is not available. Following channels are :" + str(self.getChannelLabels()))


    def getChannelTime(self, channel, startTime=0, timeDuration=np.inf) :
        if not isinstance(channel, str):
            raise TypeError("The channel argument must be of type str or unicode.")        
        
        #if not channel in self.getChannelLabels() :
        #    raise 
        
        endTime = startTime + timeDuration
        nbSamples = float(self.header.nbSamplesPerRecord[channel])
        recordTime = arange(nbSamples)/float(self.getChannelFreq(channel))
        recordDuration = nbSamples/float(self.getChannelFreq(channel))
        
        recs = self.recordsInfo
        
        if not startTime is 0:
            recs = [r for r in recs if startTime <=  r.startTime + recordDuration]     
            
        if not endTime is np.inf:
            recs = [r for r in recs if endTime >=  r.startTime]       

        time = concatenate([e.startTime + recordTime for e in recs])        
        return time[np.where((time >= startTime)*(time < endTime))]




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
            return encodeBDF(np.fromstring(samples, '<u1'))
            #return array([struct.unpack('<i', samples[noByte:(noByte+3)] + ('\0' if samples[noByte+2] < '\x80' else '\xff'))[0] for noByte in range(0, N, 3)])                    
             
        else:
            raise IOError






            
       
            
        #else:
        #    raise IOError("File type should be either 'EDF', 'EDF+', or 'BDF'."
        #                  + " Received value: " + str(self.header.fileType))





    """
     Digital to physical conversion.
    """
    def digital2physical(self, dig, channel):
        return (dig - self.header.digitalMin[channel]) * self.header.gain[channel] + self.header.physicalMin[channel]




    """
     Read the complete signal recorded by a given channel.
     
     Normally, decode is always True but can be set to False to get the raw
     record. This is useful for example when the signal is actually a 
     text string.
    """
    def readChannel(self, signalName, usePickled=False, startTime=None, decode=True):

        if not isinstance(signalName, str) :
            raise TypeError("The signalName argument must be a string. Received: " + str(type(signalName)))
        
        
        if usePickled:
            raise NotImplementedError
        else:
            with io.open(self.fileName, 'rb') as fileObj:
                fileObj.seek(self.header.headerNbBytes)        

                data = RecordedChannel()         
                data.samplingRate   = self.getChannelFreq(signalName)
                data.type           = "EEG"
                data.startTime      = self.header.startDateTime          
                
                N = self.getNbSample(signalName)
                data.signal = np.zeros(N)
                i = 0
                for noRecord in range(self.getNbRecords()):
                    bloc = self.readRawRecord(fileObj)[signalName]
                    if decode:
                        bloc = self.digital2physical(self.byteStr2integers(bloc), signalName)
                    else:
                        bloc = np.array([x for x in bloc])
                    data.signal[i:(i+len(bloc))] = bloc
                    i += len(bloc)
                    
                assert(i == N)
                return data    



    """
     Return the time associated with the next sample following startTime. Return
     None if there is no sample next to startTime.
    """
    def getNextSampleStartTime(self, signalName, startTime):
             
        X = self.sampleFromTime(signalName, startTime, mode="after", returnTime=True)             
        if X is None:
            return None
            
        sample, time = X
        return time

    def getEventSamples(self, event, channel=""):
        
        if event.channel == "":
            if channel == "":
                raise ValueError("The event is not associated to any channel and no"\
                             " value for the 'channel' attribute was provided. "\
                             "The notion of sample makes no sense for EDF "     \
                             "recordings if it cannot be associated with a "    \
                             "specific channel since different channel can have"\
                             " different sampling rate.")
        else:
            channel = event.channel
        startSample = self.sampleFromTime(channel, event.startTime)
        durationInSample = int(np.round(event.duration()*self.getChannelFreq(channel)))
        return startSample, durationInSample


    def getSampleFromTime(self, time, signalName, mode="closer", returnTime=False):
        return self.sampleFromTime(signalName, time, mode, returnTime)


    def sampleFromTime(self, signalName, startTime, mode="closer", returnTime=False):
        """
         Samples are discrete, time is continuous, thus the requested time is
         probably not falling exactly on a sample. The "mode" argument is used
         to set whether the returned sample should be the one right "before" time,
         right "after" time or the "closer" to time.
        """

        try:
            indRecInfo = self.recordsInfo.getIndexMaxStartTime(startTime)
        except ValueError:
            return None
        
        recInfo = self.recordsInfo[indRecInfo]
        recordStartTime = recInfo.startTime
        recordNbSamples = self.header.nbSamplesPerRecord[signalName]
        
        time = recordStartTime + arange(recordNbSamples)/float(self.getChannelFreq(signalName))
        #time = record.getStartTime() + arange(len(record.recordedSignals[signalName]))/record.samplingRates[signalName]
        IND = where(time >= startTime)[0]
        if len(IND) == 0 :
            IND = len(time)-1
        else:
            IND = IND[0]
            
        if mode == "after":
            pass
        elif mode == "before":
            IND = max(0, IND - 1)
        elif mode == "closer":
            if IND == 0:
                pass
            else:
                if startTime - time[IND - 1] < time[IND] - startTime:
                    IND -= 1
                else:
                    pass                
        else:
            raise ValueError("The 'mode' argument must be equal to 'after', 'before', or 'closer'.")

        baseSample = recordNbSamples*indRecInfo
        if returnTime:
            return baseSample + IND, time[IND]
        else:
            return baseSample + IND





    def timeFromSample(self, signalName, sample):
        
        nbSamplePerRec = self.header.nbSamplesPerRecord[signalName]
        recordNo = int(sample/nbSamplePerRec)
        dt = float(sample - nbSamplePerRec*recordNo)/float(nbSamplePerRec)*self.header.recordDuration
        return self.recordsInfo[recordNo].startTime + dt





    def read(self, signalNames, startTime, timeDuration, debug=False):
     
        records = []
        for noRecord, startTimeEvent in enumerate(self.recordsInfo):
            recordStartTime = startTimeEvent.startTime
            recordDuration  = startTimeEvent.duration            


            if recordStartTime > startTime + timeDuration:    
                break
            
            if recordStartTime + recordDuration >= startTime:
                record = self.readRecord(signalNames, noRecord+1)
                records.append(record)
            
                #try:
                #    assert(records[-1].getStartTime() == recordStartTime)
                #except AssertionError:
                #    print(records[-1].getStartTime(), recordStartTime)
                #    raise
    
            # Caution: The records need to be in order.

        assert(len(records)>0)

        info = deepcopy(records[0])                    
        signals = {}
        time    = {}                    
        for channel in signalNames:
            signals[channel] = concatenate([rec.recordedSignals[channel] for rec in records])
            time[channel]    = concatenate([rec.getStartTime() + arange(len(rec.recordedSignals[channel]))/float(self.getChannelFreq(channel))
                                                                        for rec in records])
                              
        info.time = time            
        info.recordedSignals = signals            

        if debug:
            channel = list(info.recordedSignals.keys())[0]
            time = info.getStartTime() + arange(len(info.recordedSignals[channel]))/float(info.samplingRates[channel])
            print(("Reader : ", len(records), startTime, timeDuration, time, info.samplingRates[channel]))




        returnData = {}  
        for channel in info.recordedSignals:
            time = info.time[channel]
            ind  = where((time >= startTime)*(time < startTime + timeDuration))[0]
            
            try:
                assert(len(ind)>0)
            except :
                print("time: ", time[0], " ... ", time[-1])
                print(info.getStartTime(), len(info.recordedSignals[channel]), info.samplingRates[channel])
                print(startTime, timeDuration)
                print(channel)
                raise

            returnData[channel]                = RecordedChannel()               
            returnData[channel].samplingRate   = info.samplingRates[channel]
            returnData[channel].type           = None
            returnData[channel].signal         = info.recordedSignals[channel][ind]
            returnData[channel].startTime      = time[ind[0]]      

        return returnData

        
    def setRecordDuration(self, duration):
        nbSamples = copy(self.header.nbSamplesPerRecord)
        if self.eventChannel in nbSamples:
            del nbSamples[self.eventChannel]
        
        # len(nbSamples) == 0 for edfa files.
        if len(nbSamples):
            dtMax = self.header.recordDuration/float(max(array(list(nbSamples.values()))))
            # We change the duration only if its greater than half the size of the 
            # smaller sampling period. Else, because of the digitization resolution,
            # the change has no effet. We remove the self.eventChannel from the computation
            # of the sampling period because the sampling period of this channel has
            # no meaning.
            if abs(duration - self.getRecordDuration()) > dtMax:
                #print dtMax, max(array(self.header.nbSamplesPerRecord.values())), self.header.recordDuration
                #print nbSamples.values()
                #print nbSamples        
                #print abs(duration - self.getRecordDuration()), duration, self.getRecordDuration()
                self.changeRecordDuration(duration)
        
        
    def getRecordDuration(self):
        return self.header.recordDuration
                
        
    def readRecord(self, channelList, recordId):
        # RecordId are numbered starting from 1, not from 0.

        if isinstance(channelList, str):
            channelList = [channelList]

        # Position the file cursor to the begin of the record no. recordId :
        recordPosition = self.header.headerNbBytes  + (recordId-1)*self.header.recordSize
        with io.open(self.fileName, 'rb') as fileObj:
            fileObj.seek(recordPosition)
            record = self.readFormatedRecord(fileObj)
            
        sigStart = record[0]        
        
        recordedSignals = {}
        samplingRates   = {}
        for channel in channelList:
            samplingRates[channel]    = self.getChannelFreq(channel)
            recordedSignals [channel] = record[1][channel]
                
        return EEGRecord(samplingRates, sigStart, recordedSignals, self.header.startDateTime)
        
        
        
        
    def getNbSample(self, channel=None):

        if channel is None:
            return max(self.header.nbSamplesPerRecord.values())*self.header.nbRecords
        
        if not isinstance(channel, str):
            raise TypeError("The channel parameter must be a string. (value: type) received : (" 
                                        + str(channel) + " : " + str(type(channel)) + ")")
                    
        assert(channel in self.getChannelLabels())
                                
        return self.header.nbSamplesPerRecord[channel]*self.header.nbRecords
               
        


      
      
      
    '''
     Read a record and return a dictionnary indexed by the signal labels
     containing arrays with raw bytes.  
    '''        
    def readRawRecord(self, fileObj):

        result = {}
        
        init_pos = fileObj.tell()
        for channel in self.header.channelLabels:
            nbBytes = self.header.nbSamplesPerRecord[channel]*self.header.nbBytes
            result[channel] = fileObj.read(nbBytes)
            if len(result[channel]) != nbBytes :
                raise EOFError  
                
        new_pos = fileObj.tell()                
        try:
            assert(new_pos - init_pos == self.header.recordSize)
        except:
            print(new_pos, init_pos, self.header.recordSize, new_pos - init_pos, 
                  self.header.channelLabels)
            raise
        
        return result


    '''
     Read a raw record, convert it to a (time, signals, events) tuple based on
     information in the header, and return it.
    '''    
    def readFormatedRecord(self, fileObj):

        record_start_pos = fileObj.tell()
        rawRecord = self.readRawRecord(fileObj)    
        
        #dig_min, phys_min, gain = self.dig_min, self.phys_min, self.gain
        #offset_seconds = self.currentRecordInd*self.header.record_length      
        #time = self.header.date_time + datetime.timedelta(0,offset_seconds)
        signals = {}
        events = []
        time = None

        for channel in rawRecord:
            if channel == self.eventChannel:
                ann = tal(rawRecord[channel])
                time = self.header.startDateTime + datetime.timedelta(0,ann[0][0])
                events.extend(ann[1:])
            elif channel == REFORMAT_CHANNEL:
                signals[channel] = rawRecord[channel]
            else:
                dig = self.byteStr2integers(rawRecord[channel])                                
                signals[channel] = self.digital2physical(dig, channel)
        
        if time is None:
            recordNo = (record_start_pos-self.header.headerNbBytes)/float(self.header.recordSize)
            time = self.header.startDateTime + datetime.timedelta(0, recordNo*self.header.recordDuration ) 
        
        return time, signals, events



    def changeRecordDuration(self, newDuration):
        #TODO: Implement.
        raise NotImplemented



    def checkFormat(self):
        
        def checkOverlapping():

            endTime   = 0
            startTime = -1
            for noRecord in range(self.getNbRecords()):   
                
                recordPosition = self.header.headerNbBytes  + noRecord*self.header.recordSize
                with io.open(self.fileName, 'rb') as fileObj:
                    fileObj.seek(recordPosition)
                
                    rawRecord = self.readRawRecord(fileObj)    

                    for channel in rawRecord:
                        if channel == self.eventChannel:
                            ann = tal(rawRecord[channel])
                            if endTime - ann[0][0] > 0.00000001 :
                                print(("Record " + str(noRecord) + ": Overlapping error. Record " \
                                + str(noRecord-1) + " ends at " + str(endTime) + " while record " \
                                + str(noRecord) + " starts at " + str(ann[0][0]) + "."))
                            
                            if ann[0][0] <= startTime :
                                print(("Record " + str(noRecord) + ": Record order error. Record " \
                                + str(noRecord-1) + " starts at " + str(startTime) + " while record " \
                                + str(noRecord) + " starts at " + str(ann[0][0]) + "."))                                
                                
                            
                            startTime = ann[0][0]
                            endTime   = ann[0][0] + self.header.recordDuration
                            

        checkOverlapping()








class EDFMultiReader(EDFBaseReader) :
    """
     This class allows reading multiple files using an interface that is mostly
     compatible with EDFBaseReader. 
    """

    def __init__(self, fileNames, readEventsOnInit=True, eventChannel = EVENT_CHANNEL): 
        self.readers = {f:EDFBaseReader(f, readEventsOnInit=readEventsOnInit, eventChannel=eventChannel) 
		                for f in fileNames}
        self.initEvents([reader.events for reader in list(self.readers.values())])

    def getFileName(self): 
        return [reader.getFileName() for reader in list(self.readers.values())]
        
    def getRecordingStartTime(self): 
        """
         If the different files have different start times, it returns a dict
         with file names as key an start times as values. 
         
         If all files have the same start time, just return that start time value.
        """
        startTimes = {fileName:reader.getRecordingStartTime() for fileName, reader in list(self.readers.items())}
        if len(np.unique(list(startTimes.values()))) == 1:
            return list(startTimes.values())[0]
        else:
            return startTimes

        
    def getNbRecords(self):    
        nbRecords = {fileName:reader.header.nbRecords for fileName, reader in list(self.readers.items())}
        if len(np.unique(list(nbRecords.values()))) == 1:
            return list(nbRecords.values())[0]
        else:
            return nbRecords


                
    def getChannelLabels(self):
        channelLabels = [reader.getChannelLabels() for reader in list(self.readers.values())]       
        
        # Flattening the list
        return [item for sublist in channelLabels for item in sublist]
        

    def getChannelFreq(self, channel): 
        if not isinstance(channel, str):
            raise TypeError        

        for reader in list(self.readers.values()):
            if channel in reader.getChannelLabels():
                return reader.getChannelFreq(channel)
                

    def getChannelTime(self, channel, *args, **kwargs) :
        if not isinstance(channel, str):
            raise TypeError        

        for reader in list(self.readers.values()):
            if channel in reader.getChannelLabels():
                return reader.getChannelTime(channel, *args, **kwargs)


    """
     Read the complete signal recorded by a given channel.
    """
    def readChannel(self, channel, usePickled=False):

        if not isinstance(channel, str) :
            print(("Wrong type:", type(channel), channel))
            raise TypeError        

        for reader in list(self.readers.values()):
            if channel in reader.getChannelLabels():
                return reader.channel(channel, usePickled=usePickled)



    def read(self, signalNames, startTime, timeDuration, debug=False):

        data = {}
        for reader in list(self.readers.values()):
            inReaderBool = np.in1d(signalNames, reader.getChannelLabels()) 
            inReaderChan = np.array(signalNames)[inReaderBool]
            if len(inReaderChan) > 0 :
                newData = reader.read(inReaderChan, startTime, timeDuration, debug=debug)
                data = dict(list(data.items()) + list(newData.items()))
        return data                




    def getRecordDuration(self):
        recordDuration = {fileName:reader.header.recordDuration for fileName, reader in list(self.readers.items())}
        if len(np.unique(list(recordDuration.values()))) == 1:
            return list(recordDuration.values())[0]
        else:
            return recordDuration                
                
        
    def getNbSample(self, channel=None):

        if channel is None:
            nbSample = {fileName:reader.getNbSample() for fileName, reader in list(self.readers.items())}
            if len(np.unique(list(nbSample.values()))) == 1:
                return list(nbSample.values())[0]
            else:
                return nbSample                
                    
        if not isinstance(channel, str):
            raise TypeError("The channel parameter must be a string. (value: type) received : (" 
                                        + str(channel) + " : " + str(type(channel)) + ")")

        for reader in list(self.readers.values()):
            if channel in reader.getChannelLabels():
                return reader.getNbSample(channel)
                
        raise ValueError
               






    ###########################################################################
    # Overridding some events management methods
    ##########################################################################

    def initEvents(self, eventLists):
        self._events = deepcopy(eventLists[0])
        for l in eventLists[1:]:
            for e in l:
                self._events.add(e)            
    
    @property
    def events(self, startTime=None, endTime=None):
        return self.getEvents(startTime, endTime)    
    
    @events.setter
    def events(self, events):
        raise NotImplementedError("The property events of a EDFMultiReader has "
                                  + "no setter because its content must be "
                                  + "derived from managed readers. Please use "
                                  + "initEvents(...) instead.")
        
                            
    def addEvent(self, event, fileName):
        self.readers[fileName].events.add(event)

    def removeEvent(self, event):
        for reader in list(self.readers.values()):        
            reader.events.remove(event)
        
    def removeEventType(self, eventType):
        for reader in list(self.readers.values()):        
            reader.events.removeType(eventType)        











   
   
   
   
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
    def __init__(self, fileName, isSplitted=None, annotationFileName=None,
                 readEventsOnInit=True, eventChannel = EVENT_CHANNEL):        
        
        self.eventChannel = eventChannel        
        
        if isSplitted is False:
            self.isSplitted = False
            self.annotationFileName = None
            
        elif isSplitted is True:         
            self.isSplitted = True            
            if not annotationFileName is None:
                self.annotationFileName = annotationFileName
            elif os.path.isfile(fileName + "a") :
                self.annotationFileName = fileName + "a"
            else:
                raise ValueError("In EDFReader(...), the parameter isSplitted is set to True but"\
                      " no annotation file name is specified using the annotationFileName"\
                      " parameter and no file " + fileName + "a is existing.")
                      
        elif isSplitted is None:
            if not annotationFileName is None:
                self.annotationFileName = annotationFileName
                self.isSplitted         = True 
            elif os.path.isfile(fileName + "a") :         
                self.annotationFileName = fileName + "a"
                self.isSplitted         = True 
            else:
                self.isSplitted         = False  
                self.annotationFileName = None
    
        else:
            raise ValueError("In EDFReader(...), the parameter isSplitted can only take the "\
                  "values True, False or None. Value " + str(isSplitted) + " used.")


        self.dataReader        = EDFBaseReader(fileName, readEventsOnInit=readEventsOnInit, eventChannel=eventChannel)
        super(EDFReader, self).__init__(self.getRecordDuration()) 

        if self.isSplitted :
            if isinstance(self.annotationFileName, list):
                if len(self.annotationFileName) == 1:
                    self.annotationFileName = self.annotationFileName[0]
            elif not isinstance(self.annotationFileName, str) :
                raise TypeError

            if isinstance(self.annotationFileName, list):
                self.annotationsReader  = EDFMultiReader(self.annotationFileName, readEventsOnInit=readEventsOnInit,
                                                         eventChannel=eventChannel)                
            else:
                self.annotationsReader  = EDFBaseReader(self.annotationFileName, readEventsOnInit=readEventsOnInit,
                                                        eventChannel=eventChannel)
            #self.events             = self.annotationsReader.events
        else:
            self.annotationsReader = None
            #self.events            = self.dataReader.events

    @property
    def events(self):
        if self.isSplitted:
            return self.annotationsReader.events
        else:
            return self.dataReader.events

    @events.setter
    def events(self, events):
        if self.isSplitted:
            self.annotationsReader.initEvents(events)
        else:
            self.dataReader.events = events


    def getDiscontinuityInfo(self):        
        return self.dataReader.getDiscontinuityInfo()     


    @property
    def headerString(self):
        return self.dataReader.header.headerString
                
        
    def anonymize(self):        
        if self.isSplitted :
            self.annotationsReader.header.anonymize()    
        self.dataReader.header.anonymize()        
        
    def changeChannelSet(self, chanDict, fileName=""):
        """
         This method can be used to rename channels and to keep only a 
         subset of the available channels. The argument is dictionnary of 
         channels where the key is the current channel name and the value
         is the new channel name. Channel
         labels that are not in the chanDict will be removed.
        """

        if fileName == "":
            """
            No file name has been specified. We therefore overwrite the current
            name file and change the current reader
            """
            changedReader = self
        else:
            changedReader = deepcopy(self)

        if self.isSplitted:
            annotationFileName = "" if fileName == "" else fileName+"a"
            changedReader.annotationsReader = \
                self.annotationsReader.changeChannelSet(chanDict,
                                                        annotationFileName)

        changedReader.dataReader = self.dataReader.changeChannelSet(chanDict,
                                                                    fileName)
        return changedReader
            
            
        
        
        
    def getEventStrings(self):
        """
        This methods is mostly for debugging. It reads and return the plain
        string representation of the event description.
        """        
        
        if self.isSplitted :
            return self.annotationsReader.getEventStrings() 
        else:
            return self.dataReader.getEventStrings() 
        
        


    def checkAnnotationSize(self):        
             
        if self.isSplitted :
            return self.annotationsReader.checkAnnotationSize() 
        else:
            return self.dataReader.checkAnnotationSize() 
           
        
        
        

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
        
    def getNbRecords(self):    
        return self.header.nbRecords          
                
    def getChannelLabels(self):
        return self.dataReader.getChannelLabels()

    def addEvent(self, event):
        self.events.add(event)   


    def reformatMontage(self, channelExpressions, saveFileName = "", annotationSaveFileName="", 
                        transducerTypes = None, prefilterings = None, units = None):

        if saveFileName == "":
            saveFileName = self.fileName[:-4] + "_reformatted" + self.fileName[-4:]


        # If the original recording was not a splited file and no file name
        # has been given for the annotation file, do not split the file.
        if annotationSaveFileName == "" and self.annotationFileName is None:
            pass

        # If the original recording was a splited file but no file name
        # has been given for the annotation file, use the standard naming 
        # for annotation file (exact same name as the data file but with
        # a "a" suffix to the file extension).        
        elif annotationSaveFileName == "" and not self.annotationFileName is None:
            annotationSaveFileName = saveFileName + "a"
            copyfile(self.annotationFileName, annotationSaveFileName)


        # If the original recording was not a splited file but a file name
        # has been given for the annotation file, split the reformated record.   
        elif annotationSaveFileName != "" and self.annotationFileName is None:
            self.splitRecord(annotationSaveFileName)
            self.save()

        # If the original recording was  a splited file and a file name
        # has been given for the annotation file, save the annotation file using
        # the requested name.           
        else: # anotationSaveFileName != "" and not self.annotationFileName is None

                    
        
            if isinstance(annotationSaveFileName, list):
                if len(annotationSaveFileName) == 1:
                    annotationSaveFileName = annotationSaveFileName[0]
            if isinstance(self.annotationFileName, list):
                if len(self.annotationFileName) == 1:
                    self.annotationFileName = self.annotationFileName[0]
                    
        
            if isinstance(annotationSaveFileName, list) and isinstance(self.annotationFileName, list):
                for f1, f2 in zip(self.annotationFileName, annotationSaveFileName):
                    if f1 != f2:
                        copyfile(f1, f2)                    
            elif isinstance(annotationSaveFileName, str) and isinstance(self.annotationFileName, str) :
                if self.annotationFileName != annotationSaveFileName:
                    copyfile(self.annotationFileName, annotationSaveFileName)
            else:
                raise TypeError

        self.dataReader.reformatMontage(channelExpressions, saveFileName, 
                        transducerTypes, prefilterings, units)        
                 
         
         
    def splitRecord(self, annotationFileName):

        self.isSplitted = True
        self.annotationFileName = annotationFileName          

        self.annotationsReader  = deepcopy(self.dataReader)
        self.events             = self.annotationsReader.events

        # TODO: Events in the self.dataReader could be removed (except for the time
        # stamps of the records which are mendatory for EDF+).
        #self.dataReader ...

        self.annotationsReader.header.headerNbBytes  = 8 + 80 + 80 + 8 + 8 + 8 + 44 + 8 + 8 + 4 + \
                                              (16 + 80+ 8 + 8 + 8 + 8 + 8 + 80 + 8 + 32) 
        
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
            
  
  
  



    def crop(self, startTime, endTime, saveFileName=None, saveAnnotationFileName=None):
        
        if saveAnnotationFileName is None and not saveFileName is None:
            saveAnnotationFileName = saveFileName + "a"

        if self.isSplitted:
            self.annotationsReader.crop(startTime, endTime, saveAnnotationFileName)
        self.dataReader.crop(startTime, endTime, saveFileName)


            



    def save(self, tempPath=None):
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


    def getChannelTime(self, *args, **kwargs) :
        return self.dataReader.getChannelTime(*args, **kwargs)   
 


    """
     Read the complete signal recorded by a given channel.
    """
    def readChannel(self, signalName, usePickled=False, decode=True):
        return self.dataReader.readChannel(signalName, usePickled, decode=decode)           



    """
     Return the time associated with the next sample following startTime. Return
     None if there is no sample next to startTime.
    """
    def getNextSampleStartTime(self, signalName, startTime):
        return self.dataReader.getNextSampleStartTime(signalName, startTime)                 
             

    def getSampleFromTime(self, time, signalName, mode="closer", returnTime=False):
        return self.dataReader.getSampleFromTime(time, signalName, mode, returnTime)


    def timeFromSample(self, signalName, sample):
        return self.dataReader.timeFromSample(signalName, sample)

    def getEventSamples(self, event, channel=""):
        return self.dataReader.getEventSamples(event, channel)                    
             
             

    def read(self, signalNames, startTime, timeDuration, debug=False):
        return self.dataReader.read(signalNames, startTime, timeDuration, debug)                   
          
        
    def setRecordDuration(self, duration):
        self.dataReader.setRecordDuration(duration)        

        
    def getRecordDuration(self):
        return self.dataReader.getRecordDuration()
                
        
    def readRecord(self, channelList, recordId):
        return self.dataReader.readRecord(channelList, recordId)        
        
        
    def getNbSample(self, channel=None):
        return self.dataReader.getNbSample(channel)     


    def checkFormat(self):
        
        if self.isSplitted :
            self.annotationsReader.checkFormat()
        self.dataReader.checkFormat()    
        
    def writeSAF(self, fileName):
        
        if self.isSplitted :
            reader = self.annotationsReader
        else:
            reader = self.dataReader     
        
        with io.open(fileName, 'wb') as fileWrite:
            for eventStr in reader.computeEDFAnnotations():               
                fileWrite.write(eventStr.encode("latin-1"))  
        

