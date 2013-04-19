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



from scipy import ceil

from EEGDatabaseReader import EEGDBReaderBase


import os
import re, datetime, logging
import numpy as np

EVENT_CHANNEL = 'EDF Annotations'
log = logging.getLogger(__name__)

class EDFEndOfData: pass

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



class EDFHeader :
    def __init__(self, f, fileName):

        assert f.tell() == 0  # check file position
        assert f.read(8) == '0       '
        
        
        # recording info)
        self.local_subject_id = f.read(80).strip()
        self.local_recording_id = f.read(80).strip()
        
        
        # parse timestamp
        (day, month, year) = [int(x) for x in re.findall('(\d+)', f.read(8))]
        (hour, minute, sec)= [int(x) for x in re.findall('(\d+)', f.read(8))]
        
          
        # Guessing that the data are less than 100 years old....
        year = year + (2000 if year + 2000 < datetime.date.today().year else 1900) 
        self.date_time = datetime.datetime(year, month, day, hour, minute, sec)
        
        
        # misc
        self.header_nbytes        = int(f.read(8))
        self.subtype              = f.read(44)[:5]
        self.contiguous           = self.subtype != 'EDF+D'
        self.n_records            = int(f.read(8))
        self.record_length        = float(f.read(8))  # in seconds
        self.n_channels           = int(f.read(4))
        
        
        # read channel info
        channels                  = range(self.n_channels)
        self.label                = [f.read(16).strip() for n in channels]
        self.transducer_type      = [f.read(80).strip() for n in channels]
        self.units                = [f.read(8).strip() for n in channels]
        self.physical_min         = np.asarray([float(f.read(8)) for n in channels])
        self.physical_max         = np.asarray([float(f.read(8)) for n in channels])
        self.digital_min          = np.asarray([float(f.read(8)) for n in channels])
        self.digital_max          = np.asarray([float(f.read(8)) for n in channels])
        self.prefiltering         = [f.read(80).strip() for n in channels]
        self.n_samples_per_record = [int(f.read(8)) for n in channels]
        f.read(32 * self.n_channels)  # reserved
          
        assert f.tell() == self.header_nbytes


        # Compute variables to locate pages
        #self.headerSize   = self.header['header_nbytes']    
        self.recordSize = sum(self.n_samples_per_record)*2; # Two bytes per sample    
  

        # Patch pour corriger le fait qu'Harmonie met -1 comme valeur à 
        # self.n_records.
        # Le nombre de bytes restant dans le fichier divisé par ns/sizeof(interger) 
        # va donner la valeur de records
        if self.n_records == -1:
            fileSize = os.path.getsize(fileName)
            self.n_records =  int((fileSize-self.header_nbytes)/self.recordSize);      
    

        # calculate ranges for rescaling
        self.dig_min = self.digital_min
        self.phys_min = self.physical_min
        phys_range = self.physical_max - self.physical_min
        dig_range = self.digital_max - self.digital_min
        assert np.all(phys_range > 0)
        assert np.all(dig_range > 0)
        self.gain = phys_range / dig_range        
        






   
class EDFReader(EEGDBReaderBase) :
    def __init__(self, fname):
        EEGDBReaderBase.__init__(self)
        
        self.fileName = fname
        self.file     = open(fname, 'rb')
        
        self.read_header()    
        
    def getNbPages(self):    
        return int(ceil(self.BaseEDFReader.header['n_records']*self.BaseEDFReader.header['record_length']/self.pageDuration))          
                
    def getChannelLabels(self):
        # Le -1 prend en compte que le signal EDF Annotation est exclus du record.   
        return self.BaseEDFReader.header['label'][0:-1]
        
    def getEvents(self):
        return        
        
    def readPage(self, signalNames, pageId):

        # Position de la page pageId :
        pagePosition = self.header['header_nbytes']  + (pageId-1)*self.recordSize
        self.file.seek(pagePosition)

        # TODO: Maybe replace the ID indexing by name indexing in dictionnary.
        
        record = self.read_record()
        recordedSignals = [record[1][signalID] for signalID in signalIDs]
        sigStart = record[0]        
        samplingRates = [self.header['n_samples_per_record'][signalID] for signalID in signalIDs]
        return (samplingRates, sigStart, recordedSignals)
        
        
        
        
        

    def __del__(self):
        self.file.close()


    def read_header(self):
        self.header = EDFHeader(self.file, self.fileName)
      
  
    def read_raw_record(self):
        '''Read a record with data and return a list containing arrays with raw
        bytes.
        '''
        result = []
        
        # We label records sequentially from 1 to N...
        #self.currentRecordInd = int((self.file.tell() - self.header.header_nbytes)/self.header.recordSize)   
        
        for nsamp in self.header.n_samples_per_record:
          samples = self.file.read(nsamp * 2)
          if len(samples) != nsamp * 2:
            raise EDFEndOfData
          result.append(samples)
          
        
        return result

    
    def read_record(self):
        '''Read a raw record, convert it to a (time, signals, events) tuple based on
        information in the header, and return it.
        '''
        raw_record = self.read_raw_record()    
        
        #dig_min, phys_min, gain = self.dig_min, self.phys_min, self.gain
            
        #offset_seconds = self.currentRecordInd*self.header.record_length      
        #time = self.header.date_time + datetime.timedelta(0,offset_seconds)
        signals = []
        events = []
        for (i, samples) in enumerate(raw_record):
          if self.header.label[i] == EVENT_CHANNEL:
            ann = tal(samples)
            time = self.header.date_time + datetime.timedelta(0,ann[0][0]) 
            events.extend(ann[1:])
          else:
            # 2-byte little-endian integers
            dig = np.fromstring(samples, '<i2').astype(float)
            phys = (dig - self.header.dig_min[i]) * self.header.gain[i] + self.header.phys_min[i]
            signals.append(phys)
        
        return time, signals, events




    def signalNames2IDs(self, signalNames):
        return [self.header['label'].index(signalName) for signalName in signalNames]
     
     

    def changeRecordDuration(self, newDuration):
        #TODO: Implement.
        return



    


    def records(self):
        '''
        Record generator.
        '''
        try:
          while True:
            yield self.read_record()
        except EDFEndOfData:
          pass


        
        