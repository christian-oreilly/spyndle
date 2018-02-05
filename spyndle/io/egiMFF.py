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
        
    ''' 
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
                    reformattedReader.header.nbSamplesPerRecord[self.eventChannel] = int(max([len(s) for s in eventStings])*1.2/reformattedReader.header.nbBytes)
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
                    changedReader.header.nbSamplesPerRecord[self.eventChannel] = int(max([len(s) for s in eventStings])*1.2/changedReader.header.nbBytes)
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
                    writeHeader.nbSamplesPerRecord[self.eventChannel] = int(max([len(s) for s in eventStings])/self.header.nbBytes*1.2)
               
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







    '''

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


    '''
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


    '''

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
            data.signal = self.preloaded[signalName].signal[begsam:endsam]
        else:
            data.signal = self.wonambiReader.return_dat([self.getChannelLabels().index(signalName)], 
                                                         begsam=begsam, 
                                                         endsam=endsam)[0]
        return data    



    '''

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
        
        time = recordStartTime + arange(recordNbSamples)/self.getChannelFreq(signalName)
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
        dt = float(sample - nbSamplePerRec*recordNo)/nbSamplePerRec*self.header.recordDuration
        return self.recordsInfo[recordNo].startTime + dt




    '''
    
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

    '''

        
    def setRecordDuration(self, duration):
        nbSamples = copy(self.header.nbSamplesPerRecord)
        if self.eventChannel in nbSamples:
            del nbSamples[self.eventChannel]
        
        # len(nbSamples) == 0 for edfa files.
        if len(nbSamples):
            dtMax = self.header.recordDuration/max(array(list(nbSamples.values())))
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
        
        
        
    '''
    def getNbSample(self, channel=None):
        return self.header[4]

      


