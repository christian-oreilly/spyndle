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

import warnings
import comtypes
import comtypes.client as cc
import os
import io  # For buffered writing
import ctypes, numpy
import numpy as np
from scipy.io import savemat
from scipy import array, ceil, zeros, append, arange, concatenate
from math import floor
from datetime import datetime, timedelta
import time
from sys import stdout
from copy import deepcopy
from spyndle.io.edf import EDFHeader

from .EEGDatabaseReader import EEGDBReaderBase, Event, RecordedChannel, EEGRecord, EEGRecordInfo


channelTypes = {1:"EEG", 3:"EMG", 7:"EOG", 8:"ECG", 10:"MIC", 23:"RSP"}


# Internal record size of Harmonie. According to Gaetan Poirier, internally, 
# Harmonie stores its data with records of 64 samples. The size of a data file
# is always a multiple of this record size and the true time is stored at the
# begining of every such record.
RECORD_NB = 64

OLE_TIME_ZERO = datetime(1899, 12, 30, 0, 0, 0)
def ole2datetime(oledt):
    return OLE_TIME_ZERO + timedelta(days=float(oledt))

def datetime2ole(dateTimeObj):
    d = dateTimeObj - OLE_TIME_ZERO
    return d.days + d.seconds/float(3600*24) + d.microseconds/float(3600*24*1000000)

# Harmonie constants
SIGNALFILE_FLAGS_READONLY =         0x00000001
SIGNALFILE_FLAGS_WRITEONLY =        0x00000002
SIGNALFILE_FLAGS_TRUNCATE =         0x00000004
SIGNALFILE_FLAGS_CALIBRATE =        0x00010000
SIGNALFILE_FLAGS_CALIBRATEASVOLTS = 0x00020000
SIGNALFILE_FLAGS_BASEINPUTCALIB =   0x00040000



def dayToTime(Days):
    
    x = (Days-floor(Days))*24
    heures = int(floor(x))    
    x = (x-floor(x))*60
    minutes = int(floor(x))
    x = (x-floor(x))*60
    secondes = int(floor(x))   
    x = (x-floor(x))*1000*1000
    us = int(floor(x))    
    
    return time(heures, minutes, secondes, us)

    
def dayToSecond(Days, startDay):
    return (Days-startDay)*24.0*60.0*60.0
    
def secondToDays(seconds, startDay):
    return seconds/24.0/60.0/60.0 + startDay
    



class HarmonieEvent(Event):
    def __init__(self, no, group, item, recordingStartTime, trueFreq): #
        self.no          = no
        self.groupName   = group.GetName()
        self.channel     = item.GetChannel()
        
        self.name = item.GetName()  
        

        
        # Names for sleep stages in Harmonie are not constrained such that
        # there is no standard nqmes. Tore's lab, Julie's lab and clinical
        # lab use diffrent namings. So we must check for all possibles names.
        if self.groupName.lower() == "stage" or self.groupName.lower() == "stade":
            if  self.name.lower() == "stage1" or self.name.lower() == "stade1" \
                                               or self.name.lower() == "std1" :   
                self.name = "Sleep stage 1"
                
            elif  self.name.lower() == "stage2" or self.name.lower() == "stade2" \
                                                 or self.name.lower() == "std2":   
                self.name = "Sleep stage 2"
 
            elif  self.name.lower() == "stage3" or self.name.lower() == "stade3" \
                                                 or self.name.lower() == "std3" :   
                self.name = "Sleep stage 3"
                  
            elif  self.name.lower() == "stage4" or self.name.lower() == "stade4" \
                                                 or self.name.lower() == "std4":   
                self.name = "Sleep stage 4"
                
            elif  self.name.lower() == "rem" or self.name.lower() == "sp" \
                                              or self.name.lower() == "mor":   
                self.name = "Sleep stage R"
                
            elif  self.name.lower() == "wake" or self.name == "\xc9veil" \
                                               or self.name == "ÉV.":   
                self.name = "Sleep stage W"
                
            elif  self.name.lower() == "unstaged" or self.name.lower() == "stdnd" \
                                                  or self.name.lower() ==  "stageu":                  
                self.name = "Sleep stage ?"                   
                
            else:   
                print((self.name))             
                self.name = "Sleep stage ?"
                
                
            self.groupName = "stage"

        
        # OLE complete date-time object giving the begining time of the event.
        self.dateTime    = ole2datetime(item.GetStartTime())   
          
        # In seconds, since the begining of the recording of the EEG file.
        self.startTime   = (self.dateTime - recordingStartTime).total_seconds()  
                      
        self.startSample = item.GetStartSample()
        self.sampleLength= item.GetSampleLength()
        self.color       = group.GetColor() 
        

        # For some reason, the time duration given by Harmnie does not correspond
        # to the sampleDuration according to the true frequency. This would be
        # fine if it was for example to get accurate stages length but stage length
        # are not accurate (for example, in a same recordin we can have [ 29.99612427  
        # 29.9961319   29.99807739  29.99808502  30.00003052 30.00003815]) durations.
        # To get it more manageable, we recompute the time length from the sample length.
        #self.timeLength  = item.GetTimeLength()          
        self.timeLength  = self.sampleLength/trueFreq      
    
        self.properties  = {}
        # Keys and properties should be in unicode
        for i in range(group.GetItemPropertyCount()):
            key = group.GetItemPropertyKey(i)
            self.properties[key] = item.GetItemPropertyValue(key)
            
    

            
    def __str__(self):
        STR =  str(self.no) + " " + str(self.groupName) + " " + str(self.channel) \
                    + " " + str(self.name) + " " + str(self.startTime) + " " + str(self.timeLength)
                
        for key in self.properties:
            STR += " " + key + ":" + str(self.properties[key])                       
                
        return(STR)





class HarmonieReader(EEGDBReaderBase):
    def __init__(self, fname, resynch = True, verbose=False):
        super(HarmonieReader, self).__init__()
        
        self.verbose = verbose
 
        #c:\tests>python makepy.py -i StlSignalFile.tlb
        #StlSignalFileLib
        # {72A34744-DDD9-11D1-BB8F-00001B4E6868}, lcid=0, major=1, minor=1
        #self.ISignalFile = Dispatch('{72A34754-DDD9-11D1-BB8F-00001B4E6868}')

        # TODO : ajouter cette vérification
        #ISignalFile.IsValid(fileName)

        self.fileName = fname

        tlb_id = comtypes.GUID("{72A34744-DDD9-11D1-BB8F-00001B4E6868}")
        SignalFile_id = comtypes.GUID('{72A34754-DDD9-11D1-BB8F-00001B4E6868}')

        try:        
            cc.GetModule((tlb_id, 1, 1))
        except WindowsError:
            print("Using the HarmonieReader module requires the installation of the Harmonie API. This API is normally install together with the Harmonie software.")
            raise
        
        
        import comtypes.gen.SignalFileLib as SignalLib
        import pythoncom
        
        pythoncom.CoInitialize()
        self.ISignalFile = cc.CreateObject(SignalFile_id, None, None, SignalLib.ISignalFile)
        

        try:
            self.ISignalFile.Open(str(fname), SIGNALFILE_FLAGS_READONLY)
        except:
            print(("An error has occured when trying to open " + str(fname)))
            self.ISignalFile    = None
            self.labels         = None
            self.channelFreqs  = None
            self.nbChannels     = None
            self.recordNbSamples  = None
            self.baseRecordNbSamples = None
            self.channelType        = None
            raise IOError(fname)
            
            
        #	Lecture de l'objet 'ISignalInfo'
        self.ISignalInfo = self.ISignalFile.GetSignalInfo()
        
        #	Lecture de l'objet 'IFileInfo'
        self.IFileInfo = self.ISignalInfo.GetFileInfo()
        
        #	Lecture de l'objet 'IPatientInfo'
        self.IPatientInfo =  self.ISignalInfo.GetPatientInfo()
        
        self.patientInfo = {"firstName" : self.IPatientInfo.GetFirstName(),
                            "lastName"  : self.IPatientInfo.GetLastName(),
                            "Id1"       : self.IPatientInfo.GetId1(),
                            "Id2"       : self.IPatientInfo.GetId2(),
                            "birthDate" : self.IPatientInfo.GetBirthDate(),
                            "gender"    : self.IPatientInfo.GetSex(),
                            "address"   : self.IPatientInfo.GetAddress(),
                            "city"      : self.IPatientInfo.GetCity(),
                            "state"     : self.IPatientInfo.GetState(),
                            "zipCode"   : self.IPatientInfo.GetZipCode(),
                            "country"   : self.IPatientInfo.GetCountry(),
                            "homePhone" : self.IPatientInfo.GetHomePhone(),
                            "workPhone" : self.IPatientInfo.GetWorkPhone(),
                            "comments"  : self.IPatientInfo.GetComments()}           
        
        #	Lecture de l'objet 'IRecordingCalibration'
        self.IRecordingCalibration = self.ISignalInfo.GetRecordingCalibration()
        
        #	Lecture de l'objet 'IPhysicalMontage'
        self.IPhysicalMontage = self.ISignalInfo.GetPhysicalMontage()
        
        #	Lecture de l'objet 'IRecordingMontage'
        #	ATTENTION : Il n'y a toujours qu'un seul montage d'enregistrement, à l'index '0'
        self.IRecordingMontage = self.IPhysicalMontage.GetRecordingMontage(0)


        ###########################################
        ## Getting the recording start time
        ISignalRecord = self.ISignalFile.CreateSignalRecord(RECORD_NB)
        ISignalRecord.SetStartSample(int(0))
        self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)  
        self.recordingStartTime  = ole2datetime(ISignalRecord.GetStartTime())              
        ###########################################

        # We now avoid to use IRecordingMontage.GetBaseSampleFrequency() and always
        # use self.IRecordingMontage.GetBaseSampleFrequency() to avoid any confusion.
        self.trueBaseFreq      = self.IRecordingCalibration.GetTrueSampleFrequency()
        
        ###########################################
        # Store events     
        for i in range(self.ISignalInfo.GetEventItemCount()):    
            IEventItem = self.ISignalInfo.GetEventItem(i)
            IEventGroup = IEventItem.GetGroup()
            self.events.add(HarmonieEvent(i, IEventGroup, IEventItem, self.recordingStartTime, self.trueBaseFreq))
            
        # Patch used because data recorder with previous versions of Harmonie
        # inconsistent time stamp due to discontinuous recording. These lines
        # should not be necessary with another reader.
        if resynch:
            if self.verbose :
                print("Resynchronization of the events...")
            self.resyncEvents()                
        ###########################################

        # We truncate instead of rounding because we must avoid a too long record
        # duration that would create records overlapping when converting to EDF        
        self.baseRecordNbSamples  = int(self.trueBaseFreq*self.recordDuration)   
        
        # We adjuste the record duration to account for the fact that it may not
        # be exactly X seconds as the true base frequency is not an integer
        # but a float.
        self.recordDuration       = self.baseRecordNbSamples/self.trueBaseFreq
        
        self.nbChannels         = self.IRecordingMontage.GetChannelCount()     
        
        self.labels             = [self.IRecordingMontage.GetChannelLabel(i)  for i in range(self.nbChannels)]    


        freqCorrectFactor  = self.trueBaseFreq/self.IRecordingMontage.GetBaseSampleFrequency()       
    
        self.channelFreqs = {}
        for i in range(self.nbChannels):
            self.channelFreqs[self.labels[i]] = self.IRecordingMontage.GetChannelSampleFrequency(i)*freqCorrectFactor
            
        self.recordNbSamples = {}
        for channel in self.channelFreqs:
            self.recordNbSamples[channel] = int(self.channelFreqs[channel]*self.recordDuration)
            
        self.channelType = {}
        for i in range(self.nbChannels):
            self.channelType[self.labels[i]] = self.IRecordingMontage.GetChannelType(i)            
            
  
  
  
  
        ##################################
        ## A .sig file can be troncated such that Harmonie reports an invalid
        # number of data record but does not detect the problem. Reading can be
        # performed up to the specified number of sample but the behavior of 
        # reading function are uncertain. These situations can be detected and 
        # corrected by looking at the true file size.  
        nbSampleRec = int(RECORD_NB*sum(array(list(self.channelFreqs.values()))/self.trueBaseFreq))

        # Each sample is recorded using 16 bytes and the time/date is recorded
        # at the begining of each record using 8 bytes.
        recSize = nbSampleRec*2 + 8 
 
        self.nbSamples = int(os.path.getsize(self.fileName)/recSize*RECORD_NB)
        if self.nbSamples != int(self.ISignalFile.GetRecordCount(RECORD_NB))*RECORD_NB:
            self.nbSamples = min(self.nbSamples, int(self.ISignalFile.GetRecordCount(RECORD_NB))*RECORD_NB)
            print("Warning: The file has been truncated and Harmonie is not reporting "\
                  "the correct number of samples anymore. The know problems associated"\
                  " with this situation have been corrected but unknown problems may " \
                  "arise. Consider restauring this possibly corrupted .sig file if you"\
                  " have a backup.")
            print(("record size:", recSize))
            print(("Expected number of record:", int(self.ISignalFile.GetRecordCount(RECORD_NB))))
            print(("Computed number of record:", os.path.getsize(self.fileName)/recSize))
            #print sum(array(self.channelFreqs.values())/self.trueBaseFreq)
            #print array(self.channelFreqs.values())/self.trueBaseFreq
        ####################################
        

        self.defineRecords()





    def getFileName(self):
        return self.fileName


    """
     Deines the records and save it in the self.record attribute. This is useful
     for managing the discontinuity. Every HarmonyRecord object of self.records
     contains it start/end sample and an indication of wheter this record ends
     with a discontinuity (i.e., it is an incomplete record).
    """
    def defineRecords(self):
    
        def appendRecordInfo(startSample, endSample, isComplete):
    
            nbSamples = endSample-startSample+1
            ISignalRecord = self.ISignalFile.CreateSignalRecord(nbSamples) 
        
            ISignalRecord.SetStartSample(startSample)        
            self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)
            startDateTime = ole2datetime(ISignalRecord.GetStartTime())
            
            startTime = (startDateTime - self.recordingStartTime).total_seconds()  
                    
            duration  = float(nbSamples)/self.trueBaseFreq
                
            recInfo = EEGRecordInfo(startSample=startSample, 
                                    endSample=endSample, 
                                    startTime=startTime,   
                                    duration=duration, 
                                    isComplete=isComplete)            
            
            self.recordsInfo.add(recInfo)
            
    
        discontinuitySample = array([e.startSample for e in self.events if e.groupName == "Discontinuity" or
                                                                     e.groupName == "Recording Start"], dtype=int)    
        sampleTransitions = sorted(concatenate((discontinuitySample, array([self.nbSamples], dtype=int))))
                
        noTransition = 0 
        startSample = 0                
        # We skip the first transition if it is at sample 0.
        if sampleTransitions[0] == 0:
            sampleTransitions = sampleTransitions[1:]
            
        while(True):
            
            # The record contains a discontinuity
            if (sampleTransitions[noTransition] > startSample and
                sampleTransitions[noTransition] <= startSample + self.baseRecordNbSamples):
                   
                   appendRecordInfo(startSample, sampleTransitions[noTransition]-1, isComplete = False)
                   startSample = sampleTransitions[noTransition]
                   
                   noTransition += 1
                   
                   # We reached the end of the file
                   if noTransition == len(sampleTransitions):
                       break;
                  
            # The record contains no discontinuity
            else:
                   appendRecordInfo(startSample, startSample + self.baseRecordNbSamples-1, isComplete = True)
                   startSample += self.baseRecordNbSamples    
                   
                   if startSample == self.nbSamples:
                       break;


        # Perform some sanity check
        starts = np.array([r.startSample for r in self.recordsInfo])
        ends   = np.array([r.endSample for r in self.recordsInfo])
        assert(np.all(ends[:-1]+1 == starts[1:]))
        assert(ends[-1]+1 == self.getNbSample())

    def getChannelLabels(self):
        return self.labels 

                
    def getEventsBySample(self, startSample, endSample) :
        return [e for e in self.events if (e.startSample >= startSample and e.startSample < endSample) or 
                         (e.startSample + e.sampleLength >= startSample and e.startSample + e.sampleLength < endSample)]       
        

    def getRecordingStartTime(self):    
        return self.recordingStartTime

        
        
    def getDuration(self):    # en secondes
        return self.recordDuration #float(self.nbSamples)/self.trueBaseFreq       

    def getNbSample(self, channel=None):
        if channel is None:
            return self.nbSamples
        else:
            if not (isinstance(channel, str) or isinstance(channel, str)):
                raise TypeError 
            return int(self.nbSamples*self.channelFreqs[channel]/self.trueBaseFreq)
            
            
    def getRecordNbSample(self, channel):
        if not (isinstance(channel, str) or isinstance(channel, str)):
            raise TypeError 
        return int(self.getDuration()*self.channelFreqs[channel])
                        
        
    def getElectrodesLabels(self):
        return self.labels
        
    def readRecord(self, signalNames, recordNo):

        # Indexing records from 1 to NbRecords
        recInfo = self.recordsInfo[recordNo-1]
        ISignalRecord = self.ISignalFile.CreateSignalRecord(recInfo.getNbSamples())
        ISignalRecord.SetStartSample(recInfo.startSample) 
        self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)
        
        sigStart = ole2datetime(ISignalRecord.GetStartTime()) #(ole2datetime(ISignalRecord.GetStartTime()) - self.recordingStartDateTime).total_seconds()          

        recordedSignals = {}
        indS = 0
        # TODO: Évaluer la possibilité d'utiliser GetRecordBuffer pour augmenter l'efficacité...
        record = ISignalRecord.GetRecordData()
        for channel in self.labels :
            nbSamples = int(recInfo.getNbSamples()/self.trueBaseFreq*self.channelFreqs[channel])  
                 
            if channel in signalNames:
                # Get the number of sample associate to this channel 
                recordedSignals[channel] = array(record[indS:(indS+nbSamples)])
                
            indS += nbSamples

        # TODO: Should return a "record" object which is common to all reader.
        return EEGRecord(self.channelFreqs, sigStart, recordedSignals, self.recordingStartTime)
     
     
     
    # Should return a DateTime object indicating the time when the recording
    # has begun. 

    def getChannelFreq(self, channel): 
        if channel in self.labels:
            return self.channelFreqs[channel] 
        return 0.0
    

    def getSamplingRate(self, channel) :
        warnings.warn("This function is depreacted.", DeprecationWarning, stacklevel=2)
        return self.getChannelFreq(channel)      
        
     
     
    def getSamplingRates(self):
        warnings.warn("This function is depreacted.", DeprecationWarning, stacklevel=2)
        return self.channelFreqs
     
    def getChannelTypes(self):
        return self.channelType


    def readChannel(self, signalName, usePickled=True):
        if usePickled:
            return self.readPickledChannel(signalName)
        else:
            return self.readComplete([signalName])[signalName]

    def readPickledChannel(self, signalName):

        fileName = self.fname + "_readComplete_" + signalName + ".mat"

        if os.path.exists(fileName):
            data = RecordedChannel()
            data.readPickledChannel(fileName)
            return data
        else:
            return None



    def pickleCompleteRecord(self, signalNames):

        fileName = self.fileName + "_readComplete_"

        signalNamesToPickle = []
        for signalName in signalNames:
            fileChannel = fileName + signalName + ".mat"
            if not os.path.exists(fileChannel):
                signalNamesToPickle.append(signalName)
                
        if len(signalNamesToPickle) == 0:
            return
        
        # We cannot process all channels at once for long nights at high
        # sampling rates. There is not enough available memory.
        # Array of this size correspond to approximatelly 610 M 
        # (zeros(80000000).nbytes/1024/1024).
        NbSampleMaxPerReadCompleteCall =  40000000
        NbChannelPerCall = int(NbSampleMaxPerReadCompleteCall/self.nbSamples)
        Indexes = list(range(0, len(signalNamesToPickle), NbChannelPerCall))
        if Indexes[-1] < len(signalNamesToPickle):
            Indexes.append(len(signalNamesToPickle))
        for i in range(len(Indexes)-1):   
            print(("Reading form .sig file for " + str(signalNamesToPickle[Indexes[i]:Indexes[i+1]]) + "..."))
            data = self.readComplete(signalNamesToPickle[Indexes[i]:Indexes[i+1]])

            for signalName in data:
                print(("Pickling data of " + signalName + " for next time..."))
                savemat(fileName + signalName + ".mat", data[signalName])




    # As read, but read the complete duration of the signals.
    def readComplete(self, signalNames):

        # For some reasons, the API refuse to create records with a too
        # large number of samples. We therefore have to split out entire
        # nights in sequential readings. 
        
        returnData = {}
        maxRecord = 500000
        nbPasses = int(ceil(self.getNbSample()/float(maxRecord)))
        for noPass in range(nbPasses):
            startSample =  int(maxRecord*noPass)             
            
            if noPass == nbPasses-1:
                nbSamples = self.getNbSample() - startSample
            else:
                nbSamples = maxRecord

            if noPass == 0: 
                ISignalRecord = self.ISignalFile.CreateSignalRecord(nbSamples)                
                ISignalRecord.SetStartSample(startSample)     
                self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)   
                
                startDateTime = ole2datetime(ISignalRecord.GetStartTime()) #.strftime("%a, %d %b %Y %H:%M:%S +0000")   
                for channel in signalNames:
                    returnData[channel]                = RecordedChannel()                
                    returnData[channel].signal         = zeros(self.getNbSample(channel))
                    returnData[channel].samplingRate   = self.channelFreqs[channel]
                    returnData[channel].type           = self.channelType[channel]
                    returnData[channel].startTime      = startDateTime   
            else:
                
                self.ISignalFile.InitSignalRecord(nbSamples, ISignalRecord)                
                ISignalRecord.SetStartSample(startSample)  
                self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)    
                              

            record = np.array(ISignalRecord.GetRecordData())
            #cBuffer = ISignalRecord.GetRecordBuffer()        
            #record = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(ctypes.addressof(cBuffer[0].contents), 8*cBuffer[1]), float) 
            
            indS = 0
            for channel in self.labels :

                channelNbSample = int(np.round(nbSamples/self.trueBaseFreq*self.channelFreqs[channel]))  
                
                if channel in signalNames:
                    channelStartSample = int(startSample/self.trueBaseFreq*self.channelFreqs[channel])   
                    # Get the number of sample associate to this channel 
                    returnData[channel].signal[channelStartSample:(channelStartSample+channelNbSample)] = array(record[indS:(indS+channelNbSample)])
                        
                indS += channelNbSample
                
            try:
                assert(indS ==  len(record))
            except AssertionError:
                print((indS, len(record)))
                raise


        return returnData




    """
     The API does not seem to provide any convenient way to get the sample corresponding
     to a time (using ISignalRecord.SetStartTime(startTime)  
     and then reading the record does not work). Thus this function performe this
     work.
    """
    def getSampleFromTime(self, time):
        
        approxTime   = 0.0
        approxSample = 0
        ISignalRecord = self.ISignalFile.CreateSignalRecord(RECORD_NB) 
        
        i = 0
        for i in range(30):
            # We would normally use 0.5 instead of 0.7 in the if below, but
            # because Harmonie times are only precise to the ms, we must
            # be more lenient to accomodate the rather large imprecision 
            # cause by this.
            if abs(time - approxTime) <= 0.7/self.trueBaseFreq:
                break
            
            approxSample += int(round((time - approxTime)*self.trueBaseFreq))     
            if approxSample > self.nbSamples:
                approxSample = self.nbSamples
            elif approxSample <= 0 :
                approxSample = 0
                
            ISignalRecord.SetStartSample(approxSample)        
            self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)
            startDateTime = ole2datetime(ISignalRecord.GetStartTime())
            
            approxTime = (startDateTime - self.recordingStartTime).total_seconds()  
            
            # For convergence in cases of signal discontinuity, we must ensure
            # that the condition approxTime - time < time is met.
            if approxTime >= 2.0*time:
                approxTime = time*1.5 
            
        
        if i == 29:
            raise "Failed to get sample from time."
        else:
            return approxSample


    # La fonctionnalité ISignalRecord.SetStartTime ne semble pas pouvoir être utilisée
    # pour la lecture.
    #def readWithTime(self, signalNames, startTime, timeDuration):
    def read(self, signalNames, startTime, timeDuration):
        
        recordNbSample = int(np.round(timeDuration*self.trueBaseFreq))

        ISignalRecord = self.ISignalFile.CreateSignalRecord(recordNbSample)

        ISignalRecord.SetStartSample(self.getSampleFromTime(startTime))   
        self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)

        returnData = {}
        indS = 0
        
        cBuffer = ISignalRecord.GetRecordBuffer()        
        npBuffer = numpy.core.multiarray.int_asbuffer(ctypes.addressof(cBuffer[0].contents), 8*cBuffer[1])
        record = numpy.frombuffer(npBuffer, float)
        
        startDateTime = ole2datetime(ISignalRecord.GetStartTime()) #.strftime("%a, %d %b %Y %H:%M:%S +0000")   

        for channel in self.labels:
            nbSamples = int(np.round(self.channelFreqs[channel]*timeDuration))
            if channel in signalNames:   
                returnData[channel]                = RecordedChannel()                
                returnData[channel].signal         = array(record[indS:int(indS+nbSamples)])
                returnData[channel].samplingRate   = self.channelFreqs[channel]
                returnData[channel].type           = self.channelType[channel]
                returnData[channel].startTime      = startDateTime                   
                
            indS += nbSamples

        return returnData 






    def getChannelTime(self, channel, startTime=0.0, timeDuration=np.inf) :
        if not (isinstance(channel, str) or isinstance(channel, str)):
            print((type(channel)))
            raise TypeError        
        
        freq = self.getChannelFreq(channel)
        if startTime != 0.0:
            startTime = float(self.getSampleFromTime(startTime))/freq  
            
        # We remove  0.1/freq to be sure that when we call arange(start, end, 1.0/freq)
        # we dont get the last point at ~end because it is not numerically an exact
        # multiplier of frec but just a bit larger        
        if timeDuration != np.inf:
            timeDuration = int(np.round(timeDuration*freq))/freq - 0.1/freq
        endTime = startTime + timeDuration  

        #if not channel in self.getChannelLabels() :
        #    raise 

        discontinuityEvent = [e for e in self.events if e.groupName in ["Discontinuity", "Discontinuit\xe9",
                                                                        "Recording Start", "D\xe9marrage de l'enreg."]]                                                   
        discontinuitySample = array([e.startSample for e in discontinuityEvent])                                                          
        # These sample are in base frequency. We must convert them in channel frequency.
        discontinuitySample = discontinuitySample*self.channelFreqs[channel]/self.trueBaseFreq 
        
        sampleTransitions = concatenate((discontinuitySample, [self.getNbSample(channel)]))
        
        if len(sampleTransitions) <= 2:
            nbSamples = self.getNbSample(channel)
            freq = self.getChannelFreq(channel)
            return arange(startTime, min(nbSamples/freq, endTime+1.0/freq), 1.0/freq)
        else:
            time = []
            
            for i in range(len(sampleTransitions)-1):
                nbSamples = sampleTransitions[i+1] - sampleTransitions[i]
                start = discontinuityEvent[i].startTime
                end   = start + nbSamples/freq
                if start > endTime:
                    return time
                if end < startTime:
                    continue                                        
                if end > endTime:
                    end = endTime
                if start < startTime:
                    start = startTime                    
                time = concatenate((time,  arange(start, end, 1.0/freq)))       
            return time








    # Patch parce qu'il semble y avoir un problème avec le temps associé
    # aux fuseaux de Gaétan
    def resyncEvents(self):
        ISignalRecord = self.ISignalFile.CreateSignalRecord(RECORD_NB)
        for event in self.events:
            ISignalRecord.SetStartSample(event.startSample)
            self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)
            event.dateTime    = ole2datetime(ISignalRecord.GetStartTime())
            # In seconds, since the begining of the recording of the EEG file.
            event.startTime   = (event.dateTime - self.recordingStartTime).total_seconds()   


    # Patch parce qu'il semble y avoir un problème avec le temps associé
    # aux fuseaux de Gaétan
    #def getSampleFromTime(self, approximativeSample, searchedDatetime):
    #    ISignalRecord = self.ISignalFile.CreateSignalRecord(1)
    #    ISignalRecord.SetStartSample(approximativeSample)
    #    self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)
    #    guessedDatetime    = ole2datetime(ISignalRecord.GetStartTime())   
    #    delta = searchedDatetime - guessedDatetime
    #    deltaSec = delta.days*24.0*3600.0 + delta.seconds + delta.microseconds/1000000.0
    #    print guessedDatetime, searchedDatetime, delta, deltaSec, int(round(deltaSec*self.trueBaseFreq) + approximativeSample)
    #    return int(round(deltaSec*self.trueBaseFreq) + approximativeSample)
        
            
    def addEvent(self, eventName, startSample, sampleLength, channel, description=""):
        if(self.ISignalInfo):
            self.ISignalInfo.AddEventItem(eventName, description, None, startSample,
                                          None, sampleLength, None, channel, None, None)

                             #BSTR bstrName, BSTR bstrDesc,
                             #DATE dateStartTime, LONG lStartSample,
                             #DOUBLE dTimeLength, LONG lSampleLength,
                             #IEventGroup* pEventGroup, BSTR bstrChannel,
                             #OLE_COLOR clrColor, EVisibilityStatus eVisibility



    def getNbRecords(self):
        return len(self.recordsInfo)



    def saveAsEDF(self, fileName=None, fileType = "EDF+", isSplitted=True, 
                       annotationFileName=None, channelList=None, 
                       annotationSuffixe="Annotations", psgSuffixe="PSG",
                       basePath=None, annotationBasePath=None, 
                       verbose=None, contiguous=False, 
                       recordStart=0, recordEnd=None, useIntegers=True):
        
        """
         Convert and save the .sig/.sts files in EDF/BDF format. 
         
         When working with large datafiles such as whole night polysomnographic 
         recording, recording the annotations in the same file as the recording
         data can be inneficient since adding new annotations might requires 
         reformatting the whole data file. In these case, it can be more interesting
         to split the original data and the annotations in two separate files. 
         
         We suggest to always use the same name steam for both files, adding only
         a space separated "PSG" or "Annotations" sufixes  
         (e.g. someFile PSG.edf/someFile Annotations.edf).
         Nevertheless, this class provide the option of using an arbitrary name 
         for the annotation file.
    
         The parameter isSplitted (boolean) can be used to specify wheter the 
         data file specified by fname is a splitted data/annotation set of files or
         a single file. If isSplitted == True, the annotationFileName
         parameter will specify the name of the annotation file.
         If not specified, fname[-4] + " Annotations" + fname[-4:] will be used.
         If isSplitted == False, no supplementary annotation file is considered
         regardless of the value of annotationFileName.
         
         fileName   : File name to use to record the converted file. If a 
                      basePath is given, the method only use the file name
                      and append it to basePath. If left to None, use the name
                      of the source file, replacing its extension.
                      
         channelList: If different than None, only the channel in that list 
                      will be saved in the converted file. 
                      
         basePath   : If different than None, record the in this path. Else,
                      use the same path as the source file or as in fileName. 
                      
         useIntegers: Harmonie works with precise float sampling frequency. This
                       cause a lot of problems to many EDF reader which expect
                       integer values for record duration, sampling frequency,
                       and so on. Be setting this parameter to True, the 
                       reader will approximate the non integer sampling frequency
                       (e.g., 256.01223456) by an integer frequency (e.g., 256).
        """    
        if verbose is None:
            verbose = self.verbose

        returnOutput = {}

        if fileName is None:
            fileSteam = ".".join(self.fileName.split(".")[:-1])
            # Use the same name as the source name, chnaging the extension.
            if isSplitted:
                    fileName = fileSteam + " " + psgSuffixe + "." + fileType[:3].lower()
                    
            else:
                fileName = fileSteam + "." + fileType[:3].lower()
        else:
            fileSteam = ".".join(fileName.split(".")[:-1])
            
        if not basePath is None:
            fileName = os.path.join(basePath, os.path.basename(fileName))
        

        if isSplitted:        
            if annotationFileName is None:
                annotationFileName = fileSteam + " " + annotationSuffixe + "." + fileType[:3].lower()            
                if not basePath is None:
                    annotationFileName = os.path.join(basePath, os.path.basename(annotationFileName))
                    
            if not annotationBasePath is None:
                annotationFileName = os.path.join(annotationBasePath, os.path.basename(annotationFileName))                
        
        
        returnOutput["fileName"] = fileName        
        returnOutput["annotationFileName"] = annotationFileName        
        
        """
         By default, every records are used in the conversion, but this behavior
         can be overide by specifying recordStart and recordEnd values. 
        """
        if recordEnd is None:
            if recordStart == 0:
                completeFile = True
            else:
                completeFile = False                    
            recordEnd = self.getNbRecords()
        else:
            completeFile = False            
        
        if recordStart:
            ISignalRecord = self.ISignalFile.CreateSignalRecord(RECORD_NB)   
            ISignalRecord.SetStartSample(self.recordsInfo[recordStart].startSample)
            self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)         
            recordingStartTime = ole2datetime(ISignalRecord.GetStartTime()) 
        else:
            recordingStartTime = self.recordingStartTime

        recordsInfo = self.recordsInfo[recordStart:(recordEnd+1)]

        #print "start, end, ", recordStart, recordEnd        
        
        
        if useIntegers:
            # Duration of stage events should always be an integer value 
            # (typically 20 or 30 seconds). Harmonie will return a float duration
            # of type 19.99999826160749. To avoid problems with edf reader, this
            # value is rounded.            
            for e in self.events:
                if e.groupName == "stage" :
                    e.timeLength = np.round(e.timeLength)
                
            self.timeFactor = np.round(self.recordDuration)/self.recordDuration
        
        #######################################################################
        # Internal methods
        #######################################################################        
        
        def adjustTime(time):
            if useIntegers:
                return time*self.timeFactor
            else:
                return time
        
        def subjectFields(field):
            return field if field else "X"
            
        # This conversion is needed because timedate.strftime("%d-%b-%Y") would
        # use local month abbreviation instead of english ones.
        def dateStr(OLEdate) :
            month = {1:"JAN", 2:"FEB", 3:"MAR", 4:"APR", 5:"MAY", 6:"JUN", 
                     7:"JUL", 8:"AUG", 9:"SEP", 10:"OCT", 11:"NOV", 12:"DEC"}
            return "%02d-%s-%04d" % (OLEdate.day, month[OLEdate.month], OLEdate.year)
        

        def prepareEventStr(noRecord, events):
            
            baseTime = recordsInfo[noRecord].startTime

            if contiguous:
                if useIntegers:
                    offsetTime = noRecord*np.round(self.recordDuration)
                else:
                    offsetTime = noRecord*self.recordDuration
            else:
                offsetTime = adjustTime(baseTime - recordsInfo[0].startTime)   
                # The rounding is necessary to avoid overlaping errors because
                # of numerical imprecision.
                offsetTime = np.round(offsetTime, 2)
            
            # Annotation channel    

            timeKeepingStr = "+" + str(offsetTime) + "\x14\x14\0"                                 
                            
                            
            for event in events:
                copiedEvent  = deepcopy(event)
                copiedEvent.startTime = adjustTime(event.startTime - baseTime) + offsetTime
                #print copiedEvent.startTime , event.startTime, baseTime, offsetTime
                if copiedEvent.startTime >= 0.0:
                    timeKeepingStr += copiedEvent.toEDFStr() #edfEventEncode(event)  

            return timeKeepingStr
 
 
        #######################################################################
        # Encoding events
        #######################################################################
 
        # Encoding events such that events are writen in their respective record.
        # This has the disaventage to use a lot of space (each records have
        # is reserved enough space to record as much annotation as the 
        # record with the most annotation) but to make complete data records.
        if verbose:
            print("Encoding events...")


        limitTimes = [rec.startTime for rec in recordsInfo]
        limitTimes.append(recordsInfo[-1].startTime + recordsInfo[-1].duration)

        nbEvents = 0      
        for noRecord in range(len(recordsInfo)):                
            try:
                indMin = self.events.getIndexMinStartTime(limitTimes[noRecord])
            except ValueError:
                indMin = len(self.events)

            try:
                indMax = self.events.getIndexMaxStartTime(limitTimes[noRecord+1], inclusive=False)
            except ValueError:
                # This exception is trigered when there is no events starting 
                # before startTime such that no corresponding index can be returned.                 
                indMax = 0


            filteredEvents = self.events[indMin:(indMax+1)]          
            recordsInfo[noRecord].eventStr = prepareEventStr(noRecord, filteredEvents)                                           
            nbEvents += len(filteredEvents)  


        ### Only a part of the sig file can be saved. In these cases, only the  
        ### events in the saved part will be saved in the EDF file. This number
        ### can be smaller than the total number of event in the sig file.
        try:
            if completeFile:
                assert(nbEvents == len(self.events))  
            else:
                assert(nbEvents <= len(self.events))         
        except AssertionError:
            print((nbEvents, len(self.events))) 
            raise





        #######################################################################
        # Building the EDF Header object
        #######################################################################

        dataHeader = EDFHeader()       


        id      = subjectFields(self.patientInfo["Id1"])
        gender  = "M" if self.patientInfo["gender"] == 1 else "F"
        date    = dateStr(ole2datetime(self.patientInfo["birthDate"])) if self.patientInfo["birthDate"] else "X"     
        fname   =  subjectFields(self.patientInfo["firstName"]) #.encode("ascii", "ignore")
        lname   = subjectFields(self.patientInfo["lastName"]) #.encode("ascii", "ignore")  

        dataHeader.subjectID = id + " " + gender + " " +  date + " " + fname + "_" + lname   
        
        startdate    = dateStr(recordingStartTime) if recordingStartTime else "X"                     
        dataHeader.recordingIR = "Startdate " + startdate + " " + "X" + " " + "X" + " " + "X"  

        dataHeader.startDateTime = recordingStartTime #+ timedelta(0, recordsInfo[0].startTime)
            
        # Building the channel labels for EDF recording
        harmonieLabels = {}
        EDFLabels = []
        for channel in self.getChannelLabels():
            if channelList is None or channel in channelList:
                if self.channelType[channel] in channelTypes:
                    if channelTypes[self.channelType[channel]] in ["EEG", "EMG", "EOG", "ECG"] :
                        prefix = channelTypes[self.channelType[channel]] + " "
                    elif channelTypes[self.channelType[channel]] == "MIC":
                        prefix = "Sound "
                    elif channelTypes[self.channelType[channel]] == "RSP":
                        prefix = "Resp "
                    else:
                        raise ValueError("Unmanaged channel type.")
                else:
                    prefix = ""
                    
                harmonieLabels[prefix + channel] = channel
                EDFLabels.append(prefix + channel)

        ns = len(EDFLabels)
        dataHeader.headerNbBytes      = 8 + 80 + 80 + 8 + 8 + 8 + 44 + 8 + 8 + 4 + (ns+1)* (16 + 80+ 8 + 8 + 8 + 8 + 8 + 80 + 8 + 32) 
        dataHeader.nbRecords          = len(recordsInfo)
        
        # We need to truncate the record duration to the ms because the starting
        # of the records extracted from Harmony seems truncated to the ms. If not
        # truncating the recordDuration, this results in sub-millisecond overlapping
        # of records wich cause problems when reading the EDF with some readers.
        if useIntegers:
            dataHeader.recordDuration     = np.round(self.recordDuration)
        else:
            dataHeader.recordDuration     = self.recordDuration #float(str(self.recordDuration)[:len('%.*f' % (3, self.recordDuration))])
            
        dataHeader.nbChannels         = ns +1
        
        dataHeader.channelLabels      = EDFLabels + ["EDF Annotations"]        
        harmonieLabels["EDF Annotations"]   = "EDF Annotations"        
        
        dataHeader.transducerType     = {}
        dataHeader.prefiltering       = {} 
        dataHeader.units              = {}
        dataHeader.physicalMin        = {}
        dataHeader.physicalMax        = {}          
        dataHeader.digitalMin         = {}
        dataHeader.digitalMax         = {}         
        dataHeader.nbSamplesPerRecord = {} 


        physicalMinMicro = {}
        physicalMaxMicro = {}
        digitalMin = {}
        digitalMax = {}
        
        dataHeader.contiguous = contiguous
        dataHeader.fileType = fileType
        if fileType[:3] == "EDF":
            nbByte = 2        
        elif fileType == "BDF":
            nbByte = 3
        else:
            raise ValueError("The fileType argument has an invalid value: " + fileType)

        annotationFieldLength = int(max(400, max(array([len(record.eventStr) for record in recordsInfo]))/nbByte*1.2))             
        
        for i, channel in enumerate(EDFLabels + ["EDF Annotations"]) :   
            
            if fileType[:3] == "EDF":
                digitalMin[channel] = -32768
                digitalMax[channel] =  32767
            elif fileType == "BDF":
                digitalMin[channel] = -8388608 
                digitalMax[channel] =  8388607    
            else:
                raise ValueError

            dataHeader.transducerType[channel]  = ""
            dataHeader.prefiltering[channel]    = "" 


            if i == ns:
                dataHeader.units[channel]       = ""                 
                dataHeader.physicalMin[channel] = -1
                dataHeader.physicalMax[channel] = 1
                dataHeader.digitalMin[channel]  = -32768    
                dataHeader.digitalMax[channel]  = 32767  
                if isSplitted:
                    dataHeader.nbSamplesPerRecord[channel] = 15
                else:
                    dataHeader.nbSamplesPerRecord[channel] = annotationFieldLength
            else:
                inputMin, outputMin, inputMax, outputMax = self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATEASVOLTS)
                physicalMaxMicro[channel]  = inputMax*(outputMax-outputMin)/(inputMax - inputMin)  
                physicalMinMicro[channel]  = inputMin*(outputMax-outputMin)/(inputMax - inputMin)                
                    
                dataHeader.units[channel]               = "uV"                 
                dataHeader.physicalMin[channel]         = physicalMinMicro[channel]
                dataHeader.physicalMax[channel]         = physicalMaxMicro[channel]
                dataHeader.digitalMin[channel]          = digitalMin[channel]
                dataHeader.digitalMax[channel]          = digitalMax[channel]
                dataHeader.nbSamplesPerRecord[channel]  = self.recordNbSamples[harmonieLabels[channel]]




        #######################################################################
        # Writing the annotation file if the EEG data are to be saved in a 
        # splitted set of data and annotation files.
        #######################################################################
       
        if isSplitted: 
            annotationHeader = deepcopy(dataHeader)
            
            annotationHeader.headerNbBytes  = 8 + 80 + 80 + 8 + 8 + 8 + 44 + 8 + 8 + 4 + (16 + 80+ 8 + 8 + 8 + 8 + 8 + 80 + 8 + 32) 
            annotationHeader.nbChannels     = 1
            annotationHeader.channelLabels  = ["EDF Annotations"] 

            annotationHeader.transducerType     = {"EDF Annotations":""}
            annotationHeader.prefiltering       = {"EDF Annotations":""} 
            annotationHeader.units              = {"EDF Annotations":""}
            annotationHeader.physicalMin        = {"EDF Annotations":-1}
            annotationHeader.physicalMax        = {"EDF Annotations":1}          
            annotationHeader.digitalMin         = {"EDF Annotations":-32768}
            annotationHeader.digitalMax         = {"EDF Annotations":32767}         
            annotationHeader.nbSamplesPerRecord = {"EDF Annotations":annotationFieldLength} 
            
                      
            with io.open(annotationFileName, 'wb') as f:    
                
                if verbose:
                    print("Writing annotation file header...")
                
                annotationHeader.write(f)

                if verbose:
                    print("Writing annotation file body...")
        
                for noRecord, record in enumerate(recordsInfo):            
                    encodedRecordStr = record.eventStr.encode("latin-1")
                    f.write(encodedRecordStr +  b"\0"*(annotationFieldLength*nbByte-len(encodedRecordStr)))         

               
       
       
       
       
       
        #######################################################################
        # Writing the data file...
        #######################################################################       
       
        # Using buffered writer
        with io.open(fileName, 'wb') as f:


            if verbose:
                print("Writing data file header...")
            dataHeader.write(f)



            if verbose:
                print("Writing date file body...")
                        
            ISignalRecord = self.ISignalFile.CreateSignalRecord(RECORD_NB) 

            for noRecord, record in enumerate(recordsInfo):    
                self.ISignalFile.InitSignalRecord(record.getNbSamples(), ISignalRecord) 

                ISignalRecord.SetStartSample(record.startSample)

                # This function call is the bottle neck of this section.
                self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)

                npRecord = np.array(ISignalRecord.GetRecordData())
                #print(len(npRecord), RECORD_NB)
                #cBuffer = ISignalRecord.GetRecordBuffer()        
                #npRecord = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(ctypes.addressof(cBuffer[0].contents), 8*cBuffer[1]), float) 

                # record is in uV

                #print noRecord+1, f.tell()
                                        
                indS = 0        
                for hChannel in self.getChannelLabels():
                    try:                    
                        edfChannel = dict((v, k) for (k, v) in list(harmonieLabels.items()))[hChannel] 
                    except KeyError:
                        edfChannel = ""
                        
                    if edfChannel == 'EDF Annotations':                 
                        continue
                    
                    # If this record is discontinuous within the .sig file, only 
                    # the part up to the discontinuity is used: the rest of 
                    # this EDF signal is padded with zeros to make it the correct
                    # size of a complete record.                    
                    if not record.isComplete:
                        # Get the number of sample associate to the incomplete records for that channel 

                        channelRecordNbSample = int(record.getNbSamples()/self.trueBaseFreq*self.channelFreqs[hChannel])  
                        
                        recordedSignal = append(array(npRecord[indS:(indS+channelRecordNbSample)]), 
                                                zeros(self.recordNbSamples[hChannel] - channelRecordNbSample))
                        indS += channelRecordNbSample
                        
                    else:
                        recordedSignal = array(npRecord[indS:(indS+self.recordNbSamples[hChannel])])
                        indS += self.recordNbSamples[hChannel]
                    

                    if edfChannel in EDFLabels:
                                       
                        # WRITE RECORDED SIGNAL....
                        physical_min = physicalMinMicro[edfChannel] 
                        physical_max = physicalMaxMicro[edfChannel] 
                        digital_min  = digitalMin[edfChannel] 
                        digital_max  = digitalMax[edfChannel] 
    
                        phys_range = physical_max - physical_min
                        dig_range = digital_max - digital_min
                        assert numpy.all(phys_range > 0)
                        assert numpy.all(dig_range > 0)
                        gain = dig_range/phys_range                          
    
                        recordedSignal = (recordedSignal - physical_min)*gain + digital_min  
                        
                        if nbByte == 2: # EDF
                            f.write(recordedSignal.astype('<i2').tostring())
                        elif nbByte == 3: #BDF
                            # Writing in a string of 32-bit integers and removing every fourth byte
                            # to obtain a string of 24-bit integers
                            recordedSignal = recordedSignal.astype('<i').tostring()
                            recordedSignal = bytes([recordedSignal[noBit] for noBit in range(len(recordedSignal)) if (noBit+1)%4])
                            # print(len(recordedSignal), self.recordNbSamples[hChannel], len(recordedSignal))
                            f.write(recordedSignal)
                                                    
               
                
                # Annotation channel            
                if isSplitted:
                     eventStr = record.eventStr.split("\x14\x14\0")[0] + "\x14\x14\0"
                else:
                     eventStr = record.eventStr
                     
                encodedRecordStr = eventStr.encode("latin-1")
    
                f.write(encodedRecordStr +  b"\x00"*(dataHeader.nbSamplesPerRecord["EDF Annotations"]*nbByte-len(encodedRecordStr)))         

                if verbose:
                    done=float(noRecord)/len(recordsInfo)*100.0
                    stdout.write(" Body writing percentage: %s%%      %s"%(done,"\r"))
                    stdout.flush()
        
        
        return returnOutput

