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

from EEGDatabaseReader import EEGDBReaderBase, Event, RecordedChannel, EEGPage, EEGPageInfo


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
    def __init__(self, no, group, item, recordingStartTime):
        self.no          = no
        self.groupName   = group.GetName()
        self.channel     = item.GetChannel()
        
        self.name = item.GetName()  
        

        
        # Names for sleep stages in Harmonie are not constrained such that
        # there is no standard nqmes. Tore's lab, Julie's lab and clinical
        # lab use diffrent namings. So we must check for all possibles names.
        if self.groupName.lower() == u"stage" or self.groupName.lower() == u"stade":
            if  self.name.lower() == u"stage1" or self.name.lower() == u"stade1" \
                                               or self.name.lower() == u"std1" :   
                self.name = u"Sleep stage 1"
                
            elif  self.name.lower() == u"stage2" or self.name.lower() == u"stade2" \
                                                 or self.name.lower() == u"std2":   
                self.name = u"Sleep stage 2"
 
            elif  self.name.lower() == u"stage3" or self.name.lower() == u"stade3" \
                                                 or self.name.lower() == u"std3" :   
                self.name = u"Sleep stage 3"
                  
            elif  self.name.lower() == u"stage4" or self.name.lower() == u"stade4" \
                                                 or self.name.lower() == u"std4":   
                self.name = u"Sleep stage 4"
                
            elif  self.name.lower() == u"rem" or self.name.lower() == u"sp" \
                                              or self.name.lower() == u"mor":   
                self.name = u"Sleep stage R"
                
            elif  self.name.lower() == u"wake" or self.name == u"\xc9veil" \
                                               or self.name == u"ÉV.":   
                self.name = u"Sleep stage W"
                
            elif  self.name.lower() == u"unstaged" or self.name.lower() == u"stdnd" \
                                                  or self.name.lower() ==  "stageu":                  
                self.name = u"Sleep stage ?"                   
                
            else:   
                print self.name             
                self.name = u"Sleep stage ?"
                
                
            self.groupName = u"stage"
                                                         

        
        # OLE complete date-time object giving the begining time of the event.
        self.dateTime    = ole2datetime(item.GetStartTime())   
          
        # In seconds, since the begining of the recording of the EEG file.
        self.startTime   = (self.dateTime - recordingStartTime).total_seconds()  
        
        self.timeLength  = item.GetTimeLength()                
        self.startSample = item.GetStartSample()
        self.sampleLength= item.GetSampleLength()
        self.color       = group.GetColor() 

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
    def __init__(self, fname):
        super(HarmonieReader, self).__init__()
 
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
            print "Using the HarmonieReader module requires the installation of the Harmonie API. This API is normally install together with the Harmonie software."
            raise
        
        
        import comtypes.gen.SignalFileLib as SignalLib
        import pythoncom
        
        pythoncom.CoInitialize()
        self.ISignalFile = cc.CreateObject(SignalFile_id, None, None, SignalLib.ISignalFile)
        

        try:
            self.ISignalFile.Open(str(fname), SIGNALFILE_FLAGS_READONLY)
        except:
            print "An error has occured when trying to open " + str(fname)
            self.ISignalFile    = None
            self.labels         = None
            self.channelFreqs  = None
            self.nbChannels     = None
            self.pageNbSamples  = None
            self.basePageNbSamples = None
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

        
        ###########################################
        # Store events     
        for i in range(self.ISignalInfo.GetEventItemCount()):    
            IEventItem = self.ISignalInfo.GetEventItem(i)
            IEventGroup = IEventItem.GetGroup()
            self.events.add(HarmonieEvent(i, IEventGroup, IEventItem, self.recordingStartTime))
            
        # Patch used because data recorder with previous versions of Harmonie
        # inconsistent time stamp due to discontinuous recording. These lines
        # should not be necessary with another reader.
        print "Resynchronization of the events..."
        self.resyncEvents()                
        ###########################################
            
            
        self.pageDuration = np.median([e.timeLength for e in self.events if e.groupName.lower() == "stage"])

        # We now avoid to use IRecordingMontage.GetBaseSampleFrequency() and always
        # use self.IRecordingMontage.GetBaseSampleFrequency() to avoid any confusion.
        self.trueBaseFreq      = self.IRecordingCalibration.GetTrueSampleFrequency()        
        self.basePageNbSamples  = int(round(self.trueBaseFreq*self.pageDuration))   
        
        # We adjuste the page duration to account for the fact that it may not
        # be exactly X seconds as the true base frequency is not an integer
        # but a float.
        self.pageDuration       = self.basePageNbSamples/self.trueBaseFreq
        
        self.nbChannels         = self.IRecordingMontage.GetChannelCount()     
        
        self.labels             = [self.IRecordingMontage.GetChannelLabel(i)  for i in range(self.nbChannels)]    


        freqCorrectFactor  = self.trueBaseFreq/self.IRecordingMontage.GetBaseSampleFrequency()       
    
        self.channelFreqs = {}
        for i in range(self.nbChannels):
            self.channelFreqs[self.labels[i]] = self.IRecordingMontage.GetChannelSampleFrequency(i)*freqCorrectFactor
            
        self.pageNbSamples = {}
        for channel in self.channelFreqs:
            self.pageNbSamples[channel] = int(self.channelFreqs[channel]*self.pageDuration)
            
        self.channelType = {}
        for i in range(self.nbChannels):
            self.channelType[self.labels[i]] = self.IRecordingMontage.GetChannelType(i)            
            
  
  
  
  
        ##################################
        ## A .sig file can be troncated such that Harmonie reports an invalid
        # number of data record but does not detect the problem. Reading can be
        # performed up to the specified number of sample but the behavior of 
        # reading function are uncertain. These situations can be detected and 
        # corrected by looking at the true file size.  
        nbSampleRec = int(RECORD_NB*sum(array(self.channelFreqs.values())/self.trueBaseFreq))

        # Each sample is recorded using 16 bytes and the time/date is recorded
        # at the begining of each record using 8 bytes.
        recSize = nbSampleRec*2 + 8 
 
        self.nbSamples = os.path.getsize(self.fileName)/recSize*RECORD_NB
        if self.nbSamples != int(self.ISignalFile.GetRecordCount(RECORD_NB))*RECORD_NB:
            self.nbSamples = min(self.nbSamples, int(self.ISignalFile.GetRecordCount(RECORD_NB))*RECORD_NB)
            print "Warning: The file has been truncated and Harmonie is not reporting "\
                  "the correct number of samples anymore. The know problems associated"\
                  " with this situation have been corrected but unknown problems may " \
                  "arise. Consider restauring this possibly corrupted .sig file if you"\
                  " have a backup."
            print "record size:", recSize
            print "Expected number of record:", int(self.ISignalFile.GetRecordCount(RECORD_NB))
            print "Computed number of record:", os.path.getsize(self.fileName)/recSize
            #print sum(array(self.channelFreqs.values())/self.trueBaseFreq)
            #print array(self.channelFreqs.values())/self.trueBaseFreq
        ####################################
        

        self.definePages()





    def getFileName(self):
        return self.fileName


    """
     Deines the pages and save it in the self.page attribute. This is useful
     for managing the discontinuity. Every HarmonyPage object of self.pages
     contains it start/end sample and an indication of wheter this page ends
     with a discontinuity (i.e., it is an incomplete page).
    """
    def definePages(self):
    
        discontinuitySample = array([e.startSample for e in self.events if e.groupName == "Discontinuity" or
                                                                     e.groupName == "Recording Start"], dtype=int)    
        sampleTransitions = sorted(concatenate((discontinuitySample, array([self.nbSamples+1], dtype=int))))
                
        noTransition = 0 
        startSample = 0                
        # We skip the first transition if it is at sample 0.
        if sampleTransitions[0] == 0:
            sampleTransitions = sampleTransitions[1:]
            
        while(True):
            
            # The page contains a discontinuity
            if (sampleTransitions[noTransition] > startSample and
                sampleTransitions[noTransition] <= startSample + self.basePageNbSamples):
                   
                   self.getInfoPages().append(EEGPageInfo(startSample, sampleTransitions[noTransition]-1, False))
                   startSample = sampleTransitions[noTransition]
                   
                   noTransition += 1
                   
                   # We reached the end of the file
                   if noTransition == len(sampleTransitions):
                       break;
                  
            # The page contains no discontinuity
            else:
                   self.getInfoPages().append(EEGPageInfo(startSample, startSample + self.basePageNbSamples, True))
                   startSample += self.basePageNbSamples    
                   
                   if startSample == self.nbSamples:
                       break;


    def getChannelLabels(self):
        return self.labels 

                
    def getEventsBySample(self, startSample, endSample) :
        return filter(lambda e: (e.startSample >= startSample and e.startSample < endSample) or 
                         (e.startSample + e.sampleLength >= startSample and e.startSample + e.sampleLength < endSample) , self.events)       
        

    def getRecordingStartTime(self):    
        return self.recordingStartTime

        
        
    def getDuration(self):    # en secondes
        return float(self.nbSamples)/self.trueBaseFreq       

    def getNbSample(self, channel=None):
        if channel is None:
            return self.nbSamples
        else:
            if not (isinstance(channel, str) or isinstance(channel, unicode)):
                raise TypeError 
            return int(self.getDuration()*self.channelFreqs[channel])
            
        
    def getElectrodesLabels(self):
        return self.labels
        
    def readPage(self, signalNames, pageNo):

        # Indexing pages from 1 to NbPages
        ISignalRecord = self.ISignalFile.CreateSignalRecord(self.getInfoPages(pageNo).getNbSamples())
        ISignalRecord.SetStartSample(self.getInfoPages(pageNo).startSample) 
        self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)
        
        sigStart = ole2datetime(ISignalRecord.GetStartTime()) #(ole2datetime(ISignalRecord.GetStartTime()) - self.recordingStartDateTime).total_seconds()          

        recordedSignals = {}
        indS = 0
        # TODO: Évaluer la possibilité d'utiliser GetRecordBuffer pour augmenter l'efficacité...
        record = ISignalRecord.GetRecordData()
        for channel in self.labels :
            nbSamples = int(self.getInfoPages(pageNo).getNbSamples()/self.trueBaseFreq*self.channelFreqs[channel])  
                 
            if channel in signalNames:
                # Get the number of sample associate to this channel 
                recordedSignals[channel] = array(record[indS:(indS+nbSamples)])
                
            indS += nbSamples

        # TODO: Should return a "page" object which is common to all reader.
        return EEGPage(self.channelFreqs, sigStart, recordedSignals, self.recordingStartTime)
     
     
     
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

        fileName = self.fname + "_readComplete_"

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
        Indexes = range(0, len(signalNamesToPickle), NbChannelPerCall)
        if Indexes[-1] < len(signalNamesToPickle):
            Indexes.append(len(signalNamesToPickle))
        for i in range(len(Indexes)-1):   
            print "Reading form .sig file for " + str(signalNamesToPickle[Indexes[i]:Indexes[i+1]]) + "..."
            data = self.readComplete(signalNamesToPickle[Indexes[i]:Indexes[i+1]])

            for signalName in data:
                print "Pickling data of " + signalName + " for next time..."
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
                              

            cBuffer = ISignalRecord.GetRecordBuffer()        
            record = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(ctypes.addressof(cBuffer[0].contents), 8*cBuffer[1]), float) 
            
            indS = 0
            for channel in self.labels :
                channelNbSample = int(nbSamples/self.trueBaseFreq*self.channelFreqs[channel])  
                
                if channel in signalNames:
                    channelStartSample = int(startSample/self.trueBaseFreq*self.channelFreqs[channel])   
                    # Get the number of sample associate to this channel 
                    returnData[channel].signal[channelStartSample:(channelStartSample+channelNbSample)] = array(record[indS:(indS+channelNbSample)])
                        
                indS += channelNbSample
                
            assert(indS ==  len(record))



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
            if abs(time - approxTime) <= 0.5/self.trueBaseFreq:
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
                
        recordNbSample = timeDuration*self.trueBaseFreq
        ISignalRecord = self.ISignalFile.CreateSignalRecord(int(recordNbSample))
        ISignalRecord.SetStartSample(self.getSampleFromTime(startTime))   
        self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)

        returnData = {}
        indS = 0
        
        cBuffer = ISignalRecord.GetRecordBuffer()        
        npBuffer = numpy.core.multiarray.int_asbuffer(ctypes.addressof(cBuffer[0].contents), 8*cBuffer[1])
        record = numpy.frombuffer(npBuffer, float)
        
        startDateTime = ole2datetime(ISignalRecord.GetStartTime()) #.strftime("%a, %d %b %Y %H:%M:%S +0000")      
        for channel in self.labels:
            nbSamples = self.channelFreqs[channel]*timeDuration       
            if channel in signalNames:   
                returnData[channel]                = RecordedChannel()                
                returnData[channel].signal         = array(record[indS:int(indS+nbSamples)])
                returnData[channel].samplingRate   = self.channelFreqs[channel]
                returnData[channel].type           = self.channelType[channel]
                returnData[channel].startTime      = startDateTime                   
                
            indS += int(nbSamples)

        return returnData 






    def getChannelTime(self, channel) :
        if not (isinstance(channel, str) or isinstance(channel, unicode)):
            print type(channel)
            raise TypeError        
        
        #if not channel in self.getChannelLabels() :
        #    raise 

        discontinuityEvent = [e for e in self.events if e.groupName == "Discontinuity" or
                                                        e.groupName == "Recording Start"]      
                                                        
        discontinuitySample = array([e.startSample for e in discontinuityEvent])                                                          
        # These sample are in base frequency. We must convert them in channel frequency.
        discontinuitySample = discontinuitySample*self.channelFreqs[channel]/self.trueBaseFreq 
                                                        
        sampleTransitions = concatenate((discontinuitySample, [self.getNbSample(channel)]))

        time = []
        samplingRate = self.getChannelFreq(channel)
        for i in range(len(sampleTransitions)-1):
            nbSamples = sampleTransitions[i+1] - sampleTransitions[i]
            time = concatenate((time, discontinuityEvent[i].startTime + arange(nbSamples)/samplingRate))       

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



    def getNbPages(self):
        return len(self.getInfoPages())



    def saveAsEDF(self, fileName=None, fileType = "EDF", isSplitted=True, 
                       annotationFileName=None, channelList=None,
                       basePath=None, annotationBasePath=None, verbose=True):
        
        """
         Convert and save the .sig/.sts files in EDF/BDF format. 
         
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
         regardless of the value of annotationFileName.
         
         fileName   : File name to use to record the converted file. If a 
                      basePath is given, the method only use the file name
                      and append it to basePath. If left to None, use the name
                      of the source file, replacing its extension.
                      
         channelList: If different than None, only the channel in that list 
                      will be saved in the converted file. 
                      
         basePath   : If different than None, record the in this path. Else,
                      use the same path as the source file or as in fileName. 
        """        

        if fileName is None:
            # Use the same name as the source name, chnaging the extension.
            fileName = ".".join(self.fileName.split(".")[:-1] + [fileType.lower()])
  
            
        if not basePath is None:
            fileName = os.path.join(basePath, os.path.basename(fileName))
        

        if isSplitted:        
            if annotationFileName is None:
                annotationFileName = fileName + "a"  
                
            if not annotationBasePath is None:
                annotationFileName = os.path.join(annotationBasePath, os.path.basename(annotationFileName))                
        

        

        
        #######################################################################
        # Internal methods
        #######################################################################        
        
        def subjectFields(field):
            return field if field else "X"
            
        # This conversion is needed because timedate.strftime("%d-%b-%Y") would
        # use local month abbreviation instead of english ones.
        def dateStr(OLEdate) :
            month = {1:"JAN", 2:"FEB", 3:"MAR", 4:"APR", 5:"MAY", 6:"JUN", 
                     7:"JUL", 8:"AUG", 9:"SEP", 10:"OCT", 11:"NOV", 12:"DEC"}
            return "%02d-%s-%04d" % (OLEdate.day, month[OLEdate.month], OLEdate.year)
        

        def prepareEventStr(timeDiff, events):
            # Annotation channel            
            timeKeepingStr = "+" + str(timeDiff) + "\x14\x14\0"       
                            
            for event in events:
                timeKeepingStr += event.toEDFStr() #edfEventEncode(event)  

            return timeKeepingStr
 
 
        #######################################################################
        # Encoding events
        #######################################################################
 
        # Encoding events such that events are writen in their respective page.
        # This has the disaventage to use a lot of space (each pages have
        # is reserved enough space to record as much annotation as the 
        # page with the most annotation) but to make complete data records.
        if verbose:
            print "Encoding events..."

        ISignalRecord = self.ISignalFile.CreateSignalRecord(RECORD_NB)   
        for page in self.getInfoPages():
              
            if self.getNbSample() - page.startSample >= RECORD_NB: 
                start = page.startSample
                ISignalRecord.SetStartSample(page.startSample)
            else:
                start = self.getNbSample() - RECORD_NB
                ISignalRecord.SetStartSample(self.getNbSample() - RECORD_NB)

            # This function call is the bottle neck of this section.
            try:             
                self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)
            except:
                print "error: ", self.getNbSample(), start
                raise

            # Annotation channel            
            page.startTime = (ole2datetime(ISignalRecord.GetStartTime()) - self.recordingStartTime).total_seconds()         

            
            # This function call is the bottle neck of this section.
            self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)

            # Annotation channel            
            page.startTime = (ole2datetime(ISignalRecord.GetStartTime()) - self.recordingStartTime).total_seconds()         

            if self.getNbSample() - page.startSample < RECORD_NB: 
               page.startTime += (page.startSample - (self.getNbSample() - RECORD_NB))/self.trueBaseFreq            




        startTime = self.getInfoPages()[1].startTime
        indMax = max(0, min(len(self.events)-1, self.events.getIndexMaxStartTime(startTime, inclusive=False)))     
        #filteredEvents = filter(lambda e: e.startTime < self.getInfoPages()[1].startTime, self.events)    
        filteredEvents = self.events[:(indMax+1)]
        self.getInfoPages()[0].eventStr = prepareEventStr(self.getInfoPages()[0].startTime, filteredEvents)
        nbEvents = len(filteredEvents)  
        #print 0, indMax, len(filteredEvents), startTime,  [e.startTime for e in filteredEvents]  
            
        for nopage in range(1, len(self.getInfoPages())-1):

            indMax = min(len(self.events)-1, self.events.getIndexMaxStartTime(self.getInfoPages()[nopage+1].startTime, inclusive=False))   
            indMin = max(0, self.events.getIndexMinStartTime(self.getInfoPages()[nopage].startTime))
            
            filteredEvents = self.events[indMin:(indMax+1)]          
            #filteredEvents = filter(lambda e: e.startTime <  self.getInfoPages()[nopage+1].startTime and 
            #                                  e.startTime >= self.getInfoPages()[nopage].startTime, self.events)   
            startTime = self.getInfoPages()[nopage].startTime
            self.getInfoPages()[nopage].eventStr = prepareEventStr(startTime, filteredEvents)                                           
            nbEvents += len(filteredEvents)  
            #print indMin, indMax, len(filteredEvents) , startTime,  [e.startTime for e in filteredEvents]  
            
        startTime = self.getInfoPages()[len(self.getInfoPages())-1].startTime
        indMin = max(0, self.events.getIndexMinStartTime(startTime))
        filteredEvents = self.events[indMin:]        
        #filteredEvents = filter(lambda e: e.startTime >= self.getInfoPages()[len(self.getInfoPages())-1].startTime, self.events)     
        self.getInfoPages()[len(self.getInfoPages())-1].eventStr = prepareEventStr(self.getInfoPages()[len(self.getInfoPages())-1].startTime, filteredEvents)        
        nbEvents += len(filteredEvents)  
        #print indMin, len(self.events)-1, len(filteredEvents), startTime,  [e.startTime for e in filteredEvents]         
        
        #print "res:", nbEvents, len(self.events)
        assert(nbEvents == len(self.events))            
            





        #######################################################################
        # Building the EDF Header object
        #######################################################################

        dataHeader = EDFHeader()       

        dataHeader.fileType = fileType
        dataHeader.contiguous = False

        id      = subjectFields(self.patientInfo["Id1"])
        gender  = "M" if self.patientInfo["gender"] == 1 else "F"
        date    = dateStr(ole2datetime(self.patientInfo["birthDate"])) if self.patientInfo["birthDate"] else "X"     
        fname   =  subjectFields(self.patientInfo["firstName"]).encode("ascii", "ignore")
        lname   = subjectFields(self.patientInfo["lastName"]).encode("ascii", "ignore")            
        dataHeader.subjectID = id + " " + gender + " " +  date + " " + fname + "_" + lname   
        
        startdate    = dateStr(self.recordingStartTime) if self.recordingStartTime else "X"                     
        dataHeader.recordingIR = "Startdate " + startdate + " " + "X" + " " + "X" + " " + "X"  

        dataHeader.startDateTime = self.recordingStartTime
            
        if channelList is None:
            acceptedChannels = self.getChannelLabels()   
        else:
            acceptedChannels = [channel for channel in self.getChannelLabels() if channel in channelList]            
            
            
        ns = len(acceptedChannels)
        dataHeader.headerNbBytes      = 8 + 80 + 80 + 8 + 8 + 8 + 44 + 8 + 8 + 4 + (ns+1)* (16 + 80+ 8 + 8 + 8 + 8 + 8 + 80 + 8 + 32) 
        dataHeader.nbRecords          = self.getNbPages()
        dataHeader.recordDuration     = self.pageDuration
        dataHeader.nbChannels         = ns +1
        
        dataHeader.channelLabels      = acceptedChannels + ["EDF Annotations"]        
        
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
        
        if fileType == "EDF":
            nbByte = 2
        elif fileType == "BDF":
            nbByte = 3

        annotationFieldLength = int(max(400, max(array([len(page.eventStr) for page in self.getInfoPages()]))/nbByte*1.2))             
        
        for i, channel in enumerate(dataHeader.channelLabels) :   
            
            if fileType == "EDF":
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
                    dataHeader.nbSamplesPerRecord[channel] = 10
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
                dataHeader.nbSamplesPerRecord[channel]  = self.pageNbSamples[channel]




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
                    print "Writing annotation file header..."
                
                annotationHeader.write(f)

                if verbose:
                    print "Writing annotation file body..."
        
                for nopage, page in enumerate(self.getInfoPages()):            
                    encodedPageStr = page.eventStr.encode("utf8")
                    f.write(encodedPageStr +  "\0"*(annotationFieldLength*nbByte-len(encodedPageStr)))         

               
       
       
       
       
       
        #######################################################################
        # Writing the data file...
        #######################################################################       
       
        # Using buffered writer
        with io.open(fileName, 'wb') as f:


            if verbose:
                print "Writing data file header..."
            dataHeader.write(f)



            if verbose:
                print "Writing date file body..."
                        
            ISignalRecord = self.ISignalFile.CreateSignalRecord(RECORD_NB) 
            for nopage, page in enumerate(self.getInfoPages()):    
                self.ISignalFile.InitSignalRecord(page.getNbSamples(), ISignalRecord) 
                
                ISignalRecord.SetStartSample(page.startSample)

                # This function call is the bottle neck of this section.
                self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)

                cBuffer = ISignalRecord.GetRecordBuffer()        
                record = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(ctypes.addressof(cBuffer[0].contents), 8*cBuffer[1]), float) 
                #record = ISignalRecord.GetRecordData()
                # record is in uV

                #print nopage+1, f.tell()
                                        
                indS = 0                
                for channel in self.labels:
                    
                    # If this page is discontinuous within the .sig file, only 
                    # the part up to the discontinuity is used: the rest of 
                    # this EDF signal is padded with zeros to make it the correct
                    # size of a complete record.                    
                    if not page.isComplete:
                        # Get the number of sample associate to the incomplete pages for that channel 

                        channelPageNbSample = int(page.getNbSamples()/self.trueBaseFreq*self.channelFreqs[channel])  
                        
                        recordedSignal = append(array(record[indS:(indS+channelPageNbSample)]), 
                                                zeros(self.pageNbSamples[channel] - channelPageNbSample))
                        indS += channelPageNbSample
                        
                    else:
                        recordedSignal = array(record[indS:(indS+self.pageNbSamples[channel])])
                        indS += self.pageNbSamples[channel]
                    

                    if channel in acceptedChannels:
                                            
                        # WRITE RECORDED SIGNAL....
                        physical_min = physicalMinMicro[channel] 
                        physical_max = physicalMaxMicro[channel] 
                        digital_min  = digitalMin[channel] 
                        digital_max  = digitalMax[channel] 
    
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
                            recordedSignal = "".join([recordedSignal[noBit] for noBit in range(len(recordedSignal)) if (noBit+1)%4])
                            f.write(recordedSignal)
                                                    
               
                
                # Annotation channel            
                if isSplitted:
                     eventStr = page.eventStr.split("\x14\x14\0")[0] + "\x14\x14\0"
                else:
                     eventStr = page.eventStr
                     
                encodedPageStr = eventStr.encode("utf8")
                f.write(encodedPageStr +  "\0"*(dataHeader.nbSamplesPerRecord["EDF Annotations"]*nbByte-len(encodedPageStr)))         

                if verbose:
                    done=float(nopage)/len(self.getInfoPages())*100.0
                    stdout.write(" Body writing percentage: %s%%      %s"%(done,"\r"))
                    stdout.flush()
        


