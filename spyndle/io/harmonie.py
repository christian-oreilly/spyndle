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
    def __init__(self, no, groupe, item, recordingStartTime):
        self.no          = no
        self.groupeName  = groupe.GetName()
        self.channel     = item.GetChannel()
        
        self.name = item.GetName()  
        

        
        # Names for sleep stages in Harmonie are not constrained such that
        # there is no standard nqmes. Tore's lab, Julie's lab and clinical
        # lab use diffrent namings. So we must check for all possibles names.
        if self.groupeName.lower() == u"stage" or self.groupeName.lower() == u"stade":
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
                
                
            self.groupeName = u"stage"
                                                         

        
        # OLE complete date-time object giving the begining time of the event.
        self.dateTime    = ole2datetime(item.GetStartTime())   
          
        # In seconds, since the begining of the recording of the EEG file.
        self.startTime   = (self.dateTime - recordingStartTime).total_seconds()  
        
        self.timeLength  = item.GetTimeLength()                
        self.startSample = item.GetStartSample()
        self.sampleLength= item.GetSampleLength()
        self.color       = groupe.GetColor() 

        self.properties  = {}
        # Keys and properties should be in unicode
        for i in range(groupe.GetItemPropertyCount()):
            key = groupe.GetItemPropertyKey(i)
            self.properties[key] = item.GetItemPropertyValue(key)
            
    

            
    def __str__(self):
        STR =  str(self.no) + " " + str(self.groupeName) + " " + str(self.channel) \
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
        self.recordingStartDateTime = ole2datetime(ISignalRecord.GetStartTime())        
        self.startTimeInDays = ISignalRecord.GetStartTime()
        self.startDay = floor(self.startTimeInDays)
        self.recordingStartTime = dayToSecond(self.startTimeInDays, self.startDay)        
        ###########################################     

        
        ###########################################
        # Store events     
        for i in range(self.ISignalInfo.GetEventItemCount()):    
            IEventItem = self.ISignalInfo.GetEventItem(i)
            IEventGroup = IEventItem.GetGroup()
            self.events.add(HarmonieEvent(i, IEventGroup, IEventItem, self.recordingStartDateTime))
            
        # Patch used because data recorder with previous versions of Harmonie
        # inconsistent time stamp due to discontinuous recording. These lines
        # should not be necessary with another reader.
        print "Resynchronization of the events..."
        self.resyncEvents()                
        ###########################################
            
            
        self.pageDuration = np.median([e.timeLength for e in self.events if e.groupeName.lower() == "stage"])

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
    
        discontinuitySample = array([e.startSample for e in self.events if e.groupeName == "Discontinuity" or
                                                                     e.groupeName == "Recording Start"], dtype=int)    
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
        return EEGPage(self.channelFreqs, sigStart, recordedSignals, self.recordingStartDateTime)
     
     
     
    # Should return a DateTime object indicating the time when the recording
    # has begun. 

    def getChannelFreq(self, channel): 
        if channel in self.labels:
            return self.channelFreqs[channel] 
        return 0.0
    

    def getSamplingRate(self, channel) :
        warnings.warn("This function is depreacted.", DeprecationWarning, stacklevel=2)
        self.getChannelFreq(channel)      
        
     
     
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
        while i < 30 and abs(time - approxTime) > 0.5/self.trueBaseFreq:
            approxSample += int(round((time - approxTime)*self.trueBaseFreq))     
            if approxSample > self.nbSamples:
                approxSample = self.nbSamples
                
            ISignalRecord.SetStartSample(approxSample)        
            self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)
            startDateTime = ole2datetime(ISignalRecord.GetStartTime())
            
            approxTime = (startDateTime - self.recordingStartDateTime).total_seconds()  
        
        if i == 30:
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

        discontinuityEvent = [e for e in self.events if e.groupeName == "Discontinuity" or
                                                        e.groupeName == "Recording Start"]      
                                                        
        discontinuitySample = array([e.startSample for e in discontinuityEvent])                                                          
        # These sample are in base frequency. We must convert them in channel frequency.
        discontinuitySample = discontinuitySample*self.channelFreqs[channel]/self.trueBaseFreq 
                                                        
        sampleTransitions = concatenate((discontinuitySample, [self.getNbSample(channel)]))

        time = []
        samplingRate = self.getSamplingRate(channel)
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
            self.startTime   = (event.dateTime - self.recordingStartDateTime).total_seconds()     


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


    """
     Convert and save the .sig/.sts files in .edf/.ann files with the .ann file
     being a custom format file used to save complex annotations.
    """
    def saveAsEDF(self, filename, fileType = "EDF", verbose=True):
        
        def subjectFields(field):
            return field if field else "X"
            
        # This conversion is needed because timedate.strftime("%d-%b-%Y") would
        # use local month abbreviation instead of english ones.
        def dateStr(OLEdate) :
            month = {1:"JAN", 2:"FEB", 3:"MAR", 4:"APR", 5:"MAY", 6:"JUN", 
                     7:"JUL", 8:"AUG", 9:"SEP", 10:"OCT", 11:"NOV", 12:"DEC"}
            return "%02d-%s-%04d" % (OLEdate.day, month[OLEdate.month], OLEdate.year)
        
        """
        # Function used to encode Harmony events in a format compatible with EDF.
        def edfEventEncode(event):
            
            # Using standard EDF staging events 
            #(see http://www.edfplus.info/specs/edftexts.html#annotation)

            if event.groupeName.lower() == "stage":
                eventStr = event.name
            else:
                eventStr = event.getXml()
            
            timeDiff = (event.dateTime - self.recordingStartDateTime).total_seconds()  
            
            return "+" + str(timeDiff) + "\x15" + str(event.timeLength) + "\x14" + eventStr + "\x14\0"    
        """

        def prepareEventStr(timeDiff, events):
            # Annotation channel            
            timeKeepingStr = "+" + str(timeDiff) + "\x14\x14\0"       
                            
            for event in events:
                timeKeepingStr += event.toEDFStr() #edfEventEncode(event)  

            return timeKeepingStr
 
 
 
 
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
            page.startTime = (ole2datetime(ISignalRecord.GetStartTime()) - self.recordingStartDateTime).total_seconds()         

            
            # This function call is the bottle neck of this section.
            self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)

            # Annotation channel            
            page.startTime = (ole2datetime(ISignalRecord.GetStartTime()) - self.recordingStartDateTime).total_seconds()         

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
            
       
        # Using buffered writer
        with io.open(filename, 'wb') as f:

            if verbose:
                print "Writing header..."

            
            #8 ascii : version of this data format (0) 
            if fileType == "EDF":
                f.write("0       ")
            elif fileType == "BDF":
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
            id      = subjectFields(self.patientInfo["Id1"])
            gender  = "M" if self.patientInfo["gender"] == 1 else "F"
            date    = dateStr(ole2datetime(self.patientInfo["birthDate"])) if self.patientInfo["birthDate"] else "X"     
            fname   =  subjectFields(self.patientInfo["firstName"]).encode("ascii", "ignore")
            lname   = subjectFields(self.patientInfo["lastName"]).encode("ascii", "ignore")     
            writeStr = id + " " + gender + " " +  date + " " + fname + "_" + lname     
            f.write(writeStr + (80-len(writeStr))*" ")
              
            # TODO: EDF need 7-bit ASCII character which cannot for example accept
            # Benoît as a valid string. We use encode("ascii", "ignore") which only
            # drops the accentuated characters but it would be better to have some
            # translation that convert Benoît to Benoit.
              
              
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
            startdate    = dateStr(self.recordingStartDateTime) if self.recordingStartDateTime else "X"             
            writeStr = "Startdate " + startdate + " " + "X" + " " + "X" + " " + "X"  
            f.write(writeStr + (80-len(writeStr))*" ")              
              
            #8 ascii : startdate of recording (dd.mm.yy)
            f.write(self.recordingStartDateTime.strftime("%d.%m.%y"))
            
            
            #8 ascii : starttime of recording (hh.mm.ss) 
            f.write(self.recordingStartDateTime.strftime("%H.%M.%S"))

            # 8 ascii : number of bytes in header record
            ns = len(self.getChannelLabels())
            headerSize = 8 + 80 + 80 + 8 + 8 + 8 + 44 + 8 + 8 + 4 + (ns+1)* (16 + 80+ 8 + 8 + 8 + 8 + 8 + 80 + 8 + 32)
            f.write("%08d" % headerSize)   
            
            # 44 ascii : reserved
            if fileType == "EDF":
                f.write("EDF+C" + " "*39)
            elif fileType == "BDF":
                f.write("24BIT" + " "*39)            
            
            
            # 8 ascii : number of data records (-1 if unknown)
            #  The 'number of data records' can only be -1 during recording. 
            # As soon as the file is closed, the correct number is known and must be entered. 
            nbPages = self.getNbPages()
            f.write("%08d" % nbPages)  
            
            # 8 ascii : duration of a data record, in seconds
            f.write(("%8.6f" % self.pageDuration)[:8] )  
            
            # 4 ascii : number of signals (ns) in data record
            f.write("%04d" % (ns +1))
            
                
            # ns * 16 ascii : ns * label (e.g. EEG Fpz-Cz or Body temp)
            for i in range(ns): 
                f.write("%16s" % self.getChannelLabels()[i])
            f.write("%16s" % "EDF Annotations")

  
            # ns * 80 ascii : ns * transducer type (e.g. AgAgCl electrode)
            for i in range(ns+1): f.write(" "*80)

            # ns * 8 ascii : ns * physical dimension (e.g. uV or degreeC)
            for i in range(ns+1): f.write("      uV")

            #print self.IRecordingCalibration.GetBaseCalibration(SIGNALFILE_FLAGS_CALIBRATEASVOLTS)   
            #print self.IRecordingCalibration.GetBaseCalibration(SIGNALFILE_FLAGS_BASEINPUTCALIB)            
            #for i in range(ns+1): 
            #     print "channel ", i, ":", self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATEASVOLTS)
            #     print "channel ", i, ":", self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATE)
            
            
            # ns * 8 ascii : ns * physical minimum (e.g. -500 or 34)
            physicalMinMicro = {}
            physicalMaxMicro = {}
            for i in range(ns): 
                inputMin, outputMin, inputMax, outputMax = self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATEASVOLTS)
                physicalMaxMicro[self.getChannelLabels()[i]]  = inputMax*(outputMax-outputMin)/(inputMax - inputMin)  
                physicalMinMicro[self.getChannelLabels()[i]]  = inputMin*(outputMax-outputMin)/(inputMax - inputMin)
                f.write(("%8.6f" %  physicalMinMicro[self.getChannelLabels()[i]] )[:8])
            f.write(("%8.6f" %  -1)[:8]) 
  
            # ns * 8 ascii : ns * physical maximum (e.g. 500 or 40)
            for i in range(ns): 
                f.write(("%8.6f" %  physicalMaxMicro[self.getChannelLabels()[i]])[:8])
            f.write(("%8.6f" %  1)[:8])  


            # ns * 8 ascii : ns * digital minimum (e.g. -2048)
            digitalMin = {}
            digitalMax = {}
            for i in range(ns): 
                if fileType == "EDF":
                    digitalMin[self.getChannelLabels()[i]] = -32768
                    digitalMax[self.getChannelLabels()[i]] =  32767
                elif fileType == "BDF":
                    digitalMin[self.getChannelLabels()[i]] = -8388608 
                    digitalMax[self.getChannelLabels()[i]] =  8388607    
                else:
                    raise ValueError
                f.write("%08d" %  digitalMin[self.getChannelLabels()[i]])
            f.write("%08d" %  -32768)        


            # ns * 8 ascii : ns * digital maximum (e.g. 2047)
            for i in range(ns): 
                f.write("%08d" %  digitalMax[self.getChannelLabels()[i]])
            f.write("%08d" %  32767)


            # ns * 80 ascii : ns * prefiltering (e.g. HP:0.1Hz LP:75Hz)
            for i in range(ns+1): f.write(" "*80)


            if fileType == "EDF":
                nbByte = 2
            elif fileType == "BDF":
                nbByte = 3


            # ns * 8 ascii : ns * nr of samples in each data record
            for channel in self.labels: 
                f.write("%08d" % self.pageNbSamples[channel])
            annotationFieldLength = int(max(400, max(array([len(page.eventStr) for page in self.getInfoPages()]))/nbByte*1.2))            
            f.write("%08d" % annotationFieldLength)
            
            
            # ns * 32 ascii : ns * reserved
            for i in range(ns+1): f.write(" "*32)

            """
             DATA RECORD
             nr of samples[1] * integer : first signal in the data record
             nr of samples[2] * integer : second signal
             ..
             ..
             nr of samples[ns] * integer : last signal
            
             N.B. Harmony files can be discontinuous at any arbitrary time
             as opposted to EDF files wich can only be discontinuous at the end
             of a record. We therfore have to put incomplete records at every 
             discontinuity to complete the EDF records. The incomplete records
             are filled with "\0" values.
            """
            
            
            #nbEvents = len(self.events)
            #noEvent  = 0
            #eventStr = edfEventEncode(self.events[noEvent])    


            if verbose:
                print "Writing body..."
                        
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

                    #print physical_min, physical_max, digital_min, digital_max
                    #inputMin, outputMin, inputMax, outputMax = self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATEASVOLTS)
                    #chanGain = (inputMax - inputMin)/(outputMax - outputMin)                         
                    #recordedSignal = (recordedSignal - outputMin)*chanGain + inputMin                      

                    #print nopage, channel, physical_min, physical_max, digital_min, digital_max, min(recordedSignal), max(recordedSignal) 
                    #print self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATEASVOLTS), self.IRecordingCalibration.GetBaseCalibration(SIGNALFILE_FLAGS_CALIBRATEASVOLTS)
                    #print self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATE), self.IRecordingCalibration.GetBaseCalibration(SIGNALFILE_FLAGS_CALIBRATE)

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
                encodedPageStr = page.eventStr.encode("utf8")
                f.write(encodedPageStr +  "\0"*(annotationFieldLength*nbByte-len(encodedPageStr)))         

                if verbose:
                    done=float(nopage)/len(self.getInfoPages())*100.0
                    stdout.write(" Body writing percentage: %s%%      %s"%(done,"\r"))
                    stdout.flush()
        
        f.closed


