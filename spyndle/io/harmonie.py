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


import comtypes
import comtypes.client as cc
import os
import ctypes, numpy
from scipy.io import savemat
from scipy import array, ceil, zeros, append, arange, concatenate
from math import floor
from datetime import datetime, timedelta
import time


from EEGDatabaseReader import EEGDBReaderBase, Event, RecordedChannel, EEGPage, EEGPageInfo

OLE_TIME_ZERO = datetime(1899, 12, 30, 0, 0, 0)
def ole2datetime(oledt):
    return OLE_TIME_ZERO + timedelta(days=float(oledt))

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
    #print     seconds/24.0/60.0/60.0, seconds/24.0/60.0/60.0 + startDay
    return seconds/24.0/60.0/60.0 + startDay
    



class HarmonieEvent(Event):
    def __init__(self, no, groupe, item, recordingStartTime):
        self.no          = no
        self.groupeName  = groupe.GetName()
        self.channel     = item.GetChannel()
        
        self.name = item.GetName()  
        if self.groupeName.lower() == "stage":
            if  self.name.lower() == "stage1":   
                self.name = "Sleep stage 1"
            elif  self.name.lower() == "stage2":   
                self.name = "Sleep stage 2"
            elif  self.name.lower() == "stage3":   
                self.name = "Sleep stage 3"
            elif  self.name.lower() == "stage4":   
                self.name = "Sleep stage 4"
            elif  self.name.lower() == "rem":   
                self.name = "Sleep stage R"
            elif  self.name.lower() == "wake":   
                self.name = "Sleep stage W"
            else:   
                self.name = "Sleep stage ?"

        
        # OLE complete date-time object giving the begining time of the event.
        self.dateTime    = ole2datetime(item.GetStartTime())   
          
        # In seconds, since the begining of the recording of the EEG file.
        self.startTime   = (self.dateTime - recordingStartTime).total_seconds()  
        
        self.timeLength  = item.GetTimeLength()                
        self.startSample = item.GetStartSample()
        self.sampleLength= item.GetSampleLength()
        self.color       = groupe.GetColor() 

        self.properties  = {}
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
        #print ISignalFile.IsValid(fileName)

        self.fname = fname

        tlb_id = comtypes.GUID("{72A34744-DDD9-11D1-BB8F-00001B4E6868}")
        SignalFile_id = comtypes.GUID('{72A34754-DDD9-11D1-BB8F-00001B4E6868}')
        cc.GetModule((tlb_id, 1, 1))
        
        import comtypes.gen.SignalFileLib as SignalLib
        
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

        self.basePageNbSamples  = int(self.IRecordingMontage.GetBaseSampleFrequency()*self.pageDuration)
        self.trueBaseFreq       = self.IRecordingCalibration.GetTrueSampleFrequency()        
        self.baseFreq           = self.IRecordingMontage.GetBaseSampleFrequency()             
        self.nbChannels         = self.IRecordingMontage.GetChannelCount()
        
        self.labels             = [self.IRecordingMontage.GetChannelLabel(i)  for i in range(self.nbChannels)]    
        
        self.channelFreqs = {}
        for i in range(self.nbChannels):
            self.channelFreqs[self.labels[i]] = self.IRecordingMontage.GetChannelSampleFrequency(i)
            
        self.pageNbSamples = {}
        for channel in self.channelFreqs:
            self.pageNbSamples[channel] = int(self.channelFreqs[channel]*self.pageDuration)
            
        self.channelType = {}
        for i in range(self.nbChannels):
            self.channelType[self.labels[i]] = self.IRecordingMontage.GetChannelType(i)            
            
  
        
        self.nbSamples = self.ISignalFile.GetRecordCount(1)   
        

        ## Getting the recording start time
        ISignalRecord = self.ISignalFile.CreateSignalRecord(int(500))
        ISignalRecord.SetStartSample(int(0))
        self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)  
        
        
        self.recordingStartDateTime = ole2datetime(ISignalRecord.GetStartTime())
        
        self.startTimeInDays = ISignalRecord.GetStartTime()
        self.startDay = floor(self.startTimeInDays)
        self.recordingStartTime = dayToSecond(self.startTimeInDays, self.startDay)        
        
        
        # Store events     
        self.events = []
        for i in range(self.ISignalInfo.GetEventItemCount()):    
            IEventItem = self.ISignalInfo.GetEventItem(i)
            IEventGroup = IEventItem.GetGroup()
            self.events.append(HarmonieEvent(i, IEventGroup, IEventItem, self.recordingStartDateTime))
            
            
            
        self.events.sort(key = lambda x: x.startSample)

        #for event in self.events :
        #    print event.name, event.startSample, event.dateTime 
        
        stageEnvents = filter(lambda e: e.groupeName == "Stage", self.events)
        for event, i in zip(stageEnvents, range(len(stageEnvents))):
            event.stageEpoch = i+1
   
        self.definePages()




    """
     Deines the pages and save it in the self.page attribute. This is useful
     for managing the discontinuity. Every HarmonyPage object of self.pages
     contains it start/end sample and an indication of wheter this page ends
     with a discontinuity (i.e., it is an incomplete page).
    """
    def definePages(self):
    
        discontinuitySample = [e.startSample for e in self.events if e.groupeName == "Discontinuity" or
                                                                     e.groupeName == "Recording Start"]    
        sampleTransitions = sorted(concatenate((discontinuitySample, [self.nbSamples+1])))
                

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



    def getPages(self):
        stageEvents = filter(lambda e: e.groupeName == "Stage", self.events)    
        return map(lambda e: e.startSample, stageEvents), map(lambda e: e.startTime , stageEvents), map(lambda e: e.dateTime , stageEvents)

                
    #TODO: Optimiser...
    def getEvents(self, startTime, endTime) :   #TODO: time
        return filter(lambda e: (e.startTime >= startTime and e.startTime < endTime) or 
                         (e.startTime + e.timeLength >= startTime and e.startTime + e.timeLength < endTime) , self.events)       
        
    def getEventsBySample(self, startSample, endSample) :
        return filter(lambda e: (e.startSample >= startSample and e.startSample < endSample) or 
                         (e.startSample + e.sampleLength >= startSample and e.startSample + e.sampleLength < endSample) , self.events)       
        

    def getRecordingStartTime(self):    
        return self.recordingStartTime

        
        
    def getDuration(self):    # en secondes
        return float(self.ISignalFile.GetRecordCount(1))/self.baseFreq       

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
            nbSamples = int(self.getInfoPages(pageNo).getNbSamples()/self.baseFreq*self.channelFreqs[channel])  
            if channel in signalNames:
                # Get the number of sample associate to this channel 
                recordedSignals[channel] = array(record[indS:(indS+nbSamples)])
                
            indS += nbSamples

        # TODO: Should return a "page" object which is common to all reader.
        return EEGPage(self.channelFreqs, sigStart, recordedSignals, self.recordingStartDateTime)
     
     
     
     

    def getSamplingRate(self, channel) :
        for i in range(len(self.labels)):
            if self.labels[i] == channel:
                return self.channelFreqs[i] 
        return 0.0
     
     
     
    def getSamplingRates(self):
        return self.channelFreqs
     
    def getChannelTypes(self):
        return self.channelType


    def readChannel(self, signalName, usePickled=True):
        if usePickled:
            return self.readPickledChannel(signalName)
        else:
            raise NotImplementedError

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
                
                startDateTime = ole2datetime(ISignalRecord.GetStartTime()).strftime("%a, %d %b %Y %H:%M:%S +0000")   
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
                channelNbSample = int(nbSamples/self.baseFreq*self.channelFreqs[channel])  
                
                if channel in signalNames:
                    channelStartSample = int(startSample/self.baseFreq*self.channelFreqs[channel])   
                    # Get the number of sample associate to this channel 
                    print noPass, nbSamples, channel, channelStartSample, channelNbSample
                    returnData[channel].signal[channelStartSample:(channelStartSample+channelNbSample)] = array(record[indS:(indS+channelNbSample)])
                        
                indS += channelNbSample
                
            assert(indS ==  len(record))



        return returnData



    # La fonctionnalité ISignalRecord.SetStartTime ne semble pas pouvoir être utilisée
    # pour la lecture.
    def readWithTime(self, signalNames, startTime, timeDuration):
                
        recordNbSample = timeDuration*self.baseFreq
        ISignalRecord = self.ISignalFile.CreateSignalRecord(int(recordNbSample))
        ISignalRecord.SetStartTime(secondToDays(startTime, self.startDay))

        #print secondToDays(startTime, self.startDay), ISignalRecord.GetStartTime(), self.startTimeInDays

        return self.read2(signalNames, timeDuration, ISignalRecord)



    def read(self, signalNames, startSample, timeDuration):


        # TODO: Implémenter le code pour obtenir la valeur de sigStart
        startSample = max(startSample, 0)  
        startSample = min(startSample, self.nbSamples-timeDuration*self.baseFreq)  
        
        recordNbSample = timeDuration*self.baseFreq
        ISignalRecord = self.ISignalFile.CreateSignalRecord(int(recordNbSample))
        ISignalRecord.SetStartSample(int(startSample))
        
        return self.read2(signalNames, timeDuration, ISignalRecord)
        
        
        
    def read2(self, signalNames, timeDuration, ISignalRecord):    
        raise DeprecationWarning
        self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)
 
    
        indChannels   = [ind for lab, ind in zip(self.labels, range(len(self.labels))) if lab in signalNames]    
        channelFreqs = [sr  for lab, sr  in zip(self.labels, self.channelFreqs)      if lab in signalNames] 
        channelTypes  = [ct  for lab, ct  in zip(self.labels, self.channelType)        if lab in signalNames] 

    
        recordedSignals = [array([]) for i in range(len(indChannels))]
        indS = 0
        
       # print indChannels, channelFreqs, channelTypes, signalNames, self.labels
                
        cBuffer = ISignalRecord.GetRecordBuffer()        
        npBuffer = numpy.core.multiarray.int_asbuffer(ctypes.addressof(cBuffer[0].contents), 8*cBuffer[1])
        record = numpy.frombuffer(npBuffer, float)
        
        for i in range(self.nbChannels):
            nbSamples = self.channelFreqs[i]*timeDuration       
            if i in indChannels:   
                recordedSignals[indChannels.index(i)] = array(record[indS:int(indS+nbSamples)])
                
            indS += int(nbSamples)

        return (channelFreqs, recordedSignals, channelTypes, ole2datetime(ISignalRecord.GetStartTime()), array(self.labels)[indChannels])






    def getChannelTime(self, channel) :
        if not (isinstance(channel, str) or isinstance(channel, unicode)):
            print type(channel)
            raise TypeError        
        
        #if not channel in self.getChannelLabels() :
        #    raise 

        correctionFactor = self.trueBaseFreq/self.baseFreq 

        samplingRate = self.channelFreqs[channel]*correctionFactor
        
        discontinuityEvent = [e for e in self.events if e.groupeName == "Discontinuity" or
                                                        e.groupeName == "Recording Start"]      
                                                        
        discontinuitySample = array([e.startSample for e in discontinuityEvent])                                                          
        # These sample are in base frequency. We must convert them in channel frequency.
        discontinuitySample = discontinuitySample*self.channelFreqs[channel]/self.baseFreq 
                                                        
        sampleTransitions = concatenate((discontinuitySample, [self.getNbSample(channel)]))

        time = []
        for i in range(len(sampleTransitions)-1):
            nbSamples = sampleTransitions[i+1] - sampleTransitions[i]
            time = concatenate((time, discontinuityEvent[i].startTime + arange(nbSamples)/samplingRate))       

        return time


    # Patch parce qu'il semble y avoir un problème avec le temps associé
    # aux fuseaux de Gaétan
    def resyncEvents(self):
        ISignalRecord = self.ISignalFile.CreateSignalRecord(2)
        for event in self.events:
            ISignalRecord.SetStartSample(event.startSample)
            self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)
            event.startTime   = dayToSecond(ISignalRecord.GetStartTime(), self.startDay)  #en secondes #dayToTime(item.GetStartTime())
            event.dateTime    = ole2datetime(ISignalRecord.GetStartTime())    


    # Patch parce qu'il semble y avoir un problème avec le temps associé
    # aux fuseaux de Gaétan
    def getSampleFromTime(self, approximativeSample, searchedDatetime):
        ISignalRecord = self.ISignalFile.CreateSignalRecord(1)
        ISignalRecord.SetStartSample(approximativeSample)
        self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)
        guessedDatetime    = ole2datetime(ISignalRecord.GetStartTime())   
        delta = searchedDatetime - guessedDatetime
        deltaSec = delta.days*24.0*3600.0 + delta.seconds + delta.microseconds/1000000.0
        print guessedDatetime, searchedDatetime, delta, deltaSec, int(round(deltaSec*self.baseFreq) + approximativeSample)
        return int(round(deltaSec*self.baseFreq) + approximativeSample)
        
            
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
    def saveAsEDF(self, filename, fileType = "EDF"):
        
        def subjectFields(field):
            return field if field else "X"
            
        # This conversion is needed because timedate.strftime("%d-%b-%Y") would
        # use local month abbreviation instead of english ones.
        def dateStr(OLEdate) :
            month = {1:"JAN", 2:"FEB", 3:"MAR", 4:"APR", 5:"MAY", 6:"JUN", 
                     7:"JUL", 8:"AUG", 9:"SEP", 10:"OCT", 11:"NOV", 12:"DEC"}
            return "%02d-%s-%04d" % (OLEdate.day, month[OLEdate.month], OLEdate.year)
        
        
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
        
        

        with open(filename, 'wb') as f:

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
                f.write(" "*44) 
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
            for i in range(ns+1): f.write(" "*8)

            #print self.IRecordingCalibration.GetBaseCalibration(SIGNALFILE_FLAGS_CALIBRATEASVOLTS)   
            #print self.IRecordingCalibration.GetBaseCalibration(SIGNALFILE_FLAGS_BASEINPUTCALIB)            
            #for i in range(ns+1): 
            #     print "channel ", i, ":", self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATEASVOLTS)
            #     print "channel ", i, ":", self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATE)
            
            
            # ns * 8 ascii : ns * physical minimum (e.g. -500 or 34)
            for i in range(ns): 
                f.write(("%8.6f" %  self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATEASVOLTS)[0])[:8])
            f.write("%08d" %  -1)     


            # ns * 8 ascii : ns * physical maximum (e.g. 500 or 40)
            for i in range(ns): 
                f.write(("%8.6f" %  self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATEASVOLTS)[2])[:8])
            f.write("%08d" %  1)     


            # ns * 8 ascii : ns * digital minimum (e.g. -2048)
            for i in range(ns): 
                f.write("%08d" %  self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATE)[0])
            f.write("%08d" %  -32768)        


            # ns * 8 ascii : ns * digital maximum (e.g. 2047)
            for i in range(ns): 
                f.write("%08d" %  self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATE)[2])
            f.write("%08d" %  32767)




            # ns * 80 ascii : ns * prefiltering (e.g. HP:0.1Hz LP:75Hz)
            for i in range(ns+1): f.write(" "*80)


            # ns * 8 ascii : ns * nr of samples in each data record
            for channel in self.labels: 
                f.write("%08d" % self.pageNbSamples[channel])
            annotationFieldLength = 400      
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
            nbEvents = len(self.events)
            noEvent  = 0
            eventStr = edfEventEncode(self.events[noEvent])    
            if fileType == "EDF":
                nbByte = 2
            elif fileType == "BDF":
                nbByte = 3
                        
            for page in self.getInfoPages():
                ISignalRecord = self.ISignalFile.CreateSignalRecord(page.getNbSamples())     
                
                ISignalRecord.SetStartSample(page.startSample)
                self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)
                
                # TODO: Évaluer la possibilité d'utiliser GetRecordBuffer pour augmenter l'efficacité...
                record = ISignalRecord.GetRecordData()       

                indS = 0                
                for i, channel in enumerate(self.labels):
                    # Get the number of sample associate to every channels 
                    channelPageNbSample = int(page.getNbSamples()/self.baseFreq*self.channelFreqs[channel])                      
                    
                    # If this page is discontinuous within the .sig file, only 
                    # the part up to the discontinuity is used: the rest of 
                    # this EDF signal is padded with zeros to make it the correct
                    # size of a complete record.                    
                    if not page.isComplete:
                        recordedSignal = append(array(record[indS:(indS+channelPageNbSample)]), 
                                                zeros(self.channelFreqs[channel]*self.pageDuration - channelPageNbSample))
                    else:
                        recordedSignal = array(record[indS:(indS+channelPageNbSample)])

                    indS += channelPageNbSample
                    
                    
                    # WRITE RECORDED SIGNAL....
                    physical_min = self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATEASVOLTS)[0]
                    physical_max = self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATEASVOLTS)[2]
                    digital_min  = self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATE)[0]
                    digital_max  = self.IRecordingCalibration.GetChannelCalibration(i, SIGNALFILE_FLAGS_CALIBRATE)[2]

                    phys_range = physical_max - physical_min
                    dig_range = digital_max - digital_min
                    assert numpy.all(phys_range > 0)
                    assert numpy.all(dig_range > 0)
                    gain = phys_range / dig_range                          
                    
                    recordedSignal = (recordedSignal - physical_min)/gain + digital_min       

                    if nbByte == 2: # EDF
                        f.write(recordedSignal.astype('<h').tostring())
                    elif nbByte == 3: #BDF
                        # Writing in a string of 32-bit integers and removing every fourth byte
                        # to obtain a string of 24-bit integers
                        recordedSignal = recordedSignal.astype('<i').tostring()
                        recordedSignal = "".join([recordedSignal[noBit] for noBit in range(len(recordedSignal)) if (noBit+1)%4])
                        f.write(recordedSignal)
                        
                
                # Annotation channel            
                timeDiff = (ole2datetime(ISignalRecord.GetStartTime()) - self.recordingStartDateTime).total_seconds()  
                timeKeepingStr = "+" + str(timeDiff) + "\x14\x14\0"       
                                
                while(noEvent < nbEvents and len(timeKeepingStr) + len(eventStr) <= annotationFieldLength*nbByte):
                    timeKeepingStr += eventStr
                    noEvent += 1
                    if noEvent < nbEvents :
                        eventStr = edfEventEncode(self.events[noEvent])  
                    
                f.write(timeKeepingStr +  "\0"*(annotationFieldLength*nbByte-len(timeKeepingStr)))         


        if noEvent < nbEvents:
            raise IOError("Not enough place in EDF Annotation signal to store all events.")
        f.closed


