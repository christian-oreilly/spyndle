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


from scipy import array, ceil, concatenate, zeros



import comtypes
import comtypes.client as cc
import os
from scipy.io import savemat, loadmat


from math import floor
from datetime import time, datetime, timedelta

import time
import ctypes, numpy

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
    
    

class EEGDBReaderBase :
    def __init__(self):
        self.pageDuration = 30 # en secondes

        
    def setPageDuration(self, duration):
        self.pageDuration = duration
    
  
class RecordedChannel:
    def __init__(self): #, ISignalFile):
        self.signal         = array([])
        self.samplingRate   = 0.0
        self.type           = ""
        self.startTime      = ""
        
    # This mimic the behavior of Python's dictionary (needed to use the 
    #savemat function).
    def items(self):
        return [("signal", self.signal), ("samplingRate", self.samplingRate), 
                ("type", self.type), ("startTime", self.startTime)]
        
        
    def readPickledChannel(self, fileName):
        data = loadmat(fileName)        
        self.signal         = data["signal"].reshape(data["signal"].size)   
        self.samplingRate   = data["samplingRate"][0][0]    
        self.type           = data["type"][0][0]  
        self.startTime      = time.strptime(data["startTime"][0], "%a, %d %b %Y %H:%M:%S +0000")

        
        
        
class Event:
    def __init__(self): #, ISignalFile):
        self.no          = ""
        self.groupeName  = ""
        self.channel     = ""
        self.name        = ""
        self.startTime   = ""
        self.dateTime    = ""      
        self.timeLength  = ""              
        self.startSample = ""
        self.sampleLength= ""
        self.color       = ""

        self.properties  = {}

    
    def sampleEnd(self):
        return self.startSample + self.sampleLength        
    def sampleStart(self):
        return self.startSample

            
    def __str__(self):
        return( str(self.no) + " " + str(self.groupeName) + " " + str(self.channel)
                + " " + str(self.name) + " " + str(self.startTime) + " " + str(self.timeLength))








class HarmonieEvent(Event):
    def __init__(self, no, groupe, item, startDay): #, ISignalFile):
        self.no          = no
        self.groupeName  = groupe.GetName()
        self.channel     = item.GetChannel()
        self.name        = item.GetName()
        self.startTime   = dayToSecond(item.GetStartTime(), startDay)  #en secondes #dayToTime(item.GetStartTime())
        self.dateTime    = ole2datetime(item.GetStartTime())        
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
        EEGDBReaderBase.__init__(self)
 
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
            # TODO: Penser à un systeme de gestion des erreurs...
            print "Erreur ouverture"
            self.ISignalFile    = None
            self.labels         = None
            self.samplingRates  = None
            self.nbChannels     = None
            self.pageNbSamples  = None
            self.basePageNbSamples = None
            self.channelType        = None
            raise
            
            
        #	Lecture de l'objet 'ISignalInfo'
        self.ISignalInfo = self.ISignalFile.GetSignalInfo()
        
        #	Lecture de l'objet 'IFileInfo'
        self.IFileInfo = self.ISignalInfo.GetFileInfo()
        
        #	Lecture de l'objet 'IPatientInfo'
        self.IPatientInfo =  self.ISignalInfo.GetPatientInfo()
        
        self.patientInfo = {"firstname":self.IPatientInfo.GetFirstName(),
                            "lastname" :self.IPatientInfo.GetLastName(),
                            "Id1"      : self.IPatientInfo.GetId1(),
                            "Id2"      : self.IPatientInfo.GetId2()}           
        
        #	Lecture de l'objet 'IRecordingCalibration'
        self.IRecordingCalibration = self.ISignalInfo.GetRecordingCalibration()
        
        #	Lecture de l'objet 'IPhysicalMontage'
        self.IPhysicalMontage = self.ISignalInfo.GetPhysicalMontage()
        
        #	Lecture de l'objet 'IRecordingMontage'
        #	ATTENTION : Il n'y a toujours qu'un seul montage d'enregistrement, à l'index '0'
        self.IRecordingMontage = self.IPhysicalMontage.GetRecordingMontage(0)

        self.nbChannels         = self.IRecordingMontage.GetChannelCount()
        self.labels             = [self.IRecordingMontage.GetChannelLabel(i)  for i in range(self.nbChannels)]        
        self.samplingRates      = [self.IRecordingMontage.GetChannelSampleFrequency(i)  for i in range(self.nbChannels)]        
        self.pageNbSamples      = array(self.samplingRates)*self.pageDuration
        self.basePageNbSamples  = self.IRecordingMontage.GetBaseSampleFrequency()*self.pageDuration
        self.channelType        = [self.IRecordingMontage.GetChannelType(i)  for i in range(self.nbChannels)] 

        #trueBaseFreq = self.IRecordingCalibration.GetTrueSampleFrequency()        
        self.baseFreq = self.IRecordingMontage.GetBaseSampleFrequency()          
        
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
            self.events.append(HarmonieEvent(i, IEventGroup, IEventItem, self.startDay))
            
        self.events.sort(key = lambda x: x.startSample)

        #for event in self.events :
        #    print event.name, event.startSample, event.dateTime 
        
        stageEnvents = filter(lambda e: e.groupeName == "Stage", self.events)
        for event, i in zip(stageEnvents, range(len(stageEnvents))):
            event.stageEpoch = i+1
   


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
        return self.ISignalFile.GetRecordCount(int(self.baseFreq))       

    def getNbSample(self):    # en secondes
        return self.nbSamples
        
    def getElectrodesLabels(self):
        return self.labels
        
    def readPage(self, signalNames, pageNo):

        ISignalRecord = self.ISignalFile.CreateSignalRecord(self.basePageNbSamples)
        ISignalRecord.SetStartSample(self.basePageNbSamples*(pageNo-1)) # On index les page de 1 à NbPages
        self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)


        # TODO: Implémenter le code pour obtenir la valeur de sigStart
        sigStart = (pageNo-1)*self.pageDuration     
    
        recordedSignals = [array([]) for i in range(self.nbChannels)]
        indS = 0
        # TODO: Évaluer la possibilité d'utiliser GetRecordBuffer pour augmenter l'efficacité...
        record = ISignalRecord.GetRecordData()
        for i in range(self.nbChannels):
            recordedSignals[i] = array(record[indS:int(indS+self.pageNbSamples[i])])
            indS += int(self.pageNbSamples[i])

        return (self.samplingRates, sigStart, recordedSignals)
     
     

    def getSamplingRate(self, channel) :
        for i in range(len(self.labels)):
            if self.labels[i] == channel:
                return self.samplingRates[i] 
        return 0.0
     
     
     
    def getSamplingRates(self):
        return self.samplingRates
     
    def getChannelTypes(self):
        return self.channelType


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






    """
        fileName = self.fname + "_readComplete.mat"

        if os.path.exists(fileName):
            print "Using pickled data..."   
            data = loadmat(fileName)
        else:
            print "Reading form .sig file..."
            data = self.readComplete(signalNames)   
            print "Pickling for next time..."
            savemat(fileName, data)
        
        return data             
    """ 

    """
    def readCompletePickle(self, signalNames):

        fileName = self.fname + "_readComplete.mat"

        if os.path.exists(fileName):
            print "Using pickled data..."   
            DATA = loadmat(fileName)
            result = self.read(signalNames, 0, 1.0)
            return (result[0], DATA["DATA"], result[2], result[3], result[4])   
        else:
            print "Reading form .sig file..."
            data = self.readComplete(signalNames)   
            print "Pickling for next time..."
            savemat(fileName, {"DATA":data[1]})
            return data             
   


    def pickleCompleteSignals(self, signalNames):

        fileName = self.fname + "_readComplete_"

        data = {}
        for signalName in signalNames:
            fileChannel = fileName + signalName + ".mat"
            if os.path.exists(fileChannel):
                print "Using pickled data for " + signalName + "..."   
                data[signalName] = loadmat(fileChannel)
            else:
                print "Reading form .sig file for " + signalName + "..."
                data[signalName] = self.readComplete(signalName)   
                print "Pickling for next time..."
                savemat(fileChannel, data[signalName])
    
        return data       

    """       
    
    def getAvailableChannels(self):
        return self.labels


    # As read, but read the complete duration of the signals.
    def readComplete(self, signalNames):

        #indChannels   = []  
        #samplingRates = []
        #channelTypes  = []
        #recordedSignals = []
        #dateTimeRecord  = None

        # For some reasons, the API refuse to create records with a too
        # large number of samples. We therefore have to split out entire
        # nights in sequential readings. 
        
        returnData = {}
        maxRecord = 500000
        nbPasses = int(ceil(self.nbSamples/float(maxRecord)))
        for noPass in range(nbPasses):
            startSample =  int(maxRecord*noPass)             
            
            if noPass == nbPasses-1:
                nbSamples = self.nbSamples - startSample
            else:
                nbSamples = maxRecord

            #print "noPass ", noPass, "/", nbPasses-1, startSample, nbSamples, self.nbSamples                
 
            if noPass == 0: 
                ISignalRecord = self.ISignalFile.CreateSignalRecord(nbSamples)                
                ISignalRecord.SetStartSample(startSample)     
                self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)   
                
                indChannels    = [ind for lab, ind in zip(self.labels, range(len(self.labels))) if lab in signalNames]    

                #returnData["samplingRates"]  = [sr  for lab, sr  in zip(self.labels, self.samplingRates)      if lab in signalNames] 
                #returnData["channelTypes"]   = [ct  for lab, ct  in zip(self.labels, self.channelType)        if lab in signalNames] 
                #returnData["dateTimeRecord"] = ole2datetime(ISignalRecord.GetStartTime())
                #returnData["signals"]        = [zeros(self.nbSamples) for i in range(len(indChannels))]
                #returnData["channels"]       = array(self.labels)[indChannels]
                
                startDateTime = ole2datetime(ISignalRecord.GetStartTime()).strftime("%a, %d %b %Y %H:%M:%S +0000")   
                for ind in indChannels:
                    channel                            = array(self.labels)[ind]
                    returnData[channel]                = RecordedChannel()                
                    returnData[channel].signal         = zeros(self.nbSamples)
                    returnData[channel].samplingRate   = self.samplingRates[ind]
                    returnData[channel].type           = self.channelType[ind]
                    returnData[channel].startTime      = startDateTime   
                    #print channel, returnData[channel].samplingRate, returnData[channel].type 
            else:
                self.ISignalFile.InitSignalRecord(nbSamples, ISignalRecord)                
                ISignalRecord.SetStartSample(startSample)  
                self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)    
                              

            cBuffer = ISignalRecord.GetRecordBuffer()        
            record = numpy.frombuffer(numpy.core.multiarray.int_asbuffer(ctypes.addressof(cBuffer[0].contents), 8*cBuffer[1]), float) 
            
            
            indS = 0
            for i in indChannels:  
                channel  = array(self.labels)[i]
                returnData[channel].signal[startSample:(startSample+nbSamples)] = array(record[indS:int(indS+nbSamples)])
                #print channel, returnData[channel].signal
                indS += int(nbSamples)

            #print (samplingRates, recordedSignals, channelTypes, dateTimeRecord)

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
        
        self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)
 
    
        indChannels   = [ind for lab, ind in zip(self.labels, range(len(self.labels))) if lab in signalNames]    
        samplingRates = [sr  for lab, sr  in zip(self.labels, self.samplingRates)      if lab in signalNames] 
        channelTypes  = [ct  for lab, ct  in zip(self.labels, self.channelType)        if lab in signalNames] 

    
        recordedSignals = [array([]) for i in range(len(indChannels))]
        indS = 0
        
       # print indChannels, samplingRates, channelTypes, signalNames, self.labels
                
        cBuffer = ISignalRecord.GetRecordBuffer()        
        npBuffer = numpy.core.multiarray.int_asbuffer(ctypes.addressof(cBuffer[0].contents), 8*cBuffer[1])
        record = numpy.frombuffer(npBuffer, float)
        
        for i in range(self.nbChannels):
            nbSamples = self.samplingRates[i]*timeDuration       
            if i in indChannels:   
                recordedSignals[indChannels.index(i)] = array(record[indS:int(indS+nbSamples)])
                
            indS += int(nbSamples)

        return (samplingRates, recordedSignals, channelTypes, ole2datetime(ISignalRecord.GetStartTime()), array(self.labels)[indChannels])




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













            
            
            
            
            
from filters import channelType



class DevuystEvent(Event):
    def __init__(self, groupeName, name, channel, startTime, timeDuration, samplingRate): #, ISignalFile):
        self.no          = None
        self.groupeName  = groupeName
        self.channel     = channel
        self.name        = name
        self.startTime   = startTime #en secondes #dayToTime(item.GetStartTime())
        self.dateTime    = None        
        self.timeLength  = timeDuration
        self.startSample = int(startTime*samplingRate)
        self.sampleLength= int(timeDuration*samplingRate)
        self.color       = None
        self.properties  = {}



class DevuystReader(EEGDBReaderBase):
    def __init__(self, fname, samplingRate):
        EEGDBReaderBase.__init__(self)
 
        #c:\tests>python makepy.py -i StlSignalFile.tlb
        #StlSignalFileLib
        # {72A34744-DDD9-11D1-BB8F-00001B4E6868}, lcid=0, major=1, minor=1
        #self.ISignalFile = Dispatch('{72A34754-DDD9-11D1-BB8F-00001B4E6868}')

        # TODO : ajouter cette vérification
        #print ISignalFile.IsValid(fileName)

        self.fname = fname

        try:
            
            with open(fname) as f:
                content = f.readlines()            
                
            content          = array(content)
            self.signal      = [content[1:].astype(float)]
            self.labels      = [content[0].split("[")[1].split("]")[0]]     
            self.nbChannels  = 1
            self.channelType = [channelType["EEG"]]  
        except:
            # TODO: Penser à un systeme de gestion des erreurs...
            print "Erreur ouverture"
            raise
            
            
        self.samplingRates  = [samplingRate]
        self.baseFreq       = self.samplingRates[0]  
        self.nbSamples      = len(self.signal[0])   
        self.events         = []
 
 

    def getAvailableChannels(self):
        return self.labels
 


    def importEvents(self, fname, eventName, groupeName="Fuseau"):

        try:
            
            with open(fname) as f:
                content = f.readlines()            

        except:
            # TODO: Penser à un systeme de gestion des erreurs...
            print "Erreur ouverture du fichier ->" + fname + "<-"
            return
                
        content = array(content[1:])
        
        for line in content :
            startTime, timeDuration = array(line.split(), float)
            self.events.append(DevuystEvent(groupeName, eventName, self.labels[0], 
                                            startTime, timeDuration, self.samplingRates[0]))



        
        
    def importHypnogram(self, fname):
        try:
            
            with open(fname) as f:
                stages = f.readlines()            
        except:
            # TODO: Penser à un systeme de gestion des erreurs...
            print "Erreur ouverture du fichier " + fname      
            return

              
        stages = array(stages[1:])

        stageLabels = ["Stage4", "Stage3", "Stage2", "Stage1", "REM", "WAKE", "Unknown"]        
        
        startTime = 0.0
        for stage in stages :
            stage = int(stage)
            if stage < 0 : stage = 6
            self.events.append(DevuystEvent("Stage", stageLabels[stage], self.labels[0], 
                                            startTime, 5.0, self.samplingRates[0]))
            startTime += 5.0
    

                
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
        return self.ISignalFile.GetRecordCount(int(self.baseFreq))       

    def getNbSample(self):    # en secondes
        return self.nbSamples
        
    def getElectrodesLabels(self):
        return self.labels
        

    def getSamplingRate(self, channel) :
        for i in range(len(self.labels)):
            if self.labels[i] == channel:
                return self.samplingRates[i] 
        return 0.0
     
     
     
    def getSamplingRates(self):
        return self.samplingRates
     
    def getChannelTypes(self):
        return self.channelType









    # La fonctionnalité ISignalRecord.SetStartTime ne semble pas pouvoir être utilisée
    # pour la lecture.
    def readWithTime(self, signalNames, startTime, timeDuration):
                
        recordNbSample = timeDuration*self.baseFreq
        ISignalRecord = self.ISignalFile.CreateSignalRecord(int(recordNbSample))
        ISignalRecord.SetStartTime(secondToDays(startTime, self.startDay))

        #print secondToDays(startTime, self.startDay), ISignalRecord.GetStartTime(), self.startTimeInDays

        return self.read2(signalNames, timeDuration, ISignalRecord)




    def readCompletePickle(self, signalNames):

        # No pickling necessary for these data...
        return self.readComplete(signalNames)
               



    # As read, but read the complete duration of the signals.
    def readComplete(self, signalNames):
        if self.labels[0] in signalNames:
            return (self.samplingRates, self.signal, self.channelType, None, self.labels)
        else:
            return ([], [], [], None, [])


    def read(self, signalNames, startSample, timeDuration=None, sampleDuration=None):
        if self.labels[0] in signalNames:

            if sampleDuration is None:
                sampleDuration = timeDuration*self.baseFreq
            elif timeDuration is None:
                timeDuration = sampleDuration/self.baseFreq 
            
            # TODO: Implémenter le code pour obtenir la valeur de sigStart
            startSample = min(startSample, self.nbSamples-sampleDuration) 
            startSample = max(startSample, 0)   
  
            return (self.samplingRates, [self.signal[0][int(startSample):int(startSample+sampleDuration)]], self.channelType, None, self.labels)
        else:
            return ([], [], [], None, [])         



    # Patch parce qu'il semble y avoir un problème avec le temps associé
    # aux fuseaux de Gaétan
    def getSampleFromTime(self, approximativeSample, searchedDatetime):
        ISignalRecord = self.ISignalFile.CreateSignalRecord(1)
        ISignalRecord.SetStartSample(approximativeSample)
        self.ISignalFile.Read(ISignalRecord, HarmonieConst.SIGNALFILE_FLAGS_CALIBRATE)
        guessedDatetime    = ole2datetime(ISignalRecord.GetStartTime())   
        delta = searchedDatetime - guessedDatetime
        deltaSec = delta.days*24.0*3600.0 + delta.seconds + delta.microseconds/1000000.0
        print guessedDatetime, searchedDatetime, delta, deltaSec, int(round(deltaSec*self.baseFreq) + approximativeSample)
        return int(round(deltaSec*self.baseFreq) + approximativeSample)
        
                        
            
            
        