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


from scipy import array, ceil, zeros, append

import comtypes
import comtypes.client as cc
import os
from scipy.io import savemat


from math import floor
from datetime import time, datetime, timedelta

import ctypes, numpy

from EEGDatabaseReader import EEGDBReaderBase, Event, RecordedChannel

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
            print "An error has occured when trying to open " + str(fname)
            self.ISignalFile    = None
            self.labels         = None
            self.samplingRates  = None
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

        self.nbChannels         = self.IRecordingMontage.GetChannelCount()
        self.labels             = [self.IRecordingMontage.GetChannelLabel(i)  for i in range(self.nbChannels)]        
        self.samplingRates      = [self.IRecordingMontage.GetChannelSampleFrequency(i)  for i in range(self.nbChannels)]        
        self.pageNbSamples      = array(self.samplingRates, dtype=numpy.int32)*self.pageDuration
        self.basePageNbSamples  = int(self.IRecordingMontage.GetBaseSampleFrequency()*self.pageDuration)
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
        return self.ISignalFile.GetRecordCount(int(self.baseFreq))       

    def getNbSample(self, channel=None):
        if channel is None:
            return self.nbSamples
        else:
            if not isinstance(channel, str):
                raise TypeError 
            return self.getDuration()*self.samplingRates[self.labels.index(channel)] 
        
    def getElectrodesLabels(self):
        return self.labels
        
    def readPage(self, signalNames, pageNo):

        ISignalRecord = self.ISignalFile.CreateSignalRecord(self.basePageNbSamples)
        ISignalRecord.SetStartSample(self.basePageNbSamples*(pageNo-1)) # On index les page de 1 à NbPages
        self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)


        # TODO: Implémenter le code pour obtenir la valeur de sigStart
        sigStart = (pageNo-1)*self.pageDuration     
    
        recordedSignals = [array([]) for i in range(len(signalNames))]
        indS = 0
        # TODO: Évaluer la possibilité d'utiliser GetRecordBuffer pour augmenter l'efficacité...
        record = ISignalRecord.GetRecordData()
        i2 = 0
        for i in range(self.nbChannels):
            if self.getChannelLabels()[i] in signalNames:
                recordedSignals[i2] = array(record[indS:int(indS+self.pageNbSamples[i])])
                i2 += 1
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



    def getNbPages(self):
        return int(ceil(self.nbSamples/float(self.basePageNbSamples)))


    """
     Convert and save the .sig/.sts files in .edf/.ann files with the .ann file
     being a custom format file used to save complex annotations.
    """
    def saveAsEDFA(self, filename, fileType = "EDF"):
        
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
                if  event.name.lower() == "stage1":   
                    eventStr = "Sleep stage 1"
                elif  event.name.lower() == "stage2":   
                    eventStr = "Sleep stage 2"
                elif  event.name.lower() == "stage3":   
                    eventStr = "Sleep stage 3"
                elif  event.name.lower() == "stage4":   
                    eventStr = "Sleep stage 4"
                elif  event.name.lower() == "REM":   
                    eventStr = "Sleep stage R"
                elif  event.name.lower() == "WAKE":   
                    eventStr = "Sleep stage W"
                else:   
                    eventStr = "Sleep stage ?"
            else:
                eventStr = event.getXml()
            
            
            
            return "+" + str(event.startTime) + "\x15" + str(event.timeLength) + "\x14" + eventStr + "\x14\0"    
        
        
        
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
            for i in range(ns): 
                f.write("%08d" % self.pageNbSamples[i])
            annotationFieldLength = 400      
            f.write("%08d" % annotationFieldLength)
            
            
            # ns * 32 ascii : ns * reserved
            for i in range(ns+1): f.write(" "*32)

            # DATA RECORD
            # nr of samples[1] * integer : first signal in the data record
            # nr of samples[2] * integer : second signal
            # ..
            # ..
            # nr of samples[ns] * integer : last signal
            
            nbEvents = len(self.events)
            noEvent  = 0
            eventStr = edfEventEncode(self.events[noEvent])    
            if fileType == "EDF":
                nbByte = 2
            elif fileType == "BDF":
                nbByte = 3
                        
            
            ISignalRecord = self.ISignalFile.CreateSignalRecord(self.basePageNbSamples)
            for noPage in range(nbPages):
                
                # The last page is not necessary complete in Harmony so we need to 
                # define a record of a correct size to read the last samples without
                # overflowing.
                if noPage  == nbPages -1:
                    ISignalRecord = self.ISignalFile.CreateSignalRecord(self.nbSamples - self.basePageNbSamples*(nbPages-1))
                    
                ISignalRecord.SetStartSample(self.basePageNbSamples*noPage)
                self.ISignalFile.Read(ISignalRecord, SIGNALFILE_FLAGS_CALIBRATE)
                
                
                indS = 0
                # TODO: Évaluer la possibilité d'utiliser GetRecordBuffer pour augmenter l'efficacité...
                record = ISignalRecord.GetRecordData()       
                
                for i in range(self.nbChannels):
                    
                    recordedSignal = array(record[indS:(indS+self.pageNbSamples[i])])
                    if noPage  == nbPages -1:
                        recordedSignal = append(recordedSignal, zeros(self.samplingRates[i]*self.pageDuration - len(recordedSignal)))

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

                    if fileType == "EDF":
                        f.write(recordedSignal.astype('<h').tostring())
                    elif fileType == "BDF": 
                        #print "writing BDF..."
                        # Writing in a string of 32-bit integers and removing every fourth byte
                        # to obtain a string of 24-bit integers
                        recordedSignal = recordedSignal.astype('<i').tostring()
                        recordedSignal = "".join([recordedSignal[noBit] for noBit in range(len(recordedSignal)) if (noBit+1)%4])
                        f.write(recordedSignal)
                    
                    #dig = np.fromstring(samples, '<i2').astype(float)
                    #phys = (dig - self.header.dig_min[i]) * self.header.gain[i] + self.header.phys_min[i]
                    #signals.append(phys)                    
                    
                    

                    indS += int(self.pageNbSamples[i])
                
                # Annotation channel
                timeDiff = ole2datetime(ISignalRecord.GetStartTime()) - self.recordingStartDateTime
                timeKeepingStr = "+" + str(timeDiff.seconds + timeDiff.days * 86400) + "\x14\x14\0"

                while(noEvent < nbEvents and len(timeKeepingStr) + len(eventStr) <= annotationFieldLength*nbByte):
                    timeKeepingStr += eventStr
                    #print eventStr
                    noEvent += 1
                    if noEvent < nbEvents :
                        eventStr = edfEventEncode(self.events[noEvent])  
                    
                f.write(timeKeepingStr +  "\0"*(annotationFieldLength*nbByte-len(timeKeepingStr)))         
            

                #print timeKeepingStr +  "\0"*(annotationFieldLength*2-len(timeKeepingStr))
                
            #ds = 0
            #for i in range(ns): 
            #     ds += self.pageNbSamples[i]*2
            #ds += annotationFieldLength*2        
    
        if noEvent < nbEvents:
            raise IOError("Not enough place in EDF Annotation signal to store all events.")
        f.closed


