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
        
                        
            