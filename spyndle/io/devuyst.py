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

import os
from datetime import datetime
from scipy import array
from scipy.io import savemat            
            
from spyndle.filters import channelType
from spyndle.io import Event, EEGDBReaderBase, RecordedChannel


class DevuystEvent(Event):
    def __init__(self, groupeName, name, channel, startTime, timeDuration, samplingRate):
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






# TODO: Abstract methods in EEGDBReaderBase must be implemented before we can
# put back the inheritence.
class DevuystReader(): #(EEGDBReaderBase):
    def __init__(self, fname, samplingRate):
        #EEGDBReaderBase.__init__(self)
 
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




    def readWithTime(self, signalNames, startTime, timeDuration):
        return self.read(signalNames, int(startTime*self.baseFreq), timeDuration=timeDuration)



    def readCompletePickle(self, signalNames):

        # No pickling necessary for these data...
        return self.readComplete(signalNames)
               


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

            print data

            for signalName in data:
                print "Pickling data of " + signalName + " for next time..."
                savemat(fileName + signalName + ".mat", data[signalName])



    # As read, but read the complete duration of the signals.
    def readComplete(self, signalNames):
        if self.labels[0] in signalNames:
            
            returnData = {}
  
            for i in range(len(self.labels)):
                returnData[self.labels[i]]                = RecordedChannel()                
                returnData[self.labels[i]].signal         = self.signal[i]
                returnData[self.labels[i]].samplingRate   = self.samplingRates[i]
                returnData[self.labels[i]].type           = self.channelType[i]
                returnData[self.labels[i]].startTime      = [datetime.now().strftime("%a, %d %b %Y %H:%M:%S +0000")]

            return returnData
        else:
            return {}



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



                        