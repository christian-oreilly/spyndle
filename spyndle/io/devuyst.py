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
import io  # For buffered writing
from datetime import datetime
from scipy import array, concatenate, zeros
from scipy.io import savemat            
            
from spyndle.miscellaneous import channelType
from spyndle.io import Event, EEGDBReaderBase, RecordedChannel


class DevuystEvent(Event):
    def __init__(self, groupName, name, channel, startTime, timeDuration, samplingRate):
        self.no          = None
        self.groupName   = groupName
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
 
    def getChannelLabels(self):
        return self.labels

    def getAvailableChannels(self): # deprecated
        return self.labels
 


    def importEvents(self, fname, eventName, groupName="Fuseau"):

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
            self.events.append(DevuystEvent(groupName, eventName, self.labels[0], 
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

        stageLabels = [u"Sleep stage 4", u"Sleep stage 3", u"Sleep stage 2", 
                       u"Sleep stage 1", u"Sleep stage R", u"Sleep stage W", 
                       u"Sleep stage ?"]        
        
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




    """
     Convert and save the file to the EDF format.
    """
    def saveAsEDF(self, filename, fileType = "EDF", verbose=True):
        
 
        def prepareEventStr(timeDiff, events):
            # Annotation channel            
            timeKeepingStr = "+" + str(timeDiff) + "\x14\x14\0"       
                            
            for event in events:
                timeKeepingStr += event.toEDFStr() #edfEventEncode(event)  

            return timeKeepingStr
 


        pageDuration = 30.0
        nbPages      = int(self.nbSamples/self.samplingRates[0]/pageDuration)

        pageEventStr = []
        filteredEvents = filter(lambda e: e.startTime < pageDuration, self.events)     
        pageEventStr.append(prepareEventStr(0.0, filteredEvents))
        nbEvents = len(filteredEvents)        
        
        for nopage in range(1, nbPages-1):
            filteredEvents = filter(lambda e: e.startTime <  pageDuration*(nopage+1) and 
                                              e.startTime >= pageDuration*nopage, self.events)     
            pageEventStr.append(prepareEventStr(pageDuration*nopage, filteredEvents))                                           
            nbEvents += len(filteredEvents)  
            
        filteredEvents = filter(lambda e: e.startTime >= pageDuration*(nbPages-1), self.events)     
        pageEventStr.append(prepareEventStr(pageDuration*(nbPages-1), filteredEvents))      
        nbEvents += len(filteredEvents)  
    
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
            id      = "X"
            gender  = "M"
            date    = "X"     
            fname   = "firstName"
            lname   = "lastName"
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
            writeStr = "Startdate " + "X" + " " + "X" + " " + "X" + " " + "X"  
            f.write(writeStr + (80-len(writeStr))*" ")              
              
            #8 ascii : startdate of recording (dd.mm.yy)
            f.write("01.01.01")
            
            
            #8 ascii : starttime of recording (hh.mm.ss) 
            f.write("01.01.01")

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
            
            f.write("%08d" % nbPages)  
            
            # 8 ascii : duration of a data record, in seconds
            f.write(("%8.6f" % pageDuration)[:8] )  
            
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
            
            physicalMinimum = -1000.0
            physicalMaximum =  1000.0
            digitalMinimum  = -8388608
            digitalMaximum  =  8388607
            # ns * 8 ascii : ns * physical minimum (e.g. -500 or 34)
            for i in range(ns): 
                f.write(("%8.6f" %  physicalMinimum)[:8])
            f.write(("%8.6f" %  -1)[:8]) 


            # ns * 8 ascii : ns * physical maximum (e.g. 500 or 40)
            for i in range(ns): 
                f.write(("%8.6f" %  physicalMaximum)[:8])
            f.write(("%8.6f" %  1)[:8])  


            # ns * 8 ascii : ns * digital minimum (e.g. -2048)
            for i in range(ns): 
                f.write("%08d" %  digitalMinimum)
            f.write("%08d" %  -32768)        


            # ns * 8 ascii : ns * digital maximum (e.g. 2047)
            for i in range(ns): 
                f.write("%08d" %  digitalMaximum)
            f.write("%08d" %  32767)




            # ns * 80 ascii : ns * prefiltering (e.g. HP:0.1Hz LP:75Hz)
            for i in range(ns+1): f.write(" "*80)


            if fileType == "EDF":
                nbByte = 2
            elif fileType == "BDF":
                nbByte = 3


            # ns * 8 ascii : ns * nr of samples in each data record
            nbSamplePerPage = []
            for channel in self.labels: 
                nbSamplePerPage.append(int(self.samplingRates[0]*pageDuration))
                f.write("%08d" % nbSamplePerPage[-1])
            annotationFieldLength = int(max(400, max(array([len(eventStr) for eventStr in pageEventStr]))/nbByte*1.2))            
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
                        
            for nopage in range(nbPages):
             
                for iChannel in range(self.nbChannels):

              
                    indStart = nopage*nbSamplePerPage[0]
                    indStop  = indStart + nbSamplePerPage[0]

                    if indStop > len(self.signal[iChannel]):  # If the page is incomplete (this should be the last page)
                        print nopage, nbPages, indStop, len(self.signal[iChannel])
                        recordedSignal = concatenate((self.signal[iChannel][indStart:len(self.signal[0])]), 
                                                zeros(indStop -len(self.signal[iChannel])))
                    else:
                        recordedSignal = self.signal[iChannel][indStart:indStop]

                    
                    # WRITE RECORDED SIGNAL....

                    phys_range = physicalMaximum - physicalMinimum
                    dig_range = digitalMaximum - digitalMinimum
                    gain = dig_range/phys_range   

                    recordedSignal = (recordedSignal - physicalMinimum)*gain + digitalMinimum       


                    if nbByte == 2: # EDF
                        f.write(recordedSignal.astype('<h').tostring())
                    elif nbByte == 3: #BDF
                        # Writing in a string of 32-bit integers and removing every fourth byte
                        # to obtain a string of 24-bit integers
                        recordedSignal = recordedSignal.astype('<i').tostring()
                        recordedSignal = "".join([recordedSignal[noBit] for noBit in range(len(recordedSignal)) if (noBit+1)%4])
                        f.write(recordedSignal)
                                                

                # Annotation channel            
                f.write((pageEventStr[nopage] +  "\0"*(annotationFieldLength*nbByte-len(pageEventStr[nopage]))).encode("utf8"))         

        f.closed



                        