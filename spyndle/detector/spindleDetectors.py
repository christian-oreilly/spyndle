
"""
    Implementation of spindle detectors. 

    Copyright (C) 2013  Christian O'Reilly

    For personnal, educationnal, and research purpose, this software is 
    provided under the GNU GPL (V.3) license: you can redistribute it and/or
    modify it under the terms of the version 3 of the GNU General Public 
    License as published by the Free Software Foundation.
          
    To use this software in commercial application, please contact the author.
    
    If this code is used for research purpose, the reference [1] should be
    cited in the derived publication to refer the reader to the description 
    of the methodology.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


 Author: Christian O'Reilly (christian.oreilly@umontreal.ca)
 Date  : January 8, 2013

 [1] O'Reilly, C. & Nielsen, T. Sleep spindle detection: Automation and 
     performance evaluation using fine temporal resolution, Submitted to
     Exert Systems with Applications, Augustl 2013. 

"""

import os, gc

from abc import ABCMeta
import numpy as np
import warnings
from datetime import timedelta

from scipy import concatenate, zeros, mean, sqrt, mod, where, fft
from scipy import array, arange, unique
from scipy.stats.mstats import mquantiles
from scipy.fftpack import fftfreq
from scipy.io import savemat, loadmat
from scipy.integrate import trapz


from spyndle import Filter
from spyndle import computeST
from spyndle.detector import ThresholdDetector, DetectedEvent
from spyndle.io import Event, DataManipulationProcess, \
    PSGNight, Channel, SpindleEvent, EventClass




# Class representing a detected spindle.
class DetectedSpindle(DetectedEvent):
    def __init__(self, channel, startTime, endTime):
        super(DetectedSpindle, self).__init__(channel, startTime, endTime)

        self.RMSAmp         = 0.0
        self.meanFreq       = 0.0
        self.modeFreq       = 0.0
        self.sleepStage     = ""
        self.slopeOrigin    = 0.0
        self.slope          = 0.0




    def __str__(self):
        return   "{channel:"      + str(self.channel)       \
               + ", startTime:"   + str(self.startTime())   \
               + ", endTime:"     + str(self.endTime())     \
               + ", RMSAmp:"      + str(self.RMSAmp)        \
               + ", meanFreq:"    + str(self.meanFreq)      \
               + ", sleepStage:"  + str(self.sleepStage)    \
               + ", sleepStage:"  + str(self.sleepStage)    \
               + ", slopeOrigin:" + str(self.slopeOrigin)   \
               + ", slope:"       + str(self.slope) + "}"




    def computeRMS(self, fmin=11, fmax=16, reader=None):
        raise UserWarning("The code of this function need to be recoded to use time rather than sample.")        
        
        if reader is None:
            reader = self.reader
        
        # The filters need the signal to be at least 3*nbTaps
        nbTaps = 1001        
        sampleDuration = max(nbTaps*3, self.endSample-self.startSample)        
        data = reader.read(self.channel, self.startSample, sampleDuration=sampleDuration)
        
        
        #print self.channel #, data
        signal      = data[1][0]
        fs          = data[0][0]                       # sampling rate                     
            
            
        # Defining EEG filters
        bandPassFilter = Filter(fs)
        bandPassFilter.create(low_crit_freq=fmin, high_crit_freq=fmax, 
                              order=nbTaps, btype="bandpass", ftype="FIR", 
                              useFiltFilt=True)          
             
        signal = bandPassFilter.applyFilter(signal)[0:(self.endSample-self.startSample)]                                
    
        self.RMSAmp = sqrt(mean(signal**2))
        

    
                
        
    def computeMeanFreq(self, reader, fmin=10, fmax=16):
        raise UserWarning("The code of this function need to be recoded to use time rather than sample.")        
        
        
        data = reader.read(self.channel, self.startSample, sampleDuration=self.endSample-self.startSample)

        samplingRate = data[0][0]
        signal       = data[1][0]

        if signal.size < 512:
            signal = concatenate((signal, zeros(512 - signal.size)))
        
        
        FFT = abs(fft(signal))
        freqs = fftfreq(signal.size, 1.0/samplingRate)        
        
        FFT   = FFT[(freqs >= fmin)*(freqs <= fmax)]
        freqs = freqs[(freqs >= fmin)*(freqs <= fmax)]
        
        self.meanFreq = sum(freqs*FFT)/sum(FFT)
        
        
        
    def computeTimeDuration(self, reader):
        raise DeprecationWarning("No need anymore to call this function.")




"""
 All spindle detector based on using some variable bearing information on
 spindle presence can be thought using a same general framework. This class
 implement this general framework. Threshold-based detector can be implemented
 by subclassing this classe.
"""
class ThresholdSpindleDetector(ThresholdDetector, metaclass=ABCMeta):
    def __init__(self, *args, **kwargs):
        super(ThresholdSpindleDetector, self).__init__(*args, **kwargs)
         

        self.detectClass = DetectedSpindle         
         
        ###############################################################################
        # Detection patameters
        ###############################################################################
        # Low-end of the spindle frequency band
        self.lowFreq  = 11.0
        
        # High-end of the spindle frequency band        
        self.highFreq    = 16.0
        
        # Minimal duration of a valid spindle
        self.minimalDuration = 0.5 # in seconds

        # Maximal duration of a valid spindle. Avoid aberrant spindle extraction
        # such spindle with would last for tens of seconds because of a 
        # thresholding problem. It also avoid some problems related to the fact
        # that spindles are taken as being of short duration (e.g. a script
        # stoping crashing because it has exhausted all the memory making
        # a S-transform of a spindle of 30 seconds).
        self.maximalDuration = 2.0 # in seconds
        
        self.computeRMS         = False
        self.computeFreq        = False       
        self.computeSlopeFreq   = False
        self.computeFreqMode    = False 




    @property
    def detectedSpindles(self):
        return self.detectedEvents
        
    @detectedSpindles.setter
    def detectedSpindles(self, events):
        self.detectedEvents = events
        
        

    def computeRMS(self, fmin=11, fmax=16):
        super(ThresholdSpindleDetector, self).computeRMS(fmin, fmax)



    # Used to save detected events in EEG data file. 
    def saveEvents(self, reader, eventName, 
                    eventGroupName="Spindle", fileName = None, dbSession=None, removeExisting=True):
              
        # TODO: This code should not repeat the code of TransientDetector.saveEvents(...)              
        if dbSession :
            
            # Create the data manipulation process record.
            dataManipObj = DataManipulationProcess(reprStr  = repr(self))
            dbSession.add(dataManipObj)                             

            # Create the event class record if none corresponding to this event class exist.                                       
            if dbSession.query(EventClass).filter_by(name=eventName).count() == 0:
                dbSession.add(EventClass(name=eventName))     

            # Create the PSG night record if none corresponding to this night exist. 
            if dbSession.query(PSGNight).filter_by(fileName=reader.fileName).count() == 0:
                dbSession.add(PSGNight(fileName=reader.fileName))                            

            # Create the channel record if none corresponding to this channel exist.                    
            for channel in unique([e.channel for e in self.detectedEvents]):
                if dbSession.query(Channel).filter_by(name=channel).count() == 0:
                    dbSession.add(Channel(name=channel))                    
            dbSession.flush()                          
             


        if dbSession :
            transientEvents = []
            spindleEvents   = []
          
        if removeExisting:   
            reader.events.removeType(eventName)
            
        for detectedEvent in self.detectedEvents:    
            event = Event(name          = eventName, 
                          groupName     = eventGroupName, 
                          channel       = detectedEvent.channel, 
                          startTime     = detectedEvent.startTime(),
                          timeLength    = detectedEvent.timeDuration , 
                          dateTime      = reader.getRecordingStartTime() + timedelta(seconds=detectedEvent.startTime()),
                          properties = {"RMSAmp"            :detectedEvent.RMSAmp,
                                        "meanFreq"          :detectedEvent.meanFreq,
                                        "modeFreq"          :detectedEvent.modeFreq,
                                        "slopeOrigin"       :detectedEvent.slopeOrigin,
                                        "slope"             :detectedEvent.slope,
                                        "stage"             :detectedEvent.sleepStage})     

            reader.addEvent(event)

            if dbSession :
                transient, spindle = SpindleEvent.fromEvent(event, reader.fileName, dataManipObj.no)
                transientEvents.append(transient)
                spindleEvents.append(spindle)

        if dbSession :
            dbSession.add_all(transientEvents)
            dbSession.flush()
            dbSession.add_all(spindleEvents)
            dbSession.flush()
        
        if fileName is None:
            reader.save()
        else:      
            reader.saveAs(fileName)  
            
        
        
        
        
    def computeMeanFreq(self, fmin=10, fmax=16):
        for spindle in self.detectedSpindles:        
            spindle.computeMeanFreq(self.reader, fmin, fmax)




    """
     This function could be more elegantly written as :
         
        def simpleSetStagesToSpindles(self, reader):
            for spindle in self.detectedSpindles:        
                spindle.setStage(reader)
        
      but such an implementation is 50 times slower at execution than the 
      proposed version.
    """    
    def setStagesToSpindles(self, reader):
        warnings.warn("This method has been deprecated. setStagesToEvent(...) "\
                      "should be used instead.", DeprecationWarning )              
        
        self.setStagesToEvents(reader)



    """
     Detect every spindle in the channels channelList of the file opened 
     by the reader.
    """
    def detectSpindles(self, channelList=None, reader=None, verbose=True) : 
        warnings.warn("This method has been deprecated. detectEvents(...) "\
                      "should be used instead.", DeprecationWarning )            
        self.detectEvents(channelList, reader, verbose)


        
        
        

    # Used to save detected spindle in EEG data file. 
    def saveSpindle(self, reader, eventName, 
                    eventGroupName="Spindle", fileName = None, dbSession=None):
              
        warnings.warn("This method has been deprecated. saveEvents(...) "\
                      "should be used instead.", DeprecationWarning )            
        self.saveEvents(reader, eventName, eventGroupName, fileName, dbSession)              
              
              


    # Save detected spindle in EEG data file. 
    def saveSpindleTxt(self, fileName):

        try:
            f = open(fileName, "w")    
        except IOError:     
            print("Error: The selected file could not be open.")
            exit()
        
        for spindle in self.detectedSpindles:            
            f.write(spindle.channel + ";" + str(spindle.startTime()) + ";" +  
                    str(spindle.endTime()) + ";" + str(spindle.RMSAmp) + ";" + str(spindle.meanFreq )
                      + ";" + str(spindle.timeDuration) + ";" + str(spindle.sleepStage) + "\n")             



    def loadSpindleTxt(self, fileName):
        import csv

        try:
            cr = csv.reader(open(fileName, "rb"), delimiter=';')   
        except IOError:     
            print("Error: The selected file could not be open.")
            exit()
        
        self.detectedSpindles = []
        for row in cr:
            spindle = DetectedSpindle(row[0], int(row[1]), int(row[2]))   
            self.detectedSpindles.append(spindle)
  

    def computeFreqSlope_atDetection(self, signal, fs, startSpinInds, stopSpinInds, newSpindles, channelTime):
        
         for startInd, stopInd, spindle in zip(startSpinInds, stopSpinInds, newSpindles):
            sig       = signal[startInd:stopInd]
   
            nbMinSamples = int(fs)
            N            = len(sig)
            if N < nbMinSamples:
                sig = np.concatenate((sig, np.zeros(nbMinSamples - N)))   
   
   
            X, fX = computeST(sig, fs, fmin=self.lowFreq-1, fmax=self.highFreq+1)  
            
            Y = abs(np.transpose(X))
                            
            regx = arange(N)/fs
            regy = []
            try:
                for i in range(N):
                    regy.append(trapz(fX*Y[:, i], fX)/trapz(Y[:, i], fX))  
            except: 
                print((fX.shape, X.shape, Y.shape, regx.shape, sig.shape))
                print((fs, self.lowFreq-1, self.highFreq+1, sig))
                print((channelTime[startInd:stopInd]))
                raise

            z = np.polyfit(regx, regy, deg=1)     

            spindle.slopeOrigin = z[1]
            spindle.slope       = z[0]







    def postDetectionComputation(self, EEGsignal, channelTime, 
                                              startInds, stopInds, newEvents, fs):

        # The RMS amplitude computed here is the RMS amplitude of the
        # signal filtered in the spindle band.
        if self.computeRMS :
            order = int(min(3003, len(EEGsignal)-3)/3)
            bandPassFilter = Filter(fs)
            bandPassFilter.create(low_crit_freq=self.lowFreq, 
                              high_crit_freq=self.highFreq, order=order, 
                              btype="bandpass", ftype="FIR", useFiltFilt=True)     
            filteredEEGSignal = bandPassFilter.applyFilter(EEGsignal)

            for startInd, stopInd, spindle in zip(startInds, stopInds, newEvents):
                spindle.RMSAmp= sqrt(mean(filteredEEGSignal[startInd:stopInd]**2))            
    
    
        if self.computeFreq :
            for startInd, stopInd, spindle in zip(startInds, stopInds, newEvents):
                sig       = EEGsignal[startInd:stopInd]

                if sig.size < 512:
                    sig = concatenate((sig, zeros(512 - sig.size)))
                
                FFT = abs(fft(sig))
                freqs = fftfreq(sig.size, 1.0/fs)        
                
                FFT   = FFT[(freqs >= self.lowFreq)*(freqs <= self.highFreq)]
                freqs = freqs[(freqs >= self.lowFreq)*(freqs <= self.highFreq)]
                
                spindle.meanFreq = sum(freqs*FFT)/sum(FFT)
                    
    
        if self.computeFreqMode :
            for startInd, stopInd, spindle in zip(startInds, stopInds, newEvents):
                sig       = EEGsignal[startInd:stopInd]

                if sig.size < fs*10:
                    sig = concatenate((sig, zeros(int(fs)*10 - sig.size)))
                
                FFT = abs(fft(sig))
                freqs = fftfreq(sig.size, 1.0/fs)        
                
                FFT   = FFT[(freqs >= self.lowFreq)*(freqs <= self.highFreq)]
                freqs = freqs[(freqs >= self.lowFreq)*(freqs <= self.highFreq)]
                
                spindle.modeFreq = freqs[np.argmax(FFT)]
                
        
        if self.computeSlopeFreq :
            self.computeFreqSlope_atDetection(EEGsignal, fs, startInds, 
                                              stopInds, newEvents, channelTime)

        #######################################################################




    def _detectMain_(self, *args, **kwargs):
        return super(ThresholdSpindleDetector, self)._detectMain_(*args, **kwargs)
        
        

class SpindleDetectorRMS(ThresholdSpindleDetector):
    
    
    def __init__(self, *args, **kwargs):
        super(SpindleDetectorRMS, self).__init__(*args, **kwargs)
        
        ###############################################################################
        # Detection patameters
        ###############################################################################
    
        # Quantile of spindle amplitude used to determine spindle activity
        self.threshold = 0.95

        # Width of the window used to computer the RMS amplitude.
        ### Choosing such that it always encompass at least two cycles of the 
        ### smaller frequency. We have a tradeoff between stability and accuracy...
        self.averagingWindowSize = 0.20  #0.508361 #2.0/lowFreq # In seconds

        ###############################################################################

        # We allow for a spindle activity region to contain samples having an
        # amplitude going bellow the quantileThreshold for a 
        # maximum of maxAllowableGapBellowThreshold seconds.
        self.maxAllowableGapBellowThreshold = 0.0 # in seconds
        ###############################################################################
        

    # Function performing any processing steps on the data to compute and return
    # the transformed signal that is an index of spindle presence.
    def preprocessing(self, signal, time=None, samplingRate=None, channel=None):

        # Defining EEG filters
        order = int(min(3003, len(signal)-3)/3)
        bandPassFilter = Filter(samplingRate)
        bandPassFilter.create(low_crit_freq=self.lowFreq, 
                              high_crit_freq=self.highFreq, order=order, 
                              btype="bandpass", ftype="FIR", useFiltFilt=True)          

        # filtering can take a lot of memory. By making sure that the 
        # garbage collector as passed just before the filtering, we
        # increase our chances to avoid a MemoryError  
        gc.collect()                      
        signal = bandPassFilter.applyFilter(signal)     

        ################################# RMS COMPUTATION #####################
        windowNbSample = int(round(self.averagingWindowSize*samplingRate))
        if mod(windowNbSample, 2) == 0 : # We need an odd number.
            windowNbSample += 1

        # For selecting samples using a quantile-based thresholds, using abs(X)
        # or X**2 to rectify the X signal will give exactly the same result
        # since X**2 eqauls abs(X)**2 (for real numbers) and  the transformation
        # from abs(X) to abs(X)**2 is monotonically increasing, meaning that 
        # rank based statistics (such as quatiles) will give exactly the same
        # result. We use the numpy implementation of abs() because it is the
        # faster alternative.
        return np.sqrt(self.averaging(signal**2.0, windowNbSample))
        #return self.averaging(np.abs(signal), windowNbSample)
    
    
    
    def getEffectiveThreshold(self, signal):
        return mquantiles(signal, self.threshold)[0]








# Mean absolute
class SpindleDetectorMA(ThresholdSpindleDetector):
    
    
    def __init__(self, *args, **kwargs):
        super(SpindleDetectorMA, self).__init__(*args, **kwargs)
        
        ###############################################################################
        # Detection patameters
        ###############################################################################
    
        # Quantile of spindle amplitude used to determine spindle activity
        self.threshold = 0.92

        # Width of the window used to computer the RMS amplitude.
        ### Choosing such that it always encompass at least two cycles of the 
        ### smaller frequency. We have a tradeoff between stability and accuracy...
        self.averagingWindowSize = 0.20  #0.508361 #2.0/lowFreq # In seconds

        ###############################################################################

        # We allow for a spindle activity region to contain samples having an
        # amplitude going bellow the quantileThreshold for a 
        # maximum of maxAllowableGapBellowThreshold seconds.
        self.maxAllowableGapBellowThreshold = 0.0 # in seconds
        ###############################################################################
        

    # Function performing any processing steps on the data to compute and return
    # the transformed signal that is an index of spindle presence.
    def preprocessing(self, signal, time=None, samplingRate=None, channel=None):

        # Defining EEG filters
        order = int(min(3003, len(signal)-3)/3)
        bandPassFilter = Filter(samplingRate)
        bandPassFilter.create(low_crit_freq=self.lowFreq, 
                              high_crit_freq=self.highFreq, order=order, 
                              btype="bandpass", ftype="FIR", useFiltFilt=True)          

        # filtering can take a lot of memory. By making sure that the 
        # garbage collector as passed just before the filtering, we
        # increase our chances to avoid a MemoryError  
        gc.collect()                      
        signal = bandPassFilter.applyFilter(signal)     

        ################################# RMS COMPUTATION #####################
        windowNbSample = int(round(self.averagingWindowSize*samplingRate))
        if mod(windowNbSample, 2) == 0 : # We need an odd number.
            windowNbSample += 1


        return self.averaging(np.abs(signal), windowNbSample)
    
    
    
    def getEffectiveThreshold(self, signal):
        return mquantiles(signal, self.threshold)[0]






# For backward compatibility.
SpindleDetectorAmp = SpindleDetectorRMS

 
 
 
class SpindleDetectorTeager(ThresholdSpindleDetector):
    
    
    def __init__(self, *args, **kwargs):
        super(SpindleDetectorTeager, self).__init__(*args, **kwargs)
        
        ###############################################################################
        # Detection patameters
        ###############################################################################
        self.threshold = 0.6

        # We allow for a spindle activity region to contain samples having an
        # amplitude going bellow the quantileThreshold for a 
        # maximum of maxAllowableGapBellowThreshold seconds.
        self.maxAllowableGapBellowThreshold = 0.0 # in seconds
        ###############################################################################


    # Function performing any processing steps on the data to compute and return
    # the transformed signal that is an index of spindle presence.
    def preprocessing(self, signal, time=None, samplingRate=None, channel=None):

        # Defining EEG filters
        order = int(min(3003, len(signal)-3)/3)
        bandPassFilter = Filter(samplingRate)
        bandPassFilter.create(low_crit_freq=self.lowFreq, high_crit_freq=self.highFreq, 
                              order=order, btype="bandpass", ftype="FIR", useFiltFilt=True)          
    
        """
         Splitting the process in chuck of 1M elements to avoid memory errors.
        """
        NMAX = 1000000.0
        N = len(signal)
        output = [[0]]
        for i in range(int(N/NMAX)+1):
            block = signal[int(max(i*NMAX-1, 0)):int(min(N, (i+1)*NMAX+1))]
                   
            block = bandPassFilter.applyFilter(block)                      
    
            ########################## Computing Teager operator ######################
            output.append(block[1:-1]**2 - block[:-2]*block[2:])
        
        output.append([0])
        return  concatenate(output)  
    
      
    def getEffectiveThreshold(self, signal):
        return mean(signal)*self.threshold


 
 
class SpindleDetectorSigma(ThresholdSpindleDetector):
    
    
    def __init__(self, *args, **kwargs):
        super(SpindleDetectorSigma, self).__init__(*args, **kwargs)        

        ###############################################################################
        # Detection patameters
        ###############################################################################
        self.threshold = 4.5

        # We allow for a spindle activity region to contain samples having an
        # amplitude going bellow the quantileThreshold for a 
        # maximum of maxAllowableGapBellowThreshold seconds.
        self.maxAllowableGapBellowThreshold = 0.1 # in seconds

        # Duration of the window used to compute S-transform
        # (in seconds)
        self.computationWindow = 4.2   

        # Overlapping of two consecutive windows used to compute the S-transform
        # (in seconds)
        self.windowOverlapping = 0.2

        ###############################################################################



    # Function performing any processing steps on the data to compute and return
    # the transformed signal that is an index of spindle presence.
    def preprocessing(self, signal, time=None, samplingRate=None, channel=None):
    
        #fileName = self.reader.fileName + "_sigmaIndex_" + channel + "_" + str(time[0])  + ".mat"
    
        N = len(signal)
        """
        if os.path.exists(fileName):
            print "Using saved sigma index..."    
            self.sigmaIndex = loadmat(fileName)["sigma"]
            self.sigmaIndex = self.sigmaIndex.reshape(self.sigmaIndex.size)
            assert(len(signal) == len(self.sigmaIndex))
        else:
        """
        self.sigmaIndex = zeros(N) 
        
        nbPad = int(self.windowOverlapping/2.0*samplingRate)
        nbWin = int(self.computationWindow*samplingRate) - 2*nbPad

        nbIter = int(N/nbWin)
        for i in range(nbIter):
            if mod(i, 1000) == 0: print((i, nbIter))
            
            if i == 0 :            # First iteration
                indexes = arange(nbPad + nbWin) 
            elif i == nbIter-1 :   # Last iteration
                indexes = arange(i*nbWin-nbPad, N)
            else:                       # Other iterations
                indexes = arange(i*nbWin-nbPad, i*nbWin + nbWin+nbPad)

            X, fX = computeST(signal[indexes], samplingRate, fmin=4.0, fmax=40.0)  
            
            if i == 0 :            # First iteration
                indexesNoPad = arange(nbWin) 
            elif i == nbIter-1 :   # Last iteration
                indexesNoPad = arange(nbPad, nbPad + N-i*nbWin)
            else:                       # Other iterations
                indexesNoPad = arange(nbPad, nbPad + nbWin)
                
            X       = abs(X[indexesNoPad])
            indexes = indexes[indexesNoPad]

            maxalpha  = X[:, (fX >= 7.5)*(fX <= 10.0)].max(axis=1)
            maxsigma  = X[:, (fX >= self.lowFreq)*(fX <= self.highFreq)].max(axis=1)
            meanlow   = X[:, (fX >= 4.0)*(fX <= 10.0)].mean(axis=1) 
            meanhigh  = X[:, (fX >= 20.0)*(fX <= 40.0)].mean(axis=1) 

            sigmaTMP = array(2*maxsigma/(meanlow + meanhigh))
            self.sigmaIndex[indexes[where(maxalpha <= maxsigma)[0]]] = sigmaTMP[where(maxalpha <= maxsigma)[0]]
        
        #    savemat(fileName, {"sigma":self.sigmaIndex})
        #    assert(len(signal) == len(self.sigmaIndex))
        return self.sigmaIndex
              

 
 
class SpindleDetectorRSP(ThresholdSpindleDetector):
    
    
    
    def __init__(self, *args, **kwargs):
        super(SpindleDetectorRSP, self).__init__(*args, **kwargs)        
        
        ###############################################################################
        # Detection patameters
        ###############################################################################
        self.threshold = 0.22

        # We allow for a spindle activity region to contain samples having an
        # amplitude going bellow the quantileThreshold for a 
        # maximum of maxAllowableGapBellowThreshold seconds.
        self.maxAllowableGapBellowThreshold = 0.0 # in seconds

        ###############################################################################


    # Detect every spindle in the channels channelList of the file opened by the reader.
    def preprocessing(self, signal, time=None, samplingRate=None, channel=None):

        #fileName = self.reader.getFileName() + "_RSP_" + channel + "_" + str(time[0])  + ".mat"

        N = len(signal)
        """
        if os.path.exists(fileName):
            print "Using saved RSP..."    
            self.RSP = loadmat(fileName)["RSP"]
            self.RSP = self.RSP.reshape(self.RSP.size)            
            assert(len(signal) == len(self.RSP))
        else: 
        """
        self.RSP = zeros(N) 

        nbPad = int(0.1*samplingRate)
        nbWin = int(4.0*samplingRate)

        nbIter = int(N/nbWin)
        for i in range(nbIter):
            if mod(i, 1000) == 0: print((i, nbIter))
            
            if i == 0 :            # First iteration
                indexes = arange(nbPad + nbWin) 
            elif i == nbIter-1 :   # Last iteration
                indexes = arange(i*nbWin-nbPad, N)
            else:                       # Other iterations
                indexes = arange(i*nbWin-nbPad, i*nbWin + nbWin+nbPad)


            #if any(stageIndicator[indexes]):    
            X, fX = computeST(signal[indexes], samplingRate, fmin=0.5, fmax=40.0)  
            
            if i == 0 :            # First iteration
                indexesNoPad = arange(nbWin) 
            elif i == nbIter-1 :   # Last iteration
                indexesNoPad = arange(nbPad, nbPad + N-i*nbWin)
            else:                       # Other iterations
                indexesNoPad = arange(nbPad, nbPad + nbWin)
                
            X       = abs(X[indexesNoPad])
            indexes = indexes[indexesNoPad]


            spindleBand = (fX >= self.lowFreq)*(fX <= self.highFreq)
            self.RSP[indexes] = trapz(X[:, spindleBand], fX[spindleBand], axis=1)/trapz(X, fX, axis=1)
            #else:
            #    self.RSP[indexes] = 0.0
    
        #    assert(len(signal) == len(self.RSP))
        #    savemat(fileName, {"RSP":self.RSP})
        return self.RSP
