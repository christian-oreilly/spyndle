
"""
    Implementation of detectors for muscular artifact. 

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
 Date  : February 15, 2015
"""

import gc

import numpy as np

from spyndle import Filter
from spyndle.detector import ThresholdDetector
#from spyndle.io import Event, DataManipulationProcess, \
#    PSGNight, Channel, SpindleEvent, EventClass





"""
 All event detectors using some variable bearing information on
 the presence of events to be detected can be thought using a same general 
 framework. This class implement this general framework. Threshold-based 
 detector can be implemented by subclassing this classe.
"""
class AlphaArtifactDetector(ThresholdDetector):

    
    def __init__(self, *args, **kwargs):
        super(AlphaArtifactDetector, self).__init__(*args, **kwargs)
        
        ###############################################################################
        # Detection patameters
        ###############################################################################
    

        self.threshold = 1.1 # alpha/delta ratio

        self.deltaLowFreq  = 1.0 #Hz
        self.deltaHighFreq = 4.0 #Hz

        self.alphaLowFreq  = 8.0 #Hz
        self.alphaHighFreq = 11.0 #Hz
           
        # The supplementary windows on each side of the detection that
        # we include in the event duration for rejection purpose
        self.rejectionPadding = 5.0 # in seconds
        
        self.alphaAveragingWindowSize = 15.0 # in seconds
        self.deltaAveragingWindowSize = 15.0 # in seconds
        

    # Function performing any processing steps on the data to compute and return
    # the transformed signal that is an index of spindle presence.
    def preprocessing(self, signal, time=None, samplingRate=None, channel=None):
        
        def hanningSquared(N):
            return np.hanning(N)**2

        # Defining EEG filters
        order = int(min(3003, len(signal)-3)/3)
        alphaBandPassFilter = Filter(samplingRate)
        alphaBandPassFilter.create(low_crit_freq=self.alphaLowFreq, 
                              high_crit_freq=self.alphaHighFreq, order=order, 
                              btype="bandpass", ftype="FIR", useFiltFilt=True)   
                              
        deltaBandPassFilter = Filter(samplingRate)
        deltaBandPassFilter.create(low_crit_freq=self.deltaLowFreq, 
                              high_crit_freq=self.deltaHighFreq, order=order, 
                              btype="bandpass", ftype="FIR", useFiltFilt=True)          

        # filtering can take a lot of memory. By making sure that the 
        # garbage collector as passed just before the filtering, we
        # increase our chances to avoid a MemoryError  
        gc.collect()                      
        alphaSignal = alphaBandPassFilter.applyFilter(signal)   
        deltaSignal = deltaBandPassFilter.applyFilter(signal)     

        ################################# RMS COMPUTATION #####################
        alphaWindowNbSample = int(round(self.alphaAveragingWindowSize*samplingRate))
        if np.mod(alphaWindowNbSample, 2) == 0 : # We need an odd number.
            alphaWindowNbSample += 1

        deltaWindowNbSample = int(round(self.deltaAveragingWindowSize*samplingRate))
        if np.mod(deltaWindowNbSample, 2) == 0 : # We need an odd number.
            deltaWindowNbSample += 1


        deltaRMS  = np.sqrt(self.averaging(deltaSignal**2, deltaWindowNbSample, hanningSquared))
        alphaRMS  = np.sqrt(self.averaging(alphaSignal**2, alphaWindowNbSample, hanningSquared))
        return alphaRMS/deltaRMS
    
    
    def postprocessing(self, EEGsignal, channelTime, 
                                              startInds, stopInds, newEvents, fs):
    
        for startInd, stopInd, event in zip(startInds, stopInds, newEvents):
            event.startTime = event.startTime - self.rejectionPadding
            event.duration = event.duration + 2.0*self.rejectionPadding
                
    