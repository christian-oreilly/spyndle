
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
class MuscularArtifactDetector(ThresholdDetector):

    
    def __init__(self, *args, **kwargs):
        super(MuscularArtifactDetector, self).__init__(*args, **kwargs)
        
        ###############################################################################
        # Detection patameters
        ###############################################################################
    

        self.threshold = 5.75 #uV

        self.lowFreq  = 20.0 #Hz
        self.highFreq = 45.0 #Hz
           
        # The supplementary windows on each side of the detection that
        # we include in the event duration for rejection purpose
        self.rejectionPadding = 3.0 # in seconds
        
        self.averagingWindowSize = 1.0 # in seconds
        

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
        if np.mod(windowNbSample, 2) == 0 : # We need an odd number.
            windowNbSample += 1


        meanX  = self.averaging(signal, windowNbSample)
        meanX2 = self.averaging(signal**2, windowNbSample)
        sdX    = np.sqrt(meanX2 - meanX**2)
        return sdX
    
    
    def postprocessing(self, EEGsignal, channelTime, 
                                              startInds, stopInds, newEvents, fs):
    
        for startInd, stopInd, event in zip(startInds, stopInds, newEvents):
            event.startTime = event.startTime - self.rejectionPadding
            event.duration = event.duration + 2.0*self.rejectionPadding
                
    