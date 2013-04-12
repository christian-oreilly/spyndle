# -*- coding: utf-8 -*-

"""
    Code assessing transient event propogation through an array of sensors
    using cross-correlation of S-transform of the signal captured by
    the different sensors.

    Copyright (C) 2012-2013  Christian O'Reilly

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.


    For personnal, educationnal, and research purpose, this software is 
    provided under the GNU GPL (V.3) license: you can redistribute it and/or
    modify it under the terms of the version 3 of the GNU General Public 
    License as published by the Free Software Foundation.
          
    To use this software in commercial application, please contact the author. 
    If used for research purpose, the reference [1] or references [2-3] should 
    be cited in the derived publication to refere the reader to the description 
    of the methodology.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


 Author: Christian O'Reilly (christian.oreilly@umontreal.ca)
 Date  : July 11, 2012

 [1] O’Reilly, C., & Nielsen, T. (2013). Assessing the propagation of EEG 
     transient activity, Proceedings of the 8th International Workshop on 
     Systems, Signal Processing and their Applications, Algiers, Algeria, 
     12-15 May 2013.
 [2] O'Reilly, C. & Nielsen, T. Assessing EEG sleep spindle propagation. 
     Part 1: Theory and proposed methodology, Submitted to Journal of 
     Neuroscience Methods, april 2013.
 [3] O'Reilly, C. & Nielsen, T. Assessing EEG sleep spindle propagation. 
     Part 2: Experimental characterization, Submitted to Journal of 
     Neuroscience Methods, april 2013.     

"""


from scipy import where, zeros, transpose, arange, shape, diff, sign
from scipy.integrate import trapz

from ..utils import diff2
#from ..filters import Filter
from ..STransform import computeMST



# Computation of the 2D cross-corellation of S-Transform spectra as described 
# in [1].
#
# Parameters:
#    readerClass    : Class of EEG reader. As of now, only HarmonieReader has 
#                     been implementer to be able to read .SIG file generated
#                     from the Harmonie software. To use this software, users 
#                     will most probably have to implement their own reader for 
#                     EEG file format they are using folowwing the example given 
#                     by the HarmonieReader class. If help is needed, contact
#                     the author (christian.oreilly@umontreal.ca). 
#    nightLst       : File name of the EEG data.
#    dataPath       : Repertory containing the EEG data listed in nightLst.
#    resPath        : Repertory where the results will be saved.
#    eventName      : Name used to label the "spindle" events.
#    channelLst     : List of the recorded channels.
#    delta          : delta_w in [1]. The observation windows will be 
#                     slid from -delta to +delta. See [1]. (in seconds)
#                     Default value: 0.5
#    beforePad      : Padding of the observation windows before the beginning 
#                     of the spindle. Labeled delta_{sp} in [1]. (in seconds)
#                     Default value: 0.5
#    afterPad       : Padding of the observation windows after the ending 
#                     of the spindle. Labeled delta_{sp} in [1]. (in seconds)
#                     Default value: 0.5
#    offset         : Should be 0.0 for computing the delays. A value of 5.0
#                     (in seconds)  can be used to compute false positive
#                     propagation. See [1].
#                     Default value: 0.0
#    eventProperties: List of string indicating the event properties that are 
#                     to be save together with the propagation delay 
#                     information.
#                     Default value: []
#   resFilesPrefix  : Prefix that is used to name the output files. 
#                     Default value: "spindleDelays"
#                     
#
def computeXCST(readerClass, nightLst, dataPath, resPath, eventName, channelLst, 
                 delta =0.5, beforePad=0.5, afterPad=0.5, offset=0.0, eventProperties = [],
                 resFilesPrefix="spindleDelays") :  
    for night in nightLst:
        try:
            print 'Loading data of the file ' + night + '...'
            reader = readerClass(dataPath + night)    
            
            ###################################################################
            # Patch used because data recorder with previous versions of Harmonie
            # inconsistent time stamp due to discontinuous recording. These lines
            # should not be necessary with another reader.
            print "Resynchronization of the events..."
            reader.resyncEvents()    
            ###################################################################
            
        except IOError:     
            print "Error: The selected file could not be open."
            exit()        
        
        
        
        for refChannel in channelLst:        
            try:
                f = open(resPath + resFilesPrefix + "_" + night + "_" + refChannel + ".txt", "w")    
            except IOError:     
                print "Error: The selected file could not be open."
                exit()
            
            events = filter(lambda e: e.name == eventName and e.channel == refChannel, reader.events)
            for i, event in zip(range(len(events)), events):
                deltaSample = int(delta*reader.baseFreq)        
                delta       = deltaSample/reader.baseFreq    
                
                start = event.startSample -  int(beforePad*reader.baseFreq)       
                duration = event.timeLength + beforePad + afterPad
            
                signalsDataCmp = reader.read(channelLst, start - deltaSample + int(offset*reader.baseFreq), duration + 2.0*delta)  
                signalsDataRef = reader.read([refChannel], start - deltaSample, duration + 2.0*delta)  
            
                X  = []        
                    
                fs          = signalsDataRef[0][0]       # sampling rate   
                chanType    = signalsDataRef[2][0]                     
    
                #highPassFilter = Filter(fs)
                #highPassFilter.create(chanType, low_crit_freq=9.0, 
                #                      high_crit_freq=None, order=4, btype="highpass", 
                #                      ftype="butter", useFiltFilt=True)                    
            
    
                #lowPassFilter = Filter(fs)
                #lowPassFilter.create(chanType, low_crit_freq=None, 
                #                      high_crit_freq=16.0, order=4, btype="lowpass",
                #                      ftype="butter", useFiltFilt=True)                  
                    
                # Spectra computation
                for i2 in range(len(channelLst)+1) :
                    if i2 == 0 :
                        signal      = signalsDataRef[1][0]  
                    else:
                        signal      = signalsDataCmp[1][i2-1]               
                    
                             
                    #signal = highPassFilter.applyFilter(signal)                       
                    #signal = lowPassFilter.applyFilter(signal)                
                    
        
                    Xtmp, fXtmp = computeMST(signal, fs, m=0.0, k=1.0, fmin=11.0, fmax=16.0)  
                            
                            
                    Xtmp = abs(transpose(Xtmp))
                        
                    if i2 == 0:    
                        Xtmp = Xtmp[:, deltaSample:-deltaSample]
           
                    
                    X.append(Xtmp)
            
             
            
                refShape = list(shape(X[0]))      
                
                maxCrosscor = zeros(len(channelLst))
                maxDeltay   = zeros(len(channelLst))
                
                crosscor    = zeros(2*deltaSample)
                time        = (arange(2*deltaSample)-deltaSample)/fs
                refSelfCor  = trapz(trapz(X[0]*X[0])) 
                
    
                for i2 in range(1,len(channelLst)+1) :
                    if channelLst[i2-1] == refChannel:
                        maxDeltay[i2-1]   = 0       
                        maxCrosscor[i2-1] = 1.0              
                    else:                
                        XCmpSquare = X[i2]*X[i2]
                        XCmpSquare_colTrapz = trapz(XCmpSquare, axis=0)        
                        for offset in range(2*deltaSample):
                            XCmp             = X[i2][:, offset:(offset+refShape[1])]
                            CmpSelfCor       = trapz(XCmpSquare_colTrapz[offset:(offset+refShape[1])])
                            crosscor[offset] = trapz(trapz(XCmp*X[0]))/max(CmpSelfCor, refSelfCor)
    
                        ind = where(crosscor == max(crosscor))[0][0]
                        
                        diffCrosscor = diff2(arange(2*deltaSample), crosscor)                
                        
                        if ind == 0 or ind == 2*deltaSample-1:
                            maxDeltay[i2-1]   = time[ind]       
                            maxCrosscor[i2-1] = crosscor[ind]                      
                        else:
                            if ind-2 > 0 :
                                indm = ind-2
                            else:
                                indm = 0
                                
                            if ind+2 < 2*deltaSample-1 :
                                indp = ind+2
                            else:
                                indp = 2*deltaSample-1
        
                            pindinc = where(diff(sign(diffCrosscor[indm:(indp+1)])) != 0)[0][0]                
                            indFloor = indm +pindinc 
                            indFrac  = (diffCrosscor[indFloor])/(diffCrosscor[indFloor] - diffCrosscor[indFloor+1])
                          
        
                            maxDeltay[i2-1]   = time[indFloor]      + indFrac*(time[indFloor+1]-time[indFloor])       
                            maxCrosscor[i2-1] = crosscor[indFloor]  + indFrac*(crosscor[indFloor+1]-crosscor[indFloor])       
                            
                        
    
            
                result = str(i) + ";" + str(event.dateTime) + ";" + str(event.timeLength) + ";" 
                for i2 in range(len(channelLst)):
                    result += str(maxCrosscor[i2]) + ";"  +  str(maxDeltay[i2]) + ";"  
                    
                # Depending on the available properties 
                for eventProperty in eventProperties:
                    result += str(event.properties[eventProperty]) + ";" 

                print night, refChannel, i

                f.write(result[:-1] + '\n') # Write a string to a file
                
        f.close()
                