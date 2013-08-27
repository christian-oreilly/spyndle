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

import os
from scipy import where, zeros, transpose, arange, shape, diff, sign, argmax, sqrt
from scipy.integrate import trapz
from scipy import fft, ifft
from scipy.signal import correlate2d
import numpy as np


from spyndle import diff2
#from ..filters import Filter
from spyndle import computeMST



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
#    fileName       : File name of the EEG data.
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
#    artifactPad    : Used to pad the ST to compute it on a larger window
#                     that is afterward trucated such as to eliminate border
#                     artifacts.
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
#   intFct          : Allows to choose the function used in integration. By
#                     default nympy.sum is used fut numpy.trapz could be used
#                     for more precision (at cost of a much slower execution).
#                     
#
def computeXCST(readerClass, fileName, resPath, eventName, channelLst, 
                 delta =0.5, beforePad=0.5, afterPad=0.5, 
                 artifactPad=0.25,
                 offset=0.0, eventProperties = [],
                 resFilesPrefix="spindleDelays", similIndexType="inf",
                 intFct=np.sum) :  

    night = os.path.basename(fileName)

    try:
        print 'Loading data of the file ' + night + '...'
        reader = readerClass(fileName)    
        
    except IOError:     
        print "Error : computeXCST : ", fileName, "could not be opened."
        raise      
    
    
    
    for refChannel in channelLst:        
        with open(resPath + resFilesPrefix + "_" + night + "_" + refChannel + ".txt", "w") as f:   
            
            print "Processing propagation:", night, refChannel   
            
            events = filter(lambda e: e.name == eventName and e.channel == refChannel, reader.events)
            
            printHeader = True
            artifactPadSample = int(artifactPad*reader.getChannelFreq(refChannel))  
            for i, event in zip(range(len(events)), events):
         
                #### These are reinitialized at each event because their balue
                # can be altered during the computation og the XCST for
                # a given event
                deltaSampleStart = int(delta*reader.getChannelFreq(refChannel))        
                deltaStart       = deltaSampleStart/reader.getChannelFreq(refChannel)    
                deltaSampleEnd   = deltaSampleStart    
                deltaEnd         = deltaStart
                #### 
                    
                startTime = event.startTime -  beforePad  - artifactPad    
                duration  = event.timeLength + beforePad + afterPad + 2*artifactPad
                
                

                """
                 The signal at startTime - delta might not be available and the 
                 reader will start reading at the next available sample following 
                 the startTime requested. If this happens we may miss some samples 
                 in this signal. We remove these samples in both synchronous and
                 asynchronous signals (for them to have the same size) by reducing  
                 the delta accordingly.
                """
                refStartTime = reader.getNextSampleStartTime(refChannel, startTime - deltaStart)
                deltaStartRef = int((refStartTime - (startTime - deltaStart ))*reader.getChannelFreq(refChannel))
                    
                if deltaStartRef > 0 : 
                    # In case that the signal is not sufficiently long signal 
                    # before the spindle to compute delays that are at least half
                    # of the limit asked by the parameter delta, we reject this spindle
                    if deltaStartRef >  deltaSampleStart/2:      
                        # The condition of this assert should never be met if 
                        # there is no bug in spindle detection and EEG data readers.
                        assert(startTime - deltaStart <= event.startTime)
                        continue
                    
                    deltaSampleStart -= deltaStartRef
                    deltaStart        = deltaSampleStart/reader.getChannelFreq(refChannel)                  


                signalsDataRef = reader.read([refChannel], startTime - deltaStart         , duration + deltaStart + deltaEnd)                  
                
                
                """
                 The signal at startTime - delta + offset might not be 
                 available and the reader will start reading at the
                 next available sample following the startTime requested. If this
                 happens for the offset signal, just adjust the offset accordingly.
                """
                cmpStartTime = reader.getNextSampleStartTime(refChannel, startTime - deltaStart + offset)
                deltaStartCmp = int((cmpStartTime - (startTime - deltaStart + offset))*reader.getChannelFreq(refChannel))

                if deltaStartCmp > 0 : 
                    signalsDataCmp = reader.read(channelLst, cmpStartTime                    , duration + deltaStart + deltaEnd)  
                else:
                    signalsDataCmp = reader.read(channelLst, startTime - deltaStart + offset , duration + deltaStart + deltaEnd)  
                    
                
                ##############

                #print signalsDataCmp[refChannel].signal.shape, signalsDataRef[refChannel].signal.shape

                #################################################################
                # Because of numerical approximation, the signal length between the
                # synchronous and the asychronous comparisons can be different 
                # of 1 sample. We need to fix them to the same size to avoid
                # comparison problems.
                if offset != 0.0:
                    if len(signalsDataCmp[refChannel].signal) > len(signalsDataRef[refChannel].signal):
                        for channel in signalsDataCmp:
                            signalsDataCmp[channel].signal = signalsDataCmp[channel].signal[:-1]
                    elif len(signalsDataCmp[refChannel].signal) < len(signalsDataRef[refChannel].signal):
                        signalsDataRef[refChannel].signal = signalsDataRef[refChannel].signal[:-1]
                ##############


                fs = signalsDataCmp[refChannel].samplingRate                   
    
                # Spectra computation
                X, fXRef = computeMST(signalsDataRef[refChannel].signal, fs, m=0.0, k=1.0, fmin=11.0, fmax=16.0)     
                X = abs(transpose(X))                
                try:
                    #Y = zeros(X.shape)
                    #Y[:, (deltaSampleStart+artifactPadSample):-(deltaSampleEnd+artifactPadSample)] = X[:, (deltaSampleStart+artifactPadSample):-(deltaSampleEnd+artifactPadSample)]
                    X = {"Ref" :X[:, (deltaSampleStart+artifactPadSample):-(deltaSampleEnd+artifactPadSample)]}
                    #X["Ref2"] = Y
                except IndexError:
                    print deltaSampleStart, deltaSampleEnd, X.shape, artifactPadSample
                    raise               
                
                for testChannel in channelLst :      
                    Xtmp, fXtmp = computeMST(signalsDataCmp[testChannel].signal, fs, m=0.0, k=1.0, fmin=11.0, fmax=16.0)         
                    X[testChannel] = abs(transpose(Xtmp))[:, artifactPadSample:-artifactPadSample]              
            
                refShape = list(shape(X["Ref"]))      
                
                maxCrosscor = {}
                maxDeltay   = {}
                
                unnormCrossCorr = zeros(deltaSampleStart + deltaSampleEnd)
                time        = (arange(deltaSampleStart + deltaSampleEnd)-deltaSampleStart)/fs
                refSelfCor  = intFct(intFct(X["Ref"]*X["Ref"])) 
                indOffsets = arange(deltaSampleStart + deltaSampleEnd)
                

                for testChannel in channelLst :
                    if testChannel == refChannel:
                        maxDeltay[testChannel]   = 0       
                        maxCrosscor[testChannel] = 1.0              
                    else:                
                        XCmpSquare = X[testChannel]*X[testChannel]
                        XCmpSquare_colTrapz = intFct(XCmpSquare, axis=0)     
                        
                        
                        """
                        equivalent to:
                        CmpSelfCor  = zeros(deltaSampleStart + deltaSampleEnd)                        
                        for indOffset in indOffsets:
                            CmpSelfCor[indOffset]       = sum(XCmpSquare_colTrapz[indOffset:(indOffset+refShape[1])])                         
                        """
                        cums = np.cumsum(XCmpSquare_colTrapz)
                        CmpSelfCor = np.concatenate(( [cums[refShape[1]-1]], 
                                                    cums[indOffsets[1:] + refShape[1]-1] - cums[indOffsets[:-1]]))

    
                                          
                        # This...                        
                        #unnormCrossCorr = correlate2d(X[testChannel], X["Ref"], mode="valid")[0][indOffsets]   

                        # and this...
                        #K = X[testChannel].shape[1] 
                        #H = X["Ref"].shape[0]
                        #N = X["Ref"].shape[1]   
                        #y = np.concatenate((X["Ref"], zeros((H, K-1))), axis=1)           
                        #FFTy = np.conjugate(np.fft.fft(y, axis=1)) 
                        #x = np.concatenate((zeros((H, N-1)), X[testChannel]), axis=1)                        
                        #unnormCrossCorr = np.real(np.fft.ifft(np.fft.fft(x, axis=1)*FFTy, axis=1)[:, (N-1):(N-1 + deltaSampleStart + deltaSampleEnd)]) 
                        #unnormCrossCorr = np.sum(unnormCrossCorr, axis=0)                        


                        # are slower than that:
                        for indOffset in indOffsets:
                            XCmp                        = X[testChannel][:, indOffset:(indOffset+refShape[1])]
                            unnormCrossCorr[indOffset]  = intFct(intFct(XCmp*X["Ref"]))            


                        if similIndexType == "inf":
                            den = np.maximum(CmpSelfCor, refSelfCor)
                        elif similIndexType == "euclidean":
                            den = sqrt(CmpSelfCor**2 + refSelfCor**2)                                
                        else: 
                            print "Error. The similarity index type " + similIndexType + "is unknown."
                            raise TypeError
                            
                        crosscor = unnormCrossCorr/den
    
                        ind          = argmax(crosscor)                    
                        diffCrosscor = diff2(arange(deltaSampleStart + deltaSampleEnd), crosscor)                
        
                        # Interpolation to find the position of the true maximum
                        # which normally falls between samples
                        maxInd = len(crosscor)-1
                        if ind == 0 or ind == maxInd:
                            maxDeltay[testChannel]   = time[ind]       
                            maxCrosscor[testChannel] = crosscor[ind]                      
                        else:
                            if ind-2 > 0 :
                                indm = ind-2
                            else:
                                indm = 0
                                
                            if ind+2 < maxInd :
                                indp = ind+2
                            else:
                                indp = maxInd
        
                            pindinc = where(diff(sign(diffCrosscor[indm:(indp+1)])) != 0)[0][0]                
                            indFloor = indm +pindinc 
                            indCeil  = indFloor + 1
                            indFrac  = diffCrosscor[indFloor]/(diffCrosscor[indFloor] - diffCrosscor[indCeil])
                          
                            maxDeltay[testChannel]   = time[indFloor]      + indFrac*(time[indCeil]-time[indFloor])       
                            maxCrosscor[testChannel] = crosscor[indFloor]  + indFrac*(crosscor[indCeil]-crosscor[indFloor])       
                        
                        
                        """
                        import pylab
                        import numpy
                        pylab.figure()
                        
                        vmin= min(numpy.min(X["Ref"]), numpy.min(X[testChannel])) 
                        vmax= max(numpy.max(X["Ref"]), numpy.max(X[testChannel]))   
                            
                        pylab.subplot(311)
                        T, F = numpy.meshgrid(numpy.arange(X["Ref"].shape[1])/fs+maxDeltay[testChannel], fXRef)
                        pylab.pcolor(T, F, X["Ref"], vmin=vmin, vmax=vmax)
                        pylab.xlim([-0.5,X[testChannel].shape[1]/fs-0.5])
                                 
                        pylab.subplot(312)
                        T, F = numpy.meshgrid(numpy.arange(X[testChannel].shape[1])/fs-0.5, fXRef)
                        pylab.pcolor(T, F, X[testChannel], vmin=vmin, vmax=vmax)
                        pylab.xlim([-0.5,X[testChannel].shape[1]/fs-0.5])
                        
                        pylab.subplot(313)
                        pylab.plot(time, crosscor)
                        pylab.plot(maxDeltay[testChannel], maxCrosscor[testChannel], 'o')
                        pylab.xlim([-0.5,X[testChannel].shape[1]/fs-0.5])
                        
                        pylab.show()
                        """
                        
                        assert(abs(maxDeltay[testChannel]) <= delta)
                        
    
                # Print the header. A header is necessary as the column of these
                # files may changes depending on avaiablement event properties.
                if printHeader:        
                    printHeader = False
                    keys = event.properties.keys()
                    
                    header = "no;startTime;duration;"
    
                    for testChannel in channelLst:                    
                        header += "similarity_" + testChannel + ";delay_" + testChannel + ";"
    
                    header += ";".join(keys) + '\n'
                    
                    f.write(header) # Write a string to a file            
    
                result = str(i) + ";" + str(event.startTime) + ";" + str(event.timeLength) + ";" 
                for testChannel in channelLst:
                    result += str(maxCrosscor[testChannel]) + ";"  +  str(maxDeltay[testChannel]) + ";"  
                    
                # Depending on the available properties 
                for eventProperty in keys:
                    if eventProperty in event.properties:
                        result += str(event.properties[eventProperty])
                    result += ";" 
    
                f.write(result[:-1] + '\n') # Write a string to a file
            
        f.close()
            