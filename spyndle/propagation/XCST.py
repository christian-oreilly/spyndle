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
def computeXCST(readerClass, fileName, resPath, eventName, channelLst, 
                 delta =0.5, beforePad=0.5, afterPad=0.5, offset=0.0, eventProperties = [],
                 resFilesPrefix="spindleDelays") :  

    night = os.path.basename(fileName)

    try:
        print 'Loading data of the file ' + night + '...'
        reader = readerClass(fileName)    

        
    except IOError:     
        print "Error : computeXCST : ", fileName, "could not be opened."
        raise      
    
    
    
    for refChannel in channelLst:        
        try:
            f = open(resPath + resFilesPrefix + "_" + night + "_" + refChannel + ".txt", "w")    
        except IOError:     
            print "Error : computeXCST :", resPath + resFilesPrefix + "_" + night + "_" + refChannel + ".txt", "could not be opened."
            raise

        deltaSample = int(delta*reader.getChannelFreq(refChannel))        
        delta       = deltaSample/reader.getChannelFreq(refChannel)    
        
        print "Processing propagation:", night, refChannel   
        
        events = filter(lambda e: e.name == eventName and e.channel == refChannel, reader.events)
        for i, event in zip(range(len(events)), events):
 
            
            startTime = event.startTime -  beforePad      
            duration  = event.timeLength + beforePad + afterPad
        
            signalsDataCmp = reader.read(channelLst,   startTime - delta + offset, duration + 2.0*delta)  
            signalsDataRef = reader.read([refChannel], startTime - delta         , duration + 2.0*delta)  

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

            fs = signalsDataCmp[refChannel].samplingRate                   

            # Spectra computation
            X, fXtmp = computeMST(signalsDataRef[refChannel].signal, fs, m=0.0, k=1.0, fmin=11.0, fmax=16.0)   
            X = abs(transpose(X))
            X = {"Ref":X[:, deltaSample:-deltaSample]}
            for testChannel in channelLst :      
                Xtmp, fXtmp = computeMST(signalsDataCmp[testChannel].signal, fs, m=0.0, k=1.0, fmin=11.0, fmax=16.0)         
                X[testChannel] = abs(transpose(Xtmp))              
        
            refShape = list(shape(X["Ref"]))      
            
            maxCrosscor = {}
            maxDeltay   = {}
            
            crosscor    = zeros(2*deltaSample)
            time        = (arange(2*deltaSample)-deltaSample)/fs
            refSelfCor  = trapz(trapz(X["Ref"]*X["Ref"])) 
            

            for testChannel in channelLst :
                if testChannel == refChannel:
                    maxDeltay[testChannel]   = 0       
                    maxCrosscor[testChannel] = 1.0              
                else:                
                    XCmpSquare = X[testChannel]*X[testChannel]
                    XCmpSquare_colTrapz = trapz(XCmpSquare, axis=0)        
                    for indOffset in range(2*deltaSample):
                        XCmp             = X[testChannel][:, indOffset:(indOffset+refShape[1])]
                        CmpSelfCor       = trapz(XCmpSquare_colTrapz[indOffset:(indOffset+refShape[1])])
                        
                        #print XCmp.shape, X["Ref"].shape
                        crosscor[indOffset] = trapz(trapz(XCmp*X["Ref"]))/max(CmpSelfCor, refSelfCor)

                    ind = where(crosscor == max(crosscor))[0][0]
                    
                    diffCrosscor = diff2(arange(2*deltaSample), crosscor)                
                    
                    if ind == 0 or ind == 2*deltaSample-1:
                        maxDeltay[testChannel]   = time[ind]       
                        maxCrosscor[testChannel] = crosscor[ind]                      
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
                      
    
                        maxDeltay[testChannel]   = time[indFloor]      + indFrac*(time[indFloor+1]-time[indFloor])       
                        maxCrosscor[testChannel] = crosscor[indFloor]  + indFrac*(crosscor[indFloor+1]-crosscor[indFloor])       
                        
                    

            # Print the header. A header is necessary as the column of these
            # files may changes depending on avaiablement event properties.
            if i == 0:        
                keys = event.properties.keys()
                
                header = "no;startTime;duration;"

                for testChannel in channelLst:                    
                    header += "similarity_" + testChannel + ";delay_" + testChannel + ";"

                header += ";".join(keys) + '\n'
                
                f.write(header) # Write a string to a file            


            #print event.toEDFStr()
            result = str(i) + ";" + str(event.startTime) + ";" + str(event.timeLength) + ";" 
            for testChannel in channelLst:
                result += str(maxCrosscor[testChannel]) + ";"  +  str(maxDeltay[testChannel]) + ";"  
                
            # Depending on the available properties 
            for eventProperty in keys:
                if eventProperty in event.properties:
                    result += str(event.properties[eventProperty])
                result += ";" 

            #print night, refChannel, i, len(events)-1

            f.write(result[:-1] + '\n') # Write a string to a file
            
    f.close()
            