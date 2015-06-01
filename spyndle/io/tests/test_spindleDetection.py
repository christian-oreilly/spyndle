# -*- coding: utf-8 -*-

import pylab
import os,sys
parentdir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,parentdir + "\\..\\..\\..") 

from scipy import zeros, where

from spyndle.io import HarmonieReader
from spyndle.io import EDFReader
from spyndle.detector import SpindleDetectorRMS
from spyndle.detector import DetectorEvaluator

print("Reading the .sig file...")
readerSIG =  HarmonieReader(parentdir + "/test.SIG")

print("Conversion: .sig -> .bdf...")
readerSIG.saveAsEDF(parentdir + "/test.BDF", "BDF")

print("Reading the .bdf file...")
readerEDF = EDFReader(parentdir + "/test.BDF")


channelList = readerEDF.getChannelLabels()[0:2]


detectorEDF = SpindleDetectorRMS(readerEDF, usePickled=False)
detectorSIG = SpindleDetectorRMS(readerSIG, usePickled=True)


evaluator = DetectorEvaluator(detectorSIG, detectorEDF)
evaluator.printEvaluation(channelList, ["Sleep stage N2", "Sleep stage 2"])


# For each pages, verify that the signals in reader and readerEDF
# are essentially the same by verifying their signal-to-noise ratio
for noPage in range(readerSIG.getNbPages()) :
    
    # We test only complete pages since incomplete pages will be 
    # different for the two formats, the EDF format padding zeros
    # at the end of the page to make complete records.
    if readerSIG.getInfoPages(noPage+1).isComplete:
        
        sigPage   = readerSIG.readPage(channelList, noPage+1)
        edfPage   = readerEDF.readPage(channelList, noPage+1)  
        timeStart = sigPage.getStartTime() 
        duration  = sigPage.getDuration()  
        
        for channel in channelList:

            print(("noPage=", noPage+1, "(of ", readerSIG.getNbPages(),")" , "channel=", channel))    
            
            sigSpins = [s for s in detectorSIG.detectedSpindles if (s.startTime() >= timeStart and s.startTime() <= timeStart + duration) or
                                        (s.endTime()   >= timeStart and s.endTime()   <= timeStart + duration)]       
            edfSpins = [s for s in detectorEDF.detectedSpindles if (s.startTime() >= timeStart and s.startTime() <= timeStart + duration) or
                                        (s.endTime()   >= timeStart and s.endTime()   <= timeStart + duration)]              
                        
            YSIG = sigPage.getSignal(channel)
            YEDF = edfPage.getSignal(channel)
            TSIG = sigPage.getSignalTime(channel)
            TEDF = edfPage.getSignalTime(channel)
            
            spindleIndSIG = zeros(len(YSIG))
            spindleIndEDF = zeros(len(YEDF))
            
            for spin in sigSpins:
                spindleIndSIG[where( (TSIG >= spin.startTime())*(TSIG <= spin.endTime()))] = max(YSIG)
            for spin in edfSpins:
                spindleIndEDF[where( (TEDF >= spin.startTime())*(TEDF <= spin.endTime()))] = max(YEDF)
            
            if sum(abs(spindleIndSIG  - spindleIndEDF)) >  max(YSIG)*10 :
                for spin in sigSpins:
                    print(("sig spin : ", spin.startTime(), spin.endTime()))
                for spin in edfSpins:
                    print(("edf spin : ", spin.startTime(), spin.endTime()))

                
                pylab.plot(TSIG, YSIG)
                pylab.plot(TEDF, YEDF)
                pylab.plot(TSIG, spindleIndSIG)
                pylab.plot(TEDF, spindleIndEDF)
                pylab.show()
    else:
        print(("Skiping incomplete page " + str(noPage+1)))




