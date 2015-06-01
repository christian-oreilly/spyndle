# -*- coding: utf-8 -*-

import os
parentdir = os.path.dirname(os.path.abspath(__file__))
from urllib.request import urlretrieve
from urllib.error import ContentTooShortError

url      = "https://bitbucket.org/christian_oreilly/spyndle/downloads/test"
fileName = os.path.join(parentdir, "test")

from spyndle.io.harmonie import HarmonieReader
from spyndle.io.edf import EDFReader
from scipy import log10, trapz
import numpy as np

import unittest


# This does not work anymore on Python 3 because of certificate verification
# Please download and put the files at the appropriate place manually
# before lauching this code.
"""
def loadFiles():
    print("Reading the .sig file...")
    if not os.path.exists(fileName + ".SIG"):
        try:
            print("Loading the .sig file from Internet.")
            urlretrieve(url + ".SIG", fileName + ".SIG")
        except ContentTooShortError:
             print("The retreived file is shorter than expected. The download as "\
                   "probably been interrupted")
             exit

    if not os.path.exists(fileName + ".sts"):
        try:
            print("Loading the .sts file from Internet.")
            urlretrieve(url + ".sts", fileName + ".sts")
        except ContentTooShortError:
             print("The retreived file is shorter than expected. The download as "\
                   "probably been interrupted")
             exit
"""                     

class conversionHarmonieEDFTests(unittest.TestCase) :
    
    
    def testConversion(self):
        
        #loadFiles()
                 
        readerSIG =  HarmonieReader(os.path.join(parentdir, "test.SIG"))

        print("Saving a copy in the BDF format...")
        readerSIG.saveAsEDF(os.path.join(parentdir, "test.bdf"), "BDF",
                            annotationFileName=os.path.join(parentdir, "test Annotations.bdf"))
        

        print("Reading the saved BDF file...")
        readerEDF = EDFReader(os.path.join(parentdir, "test.BDF"),
                              annotationFileName=os.path.join(parentdir, "test Annotations.bdf"))

        channelListSig = readerSIG.getChannelLabels()
        channelListEdf = readerEDF.getChannelLabels()
        # For each pages, verify that the signals in reader and readerEDF
        # are essentially the same by verifying their signal-to-noise ratio
        for noPage in range(readerSIG.getNbRecords()) :
            
            # We test only complete pages since incomplete pages will be 
            # different for the two formats, the EDF format padding zeros
            # at the end of the page to make complete records.
            if readerSIG.getInfoRecords(noPage+1).isComplete:
                pageSIG = readerSIG.readRecord(channelListSig, noPage+1)
                pageEDF = readerEDF.readRecord(channelListEdf, noPage+1)               
                
                for channelSig, channelEdf in zip(channelListSig, channelListEdf):
                    Y1 = pageSIG.getSignal(channelSig)
                    Y2 = pageEDF.getSignal(channelEdf)
                    
                    SNR = 10.0*log10(trapz(Y1**2)/trapz((Y1-Y2)**2))
                    print(("noPage=", noPage+1, "(of ", readerSIG.getNbRecords(),")" , "channel=", channelSig, "SNR=", SNR))
                    
                    if channelSig != "Mic-Mic" and SNR < 50.0:
                        import pylab
                        pylab.plot(list(range(len(Y1))), Y1)
                        pylab.plot(list(range(len(Y2))), Y2)
                        pylab.show()
                    self.failIf(channelSig != "Mic-Mic" and SNR < 50.0)
            else:
                print(("Skiping incomplete page " + str(noPage+1)))

        
    def testConversion_readChannel(self):
        

        #loadFiles()        

        print("Reading the .sig file...")
        readerSIG =  HarmonieReader(os.path.join(parentdir, "test.SIG"))
        
        print("Saving a copy in the BDF format...")
        readerSIG.saveAsEDF(os.path.join(parentdir, "test.bdf"), "BDF",
                              annotationFileName=os.path.join(parentdir, "test Annotations.bdf"))
        
        print("Reading the saved BDF file...")
        readerEDF = EDFReader(os.path.join(parentdir, "test.bdf"),
                              annotationFileName=os.path.join(parentdir, "test Annotations.bdf"))

        channelListSig = readerSIG.getChannelLabels()
        channelListEdf = readerEDF.getChannelLabels()

        for channelSig, channelEdf in zip(channelListSig, channelListEdf):    
            
            data1 = readerSIG.readChannel(channelSig, False)
            Y1    = data1.signal
            T1    = readerSIG.getChannelTime(channelSig)
            
            #fsSIG       = data.samplingRate

            data2 = readerEDF.readChannel(channelEdf)
            Y2    = data2.signal
            T2    = readerEDF.getChannelTime(channelEdf)
            #fsEDF       = data.samplingRate
           
            print(("1:", np.mean(abs(Y1)) , np.mean(abs(Y2))))            
            
           
            #print sum(np.in1d(T2, T1)), sum(np.in1d(T1, T2)), len(T1),  len(T2)    
           
            #SNR = 10.0*log10(trapz(Y1**2)/trapz((Y1-Y2)**2))
            #print "channel=", channel, "SNR=", SNR
            
            #if SNR < 50.0:
            #import pylab
            #pylab.plot(T1, Y1)
            #pylab.plot(T2, Y2)
            #pylab.show()
        
            #self.failIf(SNR < 50.0)    
    
    


def main():
    unittest.main()

if __name__ == '__main__':
    main()
