# -*- coding: utf-8 -*-

import os,sys
parentdir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,parentdir + "\\..\\..\\..") 


from spyndle.io.harmonie import HarmonieReader
from spyndle.io.edf import EDFReader
from scipy import log10, trapz
import numpy as np

import unittest

class conversionHarmonieEDFTests(unittest.TestCase) :
    

    def testConversion(self):
        print "Reading the .sig file..."
        readerSIG =  HarmonieReader(parentdir + "/test.SIG")
        
        print "Saving a copy in the BDF format..."
        readerSIG.saveAsEDF(parentdir + "/test.BDF", "BDF")
        
        print "Reading the saved BDF file..."
        readerEDF = EDFReader(parentdir + "/test.BDF")

        channelList = readerSIG.getChannelLabels()
        # For each pages, verify that the signals in reader and readerEDF
        # are essentially the same by verifying their signal-to-noise ratio
        for noPage in range(readerSIG.getNbPages()) :
            
            # We test only complete pages since incomplete pages will be 
            # different for the two formats, the EDF format padding zeros
            # at the end of the page to make complete records.
            if readerSIG.getInfoPages(noPage+1).isComplete:
                pageSIG = readerSIG.readPage(channelList, noPage+1)
                pageEDF = readerEDF.readPage(channelList, noPage+1)               
                
                for channel in channelList:
                    Y1 = pageSIG.getSignal(channel)
                    Y2 = pageEDF.getSignal(channel)
                    
                    SNR = 10.0*log10(trapz(Y1**2)/trapz((Y1-Y2)**2))
                    print "noPage=", noPage+1, "(of ", readerSIG.getNbPages(),")" , "channel=", channel, "SNR=", SNR
                    
                    if channel != "Mic-Mic" and SNR < 50.0:
                        import pylab
                        pylab.plot(range(len(Y1)), Y1)
                        pylab.plot(range(len(Y2)), Y2)
                        pylab.show()
                    self.failIf(channel != "Mic-Mic" and SNR < 50.0)
            else:
                print "Skiping incomplete page " + str(noPage+1)

    
    def testConversion_readChannel(self):
        print "Reading the .sig file..."
        readerSIG =  HarmonieReader(parentdir + "/test.SIG")
        
        print "Saving a copy in the BDF format..."
        readerSIG.saveAsEDF(parentdir + "/test.BDF", "BDF")
        
        print "Reading the saved BDF file..."
        readerEDF = EDFReader(parentdir + "/test.BDF")

        channelList = readerSIG.getChannelLabels()

        for channel in channelList:    
            
            data1 = readerSIG.readChannel(channel, False)
            Y1    = data1.signal
            T1    = readerSIG.getChannelTime(channel)
            
            #fsSIG       = data.samplingRate

            data2 = readerEDF.readChannel(channel)
            Y2    = data2.signal
            T2    = readerEDF.getChannelTime(channel)
            #fsEDF       = data.samplingRate
           
            print "1:", np.mean(abs(Y1)) , np.mean(abs(Y2))            
            
           
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
