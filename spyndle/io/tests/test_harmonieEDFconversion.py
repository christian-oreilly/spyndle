# -*- coding: utf-8 -*-

import os,sys
parentdir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,parentdir + "\\..\\..\\..") 


from spyndle.io.harmonie import HarmonieReader
from spyndle.io.edf import EDFReader
from scipy import log10, trapz

import unittest

class conversionHarmonieEDFTests(unittest.TestCase) :
    
    def testConversion(self):
        print "Reading the .sig file..."
        reader =  HarmonieReader(parentdir + "/test.SIG")
        
        print "Saving a copy in the BDF format..."
        reader.saveAsEDFA(parentdir + "/test.BDF", "BDF")
        
        print "Reading the saved BDF file..."
        readerEDF = EDFReader(parentdir + "/test.BDF")

        # For each pages, verify that the signals in reader and readerEDF
        # are essentially the same by verifying their signal-to-noise ratio
        for noPage in range(reader.getNbPages()) :
            for channel in reader.getChannelLabels():
                Y1  = reader.readPage([channel], noPage+1)[2][0]
                Y2  = readerEDF.readPage([channel], noPage+1)[2][0]
                SNR = 10.0*log10(trapz(Y1**2)/trapz((Y1-Y2)**2))
                print "noPage=", noPage+1, "(of ", reader.getNbPages(),")" , "channel=", channel, "SNR=", SNR
                self.failIf(SNR < 50.0)

def main():
    unittest.main()

if __name__ == '__main__':
    main()
