# -*- coding: utf-8 -*-


import os,sys, glob, unittest

parentdir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,parentdir + "\\..\\..\\..") 

from spyndle.io.harmonie import HarmonieReader
from spyndle.io.edf import EDFReader
from spyndle.detector import SpindleDectectorRMS
from spyndle.detector import SpindleDectectorAmp
from spyndle.detector import DetectorEvaluator
from spyndle.utils import setUnbufferedPrint


setUnbufferedPrint()


class conversionSaveSpindleTest(unittest.TestCase) :
    
    def testSaving1(self):
        
        print "Converting .sig file to .bdf file..."
        readerSIG =  HarmonieReader(parentdir + "/test.SIG")
        readerSIG.saveAsEDF(parentdir + "/test.BDF", "BDF", verbose=False)
        
        
        print "Reading the .bdf file..."
        reader = EDFReader(parentdir + "/test.BDF")
        
        detectionStages = ["Sleep stage N2", "Sleep stage 2"]
        listChannels    = ['F3-A1', 'F4-A1']
        
        detector = SpindleDectectorRMS(reader, usePickled=False)
        detector.setDetectionStages(detectionStages)
        detector.detectSpindles(listChannels=listChannels)
        
        nbEvents =  len(reader.events) + len(detector.detectedSpindles)
        
        detector.saveSpindle(reader, "SpindleRMS", "Spindle")
        
        reader2 = EDFReader(parentdir + "/test.BDF")
        
        self.failIf(nbEvents != len(reader2.events))

        evaluator = DetectorEvaluator()
        
        evaluator.computeStatistics(detectionStages, "SpindleRMS", 
                          "SpindleRMS",  reader, reader2, listChannels)    
        
        for channel in listChannels:
            try:
                self.assertEqual(evaluator.sensitivity(channel), 1.0)
                self.assertEqual(evaluator.specificity(channel),  1.0)
                self.assertEqual(evaluator.PPV(channel), 1.0)
                self.assertEqual(evaluator.NPV(channel), 1.0)
            except AssertionError:
                print "channel=", channel, "; sensitivity=", evaluator.sensitivity(channel), 
                print "; specificity=", evaluator.specificity(channel), "; PPV=", 
                print evaluator.PPV(channel), "; NPV=", evaluator.NPV(channel)
                raise


    def testSaving2(self):

        channel         = "C4-A1"
        fpath           = parentdir #u'C:/DATA/Vickie_spindles/'
        fileName        = "/test.bdf" #'PETCHn1.bdf'
        AmpEventName    = 'SpindleAmp'        
        detectionStages = ["Sleep stage 2", "Sleep stage N2"]

        #print "Converting .sig file to .bdf file..."
        #readerSIG =  HarmonieReader(fpath + 'PETCHn1.SIG')
        #readerSIG.saveAsEDF(fpath + fileName, "BDF", verbose=False)

        reader1 =  EDFReader(fpath + fileName)
        # Delete any present event with the name AmpEventName
        reader1.events = [e for e in reader1.events if e.name != AmpEventName]
        detector   =  SpindleDectectorAmp(reader1, usePickled=False)  
        detector.setDetectionStages(detectionStages)
        detector.detectSpindles(listChannels=[channel])
        detector.saveSpindle(reader1, AmpEventName, "Spindle")
        
        reader2 = EDFReader(fpath + fileName)
        
        events1 = [e for e in reader1.events if e.name == AmpEventName and e.channel == channel]
        events2 = [e for e in reader2.events if e.name == AmpEventName and e.channel == channel]      
        
        for event1, event2 in zip(events1, events2):
            self.assertEqual(event1, event2)
  




def main():
    unittest.main()

if __name__ == '__main__':
    main()


