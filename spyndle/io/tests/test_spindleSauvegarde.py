# -*- coding: utf-8 -*-


import os,sys, unittest, io

from shutil import copyfile

parentdir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,parentdir + "\\..\\..\\..") 

from spyndle.io.harmonie import HarmonieReader
from spyndle.io.edf import EDFReader
from spyndle.detector import SpindleDectectorAmp
from spyndle.detector import DetectorEvaluator
from spyndle.utils import setUnbufferedPrint


setUnbufferedPrint()


class conversionSaveSpindleTest(unittest.TestCase) :


    def testSaving1(self):
        
        print("Converting .sig file to .bdf file...")
        readerSIG =  HarmonieReader(parentdir + "/test.SIG")
        readerSIG.saveAsEDF(parentdir + "/test.BDF", "BDF", verbose=False)
        
        
        print("Reading the .bdf file...")
        reader = EDFReader(parentdir + "/test.BDF")
        
        detectionStages = ["Sleep stage N2", "Sleep stage 2"]
        channelList    = ['F3-A1', 'F4-A1']
        
        detector = SpindleDectectorAmp(reader, usePickled=False)
        detector.setDetectionStages(detectionStages)
        detector.detectSpindles(channelList=channelList)
        
        nbEvents =  len(reader.events) + len(detector.detectedSpindles)
        
        detector.saveSpindle(reader, "SpindleRMS", "Spindle")
        
        reader2 = EDFReader(parentdir + "/test.BDF")
        
        self.failIf(nbEvents != len(reader2.events))

        evaluator = DetectorEvaluator()
        
        evaluator.computeStatistics(detectionStages, "SpindleRMS", 
                          "SpindleRMS",  reader, reader2, channelList)    
        
        for channel in channelList:
            try:
                self.assertEqual(evaluator.sensitivity(channel), 1.0)
                self.assertEqual(evaluator.specificity(channel),  1.0)
                self.assertEqual(evaluator.PPV(channel), 1.0)
                self.assertEqual(evaluator.NPV(channel), 1.0)
            except AssertionError:
                print("channel=", channel, "; sensitivity=", evaluator.sensitivity(channel), end=' ') 
                print("; specificity=", evaluator.specificity(channel), "; PPV=", end=' ') 
                print(evaluator.PPV(channel), "; NPV=", evaluator.NPV(channel))
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
        detector.detectSpindles(channelList=[channel])
        detector.saveSpindle(reader1, AmpEventName, "Spindle")
        
        reader2 = EDFReader(fpath + fileName)
        
        events1 = [e for e in reader1.events if e.name == AmpEventName and e.channel == channel]
        events2 = [e for e in reader2.events if e.name == AmpEventName and e.channel == channel]      
        
        for event1, event2 in zip(events1, events2):
            self.assertEqual(event1, event2)




    def testDoubleDetection(self):
    
        fpath           = parentdir #u'C:/DATA/Vickie_spindles/'
        fileName        = "/test.bdf" #'PETCHn1.bdf'        
    
        channels = ["C4-A1"]

        print("Converting .sig file to .bdf file...")
        readerSIG =  HarmonieReader(fpath + fileName[:-3] + 'SIG')
        readerSIG.saveAsEDF(fpath + fileName, "BDF", verbose=False)
    
        reader =  EDFReader(fpath + fileName)

        data        = reader.readChannel(channels[0], usePickled=False)
        signal      = data.signal
        print(signal[[1, 100, 200, 300, 400, 500]])       
       
        data        = reader.readChannel(channels[0], usePickled=False)
        signal      = data.signal
        print(signal[[1, 100, 200, 300, 400, 500]])       



    def testDoubleDetection(self):
    
        AmpEventName    = 'SpindleAmp'
        RMSEventName    = 'SpindleRMS'
        
        fpath           = parentdir #u'C:/DATA/Vickie_spindles/'
        fileName        = "/test.bdf" #'PETCHn1.bdf'        
        
        detectionStages = ["Sleep stage 2", "Sleep stage N2"]
                
        channels = ["C4-A1"]

        print("Converting .sig file to .bdf file...")
        readerSIG =  HarmonieReader(fpath + fileName[:-3] + 'SIG')
        readerSIG.saveAsEDF(fpath + fileName, "BDF", verbose=False)
    
        reader =  EDFReader(fpath + fileName)

        # Delete any present event with the name AmpEventName
        #reader.events = [e for e in reader.events if e.name != AmpEventName]
        #reader.events = [e for e in reader.events if e.name != RMSEventName]
        
        detector      =  SpindleDectectorAmp(reader, usePickled=False)  
        detector.quantileThreshold = 0.925
        detector.setDetectionStages(detectionStages)
        detector.computeRMS       = False
        detector.computeFreq      = False       
        detector.computeSlopeFreq = False          
        detector.detectSpindles(channelList=channels)
        nbSpins1 = len(detector.detectedSpindles)
        detector.saveSpindle(reader, AmpEventName, "Spindle")      
         
        detector.detectSpindles(channelList=channels)
        nbSpins2 = len(detector.detectedSpindles)       
        detector.saveSpindle(reader, RMSEventName, "Spindle")  
       
        try:
            self.assertEqual(nbSpins1, nbSpins2)
        except AssertionError:
            print(nbSpins1, nbSpins2)
            raise

        evaluatorAmp = DetectorEvaluator()   
        evaluatorAmp.computeStatistics(detectionStages, RMSEventName, 
                                  AmpEventName, reader, channelList=channels) 

        for channel in channels:   
            if channel in evaluatorAmp.TP:
                self.assertEqual(evaluatorAmp.sensitivity(channel), 1.0)
                self.assertEqual(evaluatorAmp.specificity(channel), 1.0)
                self.assertEqual(evaluatorAmp.PPV(channel), 1.0)
                self.assertEqual(evaluatorAmp.NPV(channel), 1.0)
    


    def testLoadSaveLoad(self):
    
        fpath           = parentdir #u'C:/DATA/Vickie_spindles/'
        fileName        = "/test.bdf" #'PETCHn1.bdf'        

        
        print("Converting .sig file to .bdf file...")
        readerSIG =  HarmonieReader(fpath + fileName[:-3] + 'SIG')
        readerSIG.saveAsEDF(fpath + fileName, "BDF", verbose=False)
    
        reader =  EDFReader(fpath + fileName)
        reader.saveAs(fpath + "temp.bdf")
        del reader

        with io.open(fpath + fileName, 'rb') as file1:
            with io.open(fpath + "temp.bdf", 'rb') as file2:
                
                while 1:
                    bytes1 = file1.read(100)
                    bytes2 = file2.read(100)

                    if not bytes1:
                        break                
                
                    self.assertEqual(bytes1, bytes2)
        
        os.remove(fpath + "temp.bdf")

 
    def testDoubleSave(self):
    
        fpath           = parentdir #u'C:/DATA/Vickie_spindles/'
        fileName        = "/test.bdf" #'PETCHn1.bdf'        

        
        print("Converting .sig file to .bdf file...")
        readerSIG =  HarmonieReader(fpath + fileName[:-3] + 'SIG')
        readerSIG.saveAsEDF(fpath + fileName, "BDF", verbose=False)
    
        reader =  EDFReader(fpath + fileName)
        reader.saveAs(fpath + "temp.bdf")
        reader.saveAs(fpath + "temp2.bdf")
        del reader

        with io.open(fpath + "temp.bdf", 'rb') as file2:
            with io.open(fpath + "temp2.bdf", 'rb') as file3:
                    
                while 1:
                    bytes3 = file3.read(100)
                    bytes2 = file2.read(100)
                    if not bytes3:
                        break
                    self.assertEqual(bytes3, bytes2)
                
        os.remove(fpath + "temp.bdf")
        os.remove(fpath + "temp2.bdf")


    """
     Convert from SIG to BDF, detect, save, make a copy, reopen the original
     detect, save, compare the original with the copy
    """
    def testLoadDetectSaveLoad(self):
    
        fpath           = parentdir #u'C:/DATA/Vickie_spindles/'
        fileName        = "/test.bdf" #'PETCHn1.bdf'        

        AmpEventName    = 'SpindleAmp'
        RMSEventName    = 'SpindleRMS'

        detectionStages = ["Sleep stage 2", "Sleep stage N2"]
                
        channels = ["C4-A1"]

        print("Converting .sig file to .bdf file...")
        readerSIG =  HarmonieReader(fpath + fileName[:-3] + 'SIG')
        readerSIG.saveAsEDF(fpath + fileName, "BDF", verbose=False)
        
        reader =  EDFReader(fpath + fileName)

        # Delete any present event with the name AmpEventName
        reader.events = [e for e in reader.events if e.name != AmpEventName]
        reader.events = [e for e in reader.events if e.name != RMSEventName]
        
        detector      =  SpindleDectectorAmp(reader, usePickled=False)  
        detector.quantileThreshold = 0.925
        detector.setDetectionStages(detectionStages)
        detector.computeRMS       = False
        detector.computeFreq      = False       
        detector.computeSlopeFreq = False          
        detector.detectSpindles(channelList=channels)
        detector.saveSpindle(reader, AmpEventName, "Spindle")      

        copyfile(fpath + fileName, fpath + "temp.bdf")         


        reader =  EDFReader(fpath + fileName)

        # Delete any present event with the name AmpEventName
        reader.events = [e for e in reader.events if e.name != AmpEventName]
        reader.events = [e for e in reader.events if e.name != RMSEventName]
        
        detector      =  SpindleDectectorAmp(reader, usePickled=False)  
        detector.quantileThreshold = 0.925
        detector.setDetectionStages(detectionStages)
        detector.computeRMS       = False
        detector.computeFreq      = False       
        detector.computeSlopeFreq = False          
        detector.detectSpindles(channelList=channels)
        detector.saveSpindle(reader, AmpEventName, "Spindle")    
        

        with io.open(fpath + fileName, 'rb') as file1:
            with io.open(fpath + "temp.bdf", 'rb') as file2:
                
                byte1 = file1.read(1)
                byte2 = file2.read(1)
                self.assertEqual(byte1, byte2)
        
        os.remove(fpath + "temp.bdf")



def main():
    unittest.main()

if __name__ == '__main__':
    main()












