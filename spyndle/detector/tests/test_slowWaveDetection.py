# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 22:00:51 2014

@author: oreichri
"""

import os, sys

from spyndle import DevuystKC

parentdir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(parentdir, "..")) 
sys.path.append(os.path.join(parentdir, "..", "..", "io")) 

from slowWaveDetectors import MassiminiSlowWaveDetector
from edf import EDFReader 
from detectorEvaluation import DetectorEvaluator



def main():
    
    # Loading example signal from Devuyst database
    filename = os.path.join(DevuystKC._devuystDBPath_, "excerpt1.BDF")
    #if not os.path.isfile(filename):
    from spyndle.DevuystKC import convertDevuystToBDF
    convertDevuystToBDF()
      
    reader = EDFReader(filename)
    
    detector = MassiminiSlowWaveDetector(reader)
    
    detector.detectEvents()
    print((reader.events))
    print((len(detector.detectedEvents)))

    detector.saveEvents(reader, eventName="MassiminiSL", eventGroupName="SlowWave")  
    
    listDetectionStages = ["Sleep stage 2", "Sleep stage 3", "Sleep stage 4"]

    evaluator = DetectorEvaluator()
    evaluator.computeStatistics(listDetectionStages, "KComplexV1", "MassiminiSL",  reader)
    evaluator.printEvaluation(listDetectionStages, "KComplexV1", "MassiminiSL",  reader)

    evaluator.computeStatistics(listDetectionStages, "KComplexV2", "MassiminiSL",  reader)
    evaluator.printEvaluation(listDetectionStages, "KComplexV2", "MassiminiSL",  reader)
    
    evaluator.computeStatistics(listDetectionStages, "KComplexAu", "MassiminiSL",  reader)
    evaluator.printEvaluation(listDetectionStages, "KComplexAu", "MassiminiSL",  reader)
    
    evaluator.computeStatistics(listDetectionStages, "KComplexV1", "KComplexAu",  reader)
    evaluator.printEvaluation(listDetectionStages, "KComplexV1", "KComplexAu",  reader)


    evaluator.computeStatistics(listDetectionStages, "KComplexV2", "KComplexAu",  reader)
    evaluator.printEvaluation(listDetectionStages, "KComplexV2", "KComplexAu",  reader)
    

    


if __name__ == '__main__':
    main()
