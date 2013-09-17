# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 12:38:40 2012

@author: REVESTECH
"""

from spyndle.io import DevuystReader, EDFReader

from scipy import unique
import os

def convertDevuystToBDF():
    path = os.path.dirname(__file__)
    fpath        = u'excerpt'
    pathV1       = u'Visual_scoring1_excerpt'
    pathV2       = u'Visual_scoring2_excerpt'
    pathAu       = u'Automatic_detection_excerpt'
    pathHy       = u'Hypnogram_excerpt'
    
    samplingRates = [100.0, 200.0, 50.0, 200.0, 200.0, 200.0, 200.0, 200.0]

    
    for i in range(1, 9):
        try:
            print 'Loading data for subject #' + str(i) + '...'
            reader = DevuystReader(os.path.join(path, fpath + str(i) + ".txt"), samplingRates[i-1])    
            reader.importEvents(os.path.join(path, pathV1 + str(i) + ".txt"), "spindleV1")
            reader.importEvents(os.path.join(path, pathV2 + str(i) + ".txt"), "spindleV2")
            reader.importEvents(os.path.join(path, pathAu + str(i) + ".txt"), "spindleAu")
            reader.importHypnogram(os.path.join(path, pathHy + str(i) + ".txt"))
    
            print 'Saving ' + fpath + str(i) + ".BDF" + '...'        
            reader.saveAsEDF(os.path.join(path, fpath + str(i) + ".BDF"), fileType = "BDF")   
                
        except IOError:     
            print "Error: The selected file could not be open."
            exit()        
    
    
    
