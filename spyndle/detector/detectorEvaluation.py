# -*- coding: utf-8 -*-

"""
    Code allowing evaluation of spindle detectors. 

    Copyright (C) 2013  Christian O'Reilly

    For personnal, educationnal, and research purpose, this software is 
    provided under the GNU GPL (V.3) license: you can redistribute it and/or
    modify it under the terms of the version 3 of the GNU General Public 
    License as published by the Free Software Foundation.
          
    To use this software in commercial application, please contact the author.
    
    If this code is used for research purpose, the reference [1] should be
    cited in the derived publication to refer the reader to the description 
    of the methodology.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


 Author: Christian O'Reilly (christian.oreilly@umontreal.ca)
 Date  : April 24, 2013

 [1] O'Reilly, C. & Nielsen, T. Sleep spindle detection: Automation and 
     performance evaluation using fine temporal resolution, Submitted to
     IEEE Transactions on Biomedical Engineering, april 2013. 

"""

from scipy import sqrt
import numpy as np


# Compare two detections. Can be either used to compare two detector objects
# or two compare two detections already performed and retrieved with readers.
class DetectorEvaluator:
    def __init__(self):
        # goldStandard ans tested are two detectors correctly initiated
        # with their internal reader.

        self.N  = {}
        self.TP = {}
        self.TN = {}
        self.FP = {}
        self.FN = {}

        self.spindlesLoaded = False
    
    
    def detectSpindlesAndDetectStatistics(self, goldStandardDetector, testedDetector, 
                                          channelList, listSleepStages, 
                                          goldSpindleName = "goldSpindle", 
                                          testedSpindleName = "testedSpindle"):
        # goldStandardDetector ans testedDetector are two detectors correctly
        # initiated with their internal reader.
    

        # We want the stage identification to be case-insensitive
        #for i in range(len(listSleepStages)): 
        #    listSleepStages[i] = listSleepStages[i].lower()
    
        ########### PERFORM SPINDLE DETECTION
        print "Detecting spindles using the gold standard..."
        self.goldStandard.setDetectionStages(listSleepStages)
        self.goldStandard.detectSpindles(channelList) 
        self.goldStandard.saveSpindle(self.goldStandard.reader, goldSpindleName, "Spindle")


        print "Detecting spindles using the tested detector..."
        self.tested.setDetectionStages(listSleepStages)
        self.tested.detectSpindles(channelList) 
        self.tested.saveSpindle(self.tested.reader, testedSpindleName, "Spindle")

        self.computeStatistics(self, listSleepStages, goldSpindleName, testedSpindleName,
                               self.goldStandard.reader, self.tested.reader, channelList)
        
        
        
        

    def computeStatistics(self, listDetectionStages, nameEventGold, 
                          nameEventTested,  readerGold, readerTested=None, channelList=None):

        # If no tested reader is specified, we consider the that reader passed as
        # readerGold contains the events for both the gold and the tested detection
        if readerTested is None:
            readerTested = readerGold

        if isinstance(channelList, (str, unicode)):
            channelList = [channelList]
         
        # If channelList is not specified, we consider every channel present
        # in the gold reader.
        if channelList is None:
            channelList = readerGold.getChannelLabels()


        for channel in channelList:
            indStages    = readerGold.getEventIndicator(listDetectionStages, channel, globalEvent=True)            
            indEventGold = readerGold.getEventIndicator(nameEventGold, channel)
            indEventTest = readerGold.getEventIndicator(nameEventTested, channel)
            
            TP = np.logical_and(indEventGold, indEventTest)
            TN = np.logical_and(np.logical_not(indEventGold), np.logical_not(indEventTest)) 
            FP = np.logical_and(np.logical_not(indEventGold), indEventTest)
            FN = np.logical_and(indEventGold, np.logical_not(indEventTest))
            
            self.TP[channel] = np.sum(np.logical_and(indStages, TP))
            self.TN[channel] = np.sum(np.logical_and(indStages, TN))
            self.FP[channel] = np.sum(np.logical_and(indStages, FP))
            self.FN[channel] = np.sum(np.logical_and(indStages, FN))         
       
        
                  
    def printEvaluation(self, listSleepStages, nameEventGold, 
                          nameEventTested, readerGold, readerTested=None, channelList=None):

        self.computeStatistics(listSleepStages, nameEventGold, 
                               nameEventTested, readerGold, readerTested, channelList)
        for channel in self.FP:            
            print ("Channel:%s, sensitivity=%f, specificity=%f, PPV=%f, NPV=%f\n"
                   "            TP=%f, TN=%f, FP=%f, FN=%f" % (channel, 
                   self.sensitivity(channel), self.specificity(channel),
                   self.PPV(channel), self.NPV(channel), self.TP[channel], 
                   self.TN[channel], self.FP[channel], self.FN[channel]))

          
    def sensitivity(self, channel):  
        if self.TP[channel] > 0 :
            return self.TP[channel]/float(self.TP[channel] + self.FN[channel])    
        else:
            return 0.0        


    def specificity(self, channel):  
        if self.TN[channel] > 0 :
            return self.TN[channel]/float(self.FP[channel] + self.TN[channel])
        else:
            return 0.0 


    def accuracy(self, channel):  
        if self.TN[channel] > 0 :
            N = self.FN[channel] + self.TN[channel] +  self.FP[channel] + self.TP[channel]
            return (self.TN[channel] + self.TP[channel])/N
        else:
            return 0.0                
        
        
    def PPV(self, channel):    # Positive predictive value                 
        if self.TP[channel] > 0 :
            return self.TP[channel]/float(self.TP[channel] + self.FP[channel])      
        else:
            return 0.0
                  
                  
    def NPV(self, channel):    # Negative predictive value                        
        if self.TN[channel] > 0 :
            return self.TN[channel]/float(self.TN[channel] + self.FN[channel])
        else:
            return 0.0 
        
                          
    def MCC(self, channel):    # Negative predictive value                        
        if self.TN[channel] > 0 :
            num = float(self.TP[channel]*self.TN[channel] - self.FP[channel]*self.FN[channel])
            P   = self.TP[channel] + self.FN[channel]
            P2  = self.TP[channel] + self.FP[channel]
            N   = self.FP[channel] + self.TN[channel]
            N2  = self.FN[channel] + self.TN[channel]
            den = sqrt(P*P2*N*N2)
            return num/den
        else:
            return 0.0 
        
        
                
                
                          
    def randomAggreProb(self, channel):    # Negative predictive value                        
        if self.TN[channel] > 0 :
            P   = self.TP[channel] + self.FN[channel]
            P2  = self.TP[channel] + self.FP[channel]
            N   = self.FP[channel] + self.TN[channel]
            N2  = self.FN[channel] + self.TN[channel]
            return float((P2*P + N2*N))/float((P+N)**2)
        else:
            return 0.0 
        
        
    def cohenk(self, channel):    # Negative predictive value                        
        if self.TN[channel] > 0 :
            Pe  = self.randomAggreProb(channel)
            acc = self.accuracy(channel) 
            return (acc - Pe)/(1-Pe)
        else:
            return 0.0 
                        
        

            
            
            
            