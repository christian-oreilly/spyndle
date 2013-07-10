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

from scipy import concatenate
from operator import attrgetter


class BadTransitionSequence(Exception): 
    
    def __init__(self, transition, state):
        if transition.type == "A":
            self.transition = "beginStage"
        elif transition.type == "B":
            self.transition = "endStage"
        elif transition.type == "C":
            self.transition = "beginSpindleGold"
        elif transition.type == "D":
            self.transition = "endSpindleGold"
        elif transition.type == "E":
            self.transition = "beginSpindleTest"
        elif transition.type == "F":
            self.transition = "endSpindleTest"
        
        self.state = state
        self.time  = transition.time        
        
    def __str__(self):
        return ("Encountered invalid transition sequence. Trying to make the transition "
                            + self.transition + " while in state " + self.state + 
                            " at time " + str(self.time) + ".")




class transition:
    def __init__(self, time, type):
        self.time = time
        
        """
        Type legend:
            A: beginStage
            B: endStage
            C: beginSpindleGold
            D: endSpindleGold
            E: beginSpindleTest
            F: endSpindleTest    
        """
        self.type = type

    def __str__(self):
        return str(self.time) + ":" + self.type



class State:
    pass

class StateITN(State):
    def applyTransition(self, transition, machine):
        if transition.type == "A":
            machine.lastTime = transition.time
            return StateTN()
        elif transition.type == "E":
            return StateIFP()
        elif transition.type == "C":
            return StateIFN()
        else:
            raise BadTransitionSequence(transition, "ITN")
                
class StateIFN(State):
    def applyTransition(self, transition, machine):
        if transition.type == "A":
            machine.lastTime = transition.time
            return StateFN()
        elif transition.type == "E":
            return StateITP()
        elif transition.type == "D":
            return StateITN()
        else:
            raise BadTransitionSequence(transition, "IFN")
                
                
class StateITP(State):
    def applyTransition(self, transition, machine):
        
        if transition.type == "A":
            machine.lastTime = transition.time
            return StateTP()
        elif transition.type == "F":
            return StateIFN()
        elif transition.type == "D":
            return StateIFP()
        else:
            raise BadTransitionSequence(transition, "ITP")
                
                
class StateIFP(State):
    def applyTransition(self, transition, machine):

        if transition.type == "A":
            machine.lastTime = transition.time
            return StateFP()
        elif transition.type == "C":
            return StateITP()
        elif transition.type == "F":
            return StateITN()
        else:
            raise BadTransitionSequence(transition, "IFP")
                

                
class StateTN(State):
    def applyTransition(self, transition, machine):
        
        if transition.type   == "B":  state = StateITN()
        elif transition.type == "E":  state = StateFP()
        elif transition.type == "C":  state = StateFN()
        else:
            raise BadTransitionSequence(transition, "TN")
                
        machine.TN += transition.time - machine.lastTime
        machine.lastTime = transition.time
        return state    

            
class StateFN(State):
    def applyTransition(self, transition, machine):
        
        if transition.type   == "B":  state = StateIFN()
        elif transition.type == "D":  state = StateTN()
        elif transition.type == "E":  state = StateTP()
        else:
            raise BadTransitionSequence(transition, "FN")
                
        machine.FN += transition.time - machine.lastTime
        machine.lastTime = transition.time
        return state

            
class StateTP(State):
    def applyTransition(self, transition, machine):
        
        if transition.type   == "B":  state = StateITP()
        elif transition.type == "F":  state = StateFN()
        elif transition.type == "D":  state = StateFP()
        else:
            raise BadTransitionSequence(transition, "TP")
                
        machine.TP += transition.time - machine.lastTime
        machine.lastTime = transition.time
        return state

            
class StateFP(State):
    def applyTransition(self, transition, machine):
        
        if transition.type   == "B":  state = StateIFP()
        elif transition.type == "F":  state = StateTN()
        elif transition.type == "C":  state = StateTP()
        else:
            raise BadTransitionSequence(transition, "FP")
                
        machine.FP += transition.time - machine.lastTime
        machine.lastTime = transition.time
        return state




class StateMachine:
    def __init__(self):
        self.TN         = 0.0
        self.FN         = 0.0
        self.TP         = 0.0
        self.FP         = 0.0
        self.lastTime   = None
    
    def run(self, transitionList):
        
        state = StateITN()
        for transition in transitionList:
            state = state.applyTransition(transition, self)
            


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
                                          listChannels, listSleepStages, 
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
        self.goldStandard.detectSpindles(listChannels) 
        self.goldStandard.saveSpindle(self.goldStandard.reader, goldSpindleName, "Spindle")


        print "Detecting spindles using the tested detector..."
        self.tested.setDetectionStages(listSleepStages)
        self.tested.detectSpindles(listChannels) 
        self.tested.saveSpindle(self.tested.reader, testedSpindleName, "Spindle")

        self.computeStatistics(self, listChannels, listSleepStages, self.goldStandard.reader, 
                               goldSpindleName, self.tested.reader, testedSpindleName)
        
        
        
        

    def computeStatistics(self, listSleepStages, nameEventGold, 
                          nameEventTested,  readerGold, readerTested=None, listChannels=None):

        # If no tested reader is specified, we consider the that reader passed as
        # readerGold contains the events for for both the gold and the tested detection
        if readerTested is None:
            readerTested = readerGold
         
        # If listChannels is not specified, we consider every channel present
        # in the gold reader.
        if listChannels is None:
            listChannels = readerGold.getChannelLabels()
         
        def getSpindleTransitions(reader, eventName, channel, isGold):
            channelSpindles = filter(lambda s: s.channel == channel and s.name == eventName, reader.events)     
            if len(channelSpindles) :
                if isGold :
                    #for s in channelSpindles:
                    #    print "Gold:", s                    
                    return concatenate([[transition(spindle.timeStart(), "C"), 
                                              transition(spindle.timeEnd(), "D")]  for spindle in channelSpindles])   
                else:
                    #for s in channelSpindles:
                    #    print "Test:", s
                    return concatenate([[transition(spindle.timeStart(), "E"), 
                                              transition(spindle.timeEnd(), "F")]  for spindle in channelSpindles])          
            return []
            
       
        # COMPUTING THE TRANSITIONS ASSOCIATED TO SLEEP STAGES    
        # computing arrays of sample index associated to every sleep stages for
        # which we want to evaluate the spindle detection. We only use the changes
        # of stage for this computation since we do not want to have, for example,
        # a serie of endStage2, beginStage2, endStage2, beginStage2 creating artigicial
        # discontinuities at each change of page
        stageEventNames =  [e for e in readerGold.events if e.groupeName == "Stage"]  
        
        if stageEventNames[0].name in listSleepStages:
            self.stageTransitions = [transition(stageEventNames[0].timeStart(), "A")]
        else:             
            self.stageTransitions = []
            
        for e1, e2 in zip(stageEventNames[:-1], stageEventNames[1:]):
            if e1 != e2:
                if e1.name in listSleepStages and not(e2.name in listSleepStages) :
                    self.stageTransitions.append(transition(e1.timeEnd(), "B"))
                if not(e1.name in listSleepStages) and e2.name in listSleepStages :
                    self.stageTransitions.append(transition(e2.timeStart(), "A"))
                
        if stageEventNames[-1].name in listSleepStages:
            self.stageTransitions.append(transition(stageEventNames[-1].timeEnd(), "B"))                      

        #goldStageEvents = filter(lambda e: e.name in listSleepStages, readerGold.events)  
        #if len(goldStageEvents) :                
        #    self.stageTransitions = concatenate([[transition(event.timeStart(), "A"), 
        #                                          transition(event.timeEnd(), "B")]  for event in goldStageEvents])  
        #else:
        #    self.stageTransitions = []
            
        for channel in listChannels:
            
            goldSpindleTransitions   = getSpindleTransitions(readerGold, nameEventGold, channel, True)
            testedSpindleTransitions = getSpindleTransitions(readerTested, nameEventTested, channel, False)
            transitions = concatenate([self.stageTransitions, goldSpindleTransitions, testedSpindleTransitions])

            if len(goldSpindleTransitions) != 0 or len(testedSpindleTransitions) != 0 :

                msa = StateMachine()
                                           
                # Sort transitions by time of occurence and apply finite state 
                # machine to compute TP, TN, FP, FN
                try:
                    msa.run(sorted(transitions, key=attrgetter('time')))
                except BadTransitionSequence:
                    print "channel:", channel
                    for trans in sorted(transitions, key=attrgetter('time')) :
                        print trans
                    raise                    
                    
                
    
                self.TP[channel] = msa.TP
                self.TN[channel] = msa.TN
                self.FP[channel] = msa.FP
                self.FN[channel] = msa.FN
                  
                  
    def printEvaluation(self, listSleepStages, nameEventGold, 
                          nameEventTested, readerGold, readerTested=None, listChannels=None):

        self.computeStatistics(listSleepStages, nameEventGold, 
                               nameEventTested, readerGold, readerTested, listChannels)
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
        return (self.TN[channel] + self.TP[channel])/self.N[channel]            
        
        
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
        
        
                
        

            
            
            
            