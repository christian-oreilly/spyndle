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
            


# Class representing a detected spindle.
class DetectorEvaluator:
    def __init__(self, goldStandard, tested):

        self.goldStandard   = goldStandard
        self.tested         = tested
        
        self.N  = {}
        self.TP = {}
        self.TN = {}
        self.FP = {}
        self.FN = {}

    
    def printEvaluation(self, listChannels, listSleepStages):
        self.computeStatistics(listChannels, listSleepStages)
        for channel in listChannels:            
            print ("Channel:%s, sensitivity=%f, specificity=%f, PPV=%f, NPV=%f\n"
                   "            TP=%d, TN=%d, FP=%d, FN=%d" % (channel, 
                   self.sensitivity(channel), self.specificity(channel),
                   self.PPV(channel), self.NPV(channel), self.TP[channel], 
                   self.TN[channel], self.FP[channel], self.FN[channel]))
    
    
    
    
    
    def computeStatistics(self, listChannels, listSleepStages):
    
        print "Detecting spindles using the gold standard..."
        self.goldStandard.setDetectionStages(listSleepStages)
        self.goldStandard.detectSpindles(listChannels) 

        print "Detecting spindles using the tested detector..."
        self.tested.setDetectionStages(listSleepStages)
        self.tested.detectSpindles(listChannels) 

        # We want the stage identification to be case-insensitive
        for i in range(len(listSleepStages)): 
            listSleepStages[i] = listSleepStages[i].lower()


        # computing arrays of sample index associated to every sleep stages for
        # which we want to evaluate the spindle detection.
        goldStageEvents = filter(lambda e: e.name.lower() in listSleepStages, self.goldStandard.reader.events)  
        if len(goldStageEvents) :                
            self.stageTransitions = concatenate([[transition(event.timeStart(), "A"), 
                                                  transition(event.timeEnd(), "B")]  for event in goldStageEvents])  
        else:
            self.stageTransitions = []
         
         
        def getSpindleTransitions(detector, channel, isGold):
            channelSpindles = filter(lambda s: s.channel == channel, detector.detectedSpindles)     
            if len(channelSpindles) :
                if isGold :
                    return concatenate([[transition(spindle.startTime(), "C"), 
                                              transition(spindle.endTime(), "D")]  for spindle in channelSpindles])   
                else:
                    return concatenate([[transition(spindle.startTime(), "E"), 
                                              transition(spindle.endTime(), "F")]  for spindle in channelSpindles])              
            return []
            
            
            
        for channel in listChannels:
            transitions = concatenate([self.stageTransitions, 
                                       getSpindleTransitions(self.goldStandard, channel, True),
                                       getSpindleTransitions(self.tested, channel, False)])

            msa = StateMachine()
                                       
            # Sort transitions by time of occurence and apply finite state 
            # machine to compute TP, TN, FP, FN
            msa.run(sorted(transitions, key=attrgetter('time')))
            

            self.TP[channel] = msa.TP
            self.TN[channel] = msa.TN
            self.FP[channel] = msa.FP
            self.FN[channel] = msa.FN
                  
                  
                  
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
        
        
                
        

            
            
            
            