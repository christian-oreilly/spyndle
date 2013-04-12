# -*- coding: utf-8 -*-

"""
    Code assessing transient event propogation through an array of sensors
    using cross-correlation of S-transform of the signal captured by
    the different sensors.

    Copyright (C) 2012-2013  Christian O'Reilly

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.


    For personnal, educationnal, and research purpose, this software is 
    provided under the GNU GPL (V.3) license: you can redistribute it and/or
    modify it under the terms of the version 3 of the GNU General Public 
    License as published by the Free Software Foundation.
          
    To use this software in commercial application, please contact the author. 
    If used for research purpose, the reference [1] or references [2-3] should 
    be cited in the derived publication to refere the reader to the description 
    of the methodology.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


 Author: Christian O'Reilly (christian.oreilly@umontreal.ca)
 Date  : June 27, 2012

 [1] O’Reilly, C., & Nielsen, T. (2013). Assessing the propagation of EEG 
     transient activity, Proceedings of the 8th International Workshop on 
     Systems, Signal Processing and their Applications, Algiers, Algeria, 
     12-15 May 2013.
 [2] O'Reilly, C. & Nielsen, T. Assessing EEG sleep spindle propagation. 
     Part 1: Theory and proposed methodology, Submitted to Journal of 
     Neuroscience Methods, april 2013.
 [3] O'Reilly, C. & Nielsen, T. Assessing EEG sleep spindle propagation. 
     Part 2: Experimental characterization, Submitted to Journal of 
     Neuroscience Methods, april 2013.     

"""

from copy import deepcopy

class cycleDefinitions():
     def __init__(self):
        self.setMinimal()


     def setFeinberg(self):
        self.minStage1ForSleep      = 0
        self.sleepDeterminingStages = ["Stage2", "Stage3", "Stage4", "REM"]
        self.sleepStages            = ["Stage2", "Stage3", "Stage4", "REM"]
        self.minTimeNREM            = 15
        self.minTimeLastNREM        = 5
        self.minTimeREMExceptFirst  = 5
        self.minTimeLastREM         = 5
        self.maxPeriodWithoutREM    = 15
        self.cylesEndWithNextNREM   = True
        self.type                   = "Feinberg & Floyd"


     def setSchulz(self):
        self.minStage1ForSleep      = 0
        self.sleepDeterminingStages = ["Stage2", "Stage3", "Stage4", "REM"]
        self.sleepStages            = ["Stage2", "Stage3", "Stage4", "REM"]
        self.minTimeNREM            = 15
        self.minTimeLastNREM        = 5
        self.minTimeREMExceptFirst  = 0
        self.minTimeLastREM         = 5
        self.maxPeriodWithoutREM    = 15
        self.cylesEndWithNextNREM   = True
        self.type                   = "Schulz"
    

     def setAeschbach(self):
        self.minStage1ForSleep      = 0
        self.sleepDeterminingStages = ["Stage2", "Stage3", "Stage4", "REM"]
        self.sleepStages            = ["Stage2", "Stage3", "Stage4", "REM"]
        self.minTimeNREM            = 15
        self.minTimeLastNREM        = 5
        self.minTimeREMExceptFirst  = 5
        self.minTimeLastREM         = 0
        self.maxPeriodWithoutREM    = 15
        self.cylesEndWithNextNREM   = True
        self.type                   = "Aeschbach"
        

     def setMinimal(self):
        self.minStage1ForSleep      = 1
        self.sleepDeterminingStages = ["Stage1", "Stage2", "Stage3", "Stage4", "REM"]
        self.sleepStages            = ["Stage1", "Stage2", "Stage3", "Stage4", "REM"]
        self.minTimeNREM            = 0
        self.minTimeLastNREM        = 0
        self.minTimeREMExceptFirst  = 0
        self.minTimeLastREM         = 0
        self.maxPeriodWithoutREM    = 15
        self.cylesEndWithNextNREM   = False
        self.type                   = "Minimal"
    
    
    
    
    
    
    
    
    
# Dream cycles are composed of one NREM period follow by a REM period.
class DreamCycle :
    def __init__(self):    
        self.sampleStartNREM    = 0.0
        self.timeStartNREM      = 0.0
        self.timeDurationNREM   = 0.0
        self.sampleDurationNREM = 0.0
        self.NREMpages          = []
        
        self.sampleStartREM     = 0.0
        self.timeStartREM       = 0.0
        self.timeDurationREM    = 0.0
        self.sampleDurationREM  = 0.0
        self.REMpages           = []

        self.completeCycle      = True

    def sampleDuration(self):    
        return self.sampleDurationNREM + self.sampleDurationREM

    def sampleStart(self):    
        return self.sampleStartNREM
        
    def sampleEnd(self):    
        return self.sampleStart() + self.sampleDuration()





def computeDreamCycles(events, cyclesDefinition):
    lastingEvents = deepcopy(events)
    cycles = []
    sleeping = False
    while(len(lastingEvents)):
        cycle, lastingEvents, sleeping = computeUneDreamCyle(lastingEvents, cyclesDefinition, sleeping, len(cycles)+1) 
        if not cycle is None : 
            cycles.append(cycle)

    if len(cycles):
        # Vérifier rêgles pour le dernier cycle
        if (cycles[-1].timeDurationREM   < cyclesDefinition.minTimeLastREM*60.0  or 
            cycles[-1].timeDurationNREM  < cyclesDefinition.minTimeLastNREM*60.0 ) :
               cycles = cycles[0:-1]
    
    
        # Gérer "mettre fin au début de la période NREM suivante"
        if cyclesDefinition.cylesEndWithNextNREM :
            for i in range(len(cycles)-1):
                cycles[i].timeDurationREM    =  cycles[i+1].timeStartNREM  - cycles[i].timeStartREM 
                cycles[i].sampleDurationREM  =  cycles[i+1].sampleStartNREM  - cycles[i].sampleStartREM 
    
        for cycle in cycles:
            cycle.REMpages = filter(lambda e: e.startSample >= cycle.sampleStartREM and e.startSample < cycle.sampleStartREM + cycle.sampleDurationREM, events) 
            cycle.NREMpages = filter(lambda e: e.startSample >= cycle.sampleStartNREM and e.startSample < cycle.sampleStartNREM + cycle.sampleDurationNREM, events) 

    return cycles
    
def computeUneDreamCyle(events, cycleDefinition, sleeping, noCycle): 
    
    cycle = DreamCycle()
    cycle, lastingEvents, sleeping = computeNREM(cycle, events, cycleDefinition, sleeping, noCycle)
    cycle, lastingEvents, sleeping = computeREM(cycle, lastingEvents, cycleDefinition, sleeping, noCycle)
    
    return cycle, lastingEvents, sleeping    
    
    
def computeNREM(cycle, events, cycleDefinition, sleeping, noCycle):
    
    possibleStartTime = None
    possibleStartPage = None
    startStage1       = None
    for i, event in zip(range(len(events)), events):
        #print sleeping, possibleStartTime, event.name, cycleDefinition.sleepDeterminingStages     
        if not sleeping:               
            if event.name in cycleDefinition.sleepDeterminingStages:
                if event.name == "Stage1":
                    if startStage1 is None :
                        startStage1 = event.startTime
                        possibleStartTime = event.startTime
                        possibleStartPage = event.startSample
                    elif event.startTime - startStage1 + event.timeLength >= cycleDefinition.minStage1ForSleep*60.0 :
                        sleeping = True
                elif event.name == "REM":
                    sleeping = True
                    if noCycle == 1:
                        enterREM = True
                    else:
                        enterREM = goInREM(events, cycleDefinition, i)                    

                    if enterREM:                    
                        # Période pas assez longue...
                        # Le cycle sera incomplet avec aucune période NREM
                        if possibleStartTime is None:
                            possibleStartTime   = event.startTime 
                            possibleStartPage   = event.startSample 
                        
                        cycle.sampleStartNREM    = possibleStartPage
                        cycle.timeStartNREM      = possibleStartTime
                        cycle.timeDurationNREM   = event.startTime - possibleStartTime
                        cycle.sampleDurationNREM = event.startSample - possibleStartPage
                        cycle.completeCycle      = False

                        return cycle, events[i:len(events)], sleeping
                else:
                    sleeping = True
                    if possibleStartTime is None :
                        possibleStartTime = event.startTime
                        possibleStartPage = event.startSample     

        else:
            if event.name in cycleDefinition.sleepStages:    
                if possibleStartTime is None:
                    possibleStartTime = event.startTime
                    possibleStartPage = event.startSample    
                
                if event.name == "REM":
                    if noCycle == 1:
                        enterREM = True
                    else:
                        enterREM = goInREM(events, cycleDefinition, i)                    

                    if enterREM:
                        if possibleStartTime is None:
                            possibleStartTime   = event.startTime 
                            possibleStartPage   = event.startSample                         
                        
                        if event.startTime - possibleStartTime >= cycleDefinition.minTimeNREM :                       
                            # On a notre période NREM!
                            cycle.sampleStartNREM    = possibleStartPage
                            cycle.timeStartNREM      = possibleStartTime
                            cycle.timeDurationNREM   = event.startTime - possibleStartTime
                            cycle.sampleDurationNREM = event.startSample - possibleStartPage

                            return cycle, events[i:len(events)], sleeping
                        else:
                            # Période pas assez longue...
                            # Le cycle sera incomplet avec aucune période NREM
                            cycle.sampleStartNREM    = possibleStartPage
                            cycle.timeStartNREM      = possibleStartTime
                            cycle.timeDurationNREM   = event.startTime - possibleStartTime
                            cycle.sampleDurationNREM = event.startSample - possibleStartPage
                            cycle.completeCycle      = False

                            return cycle, events[i:len(events)], sleeping

    
    return None, [], sleeping



def goInREM(events, cycleDefinition, i):
    startTime = events[i].startTime
    
    for event in events[(i+1):len(events)] :
        if event.name != "REM":
            return event.startTime - startTime >= cycleDefinition.minTimeREMExceptFirst*60.0
    


def computeREM(cycle, events, cycleDefinition, sleeping, noCycle):
 
    if not len(events) : 
        return None, [], sleeping
        
    if events[0].name != "REM":
        cycle.completeCycle      = False
        return cycle, events, sleeping
        
    startTime   = events[0].startTime        
    startSample = events[0].startSample       
        

    for event, i in zip(events, range(len(events))):
        if event.name == "REM":
            finPossibleTime     = None
            finPossibleSample   = None
            finI                = None
        else:
            if finPossibleTime is None:
                finPossibleTime   = event.startTime
                finPossibleSample = event.startSample
                finI              = i
                
            if event.startTime - finPossibleTime + event.timeLength >= cycleDefinition.maxPeriodWithoutREM*60.0:
                # Fin de la période REM
                cycle.sampleStartREM     = startSample
                cycle.timeStartREM       = startTime
                cycle.timeDurationREM    = finPossibleTime - startTime
                cycle.sampleDurationREM  = finPossibleSample - startSample              
                return cycle, events[finI:len(events)], events[finI].name in cycleDefinition.sleepStages  
     
     
     
    # Fin de la nuit    
    if finPossibleTime is None:
        finPossibleTime   = events[-1].startTime + events[-1].timeLength  
        finPossibleSample = events[-1].startSample + events[-1].sampleLength 
        
    cycle.sampleStartREM     = startSample
    cycle.timeStartREM       = startTime
    cycle.timeDurationREM    = finPossibleTime - startTime
    cycle.sampleDurationREM  = finPossibleSample - startSample              
    return cycle, [], False    



