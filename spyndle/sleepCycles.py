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
        self.sleepDeterminingStages = ["Sleep stage 2", "Sleep stage 3", 
                                       "Sleep stage 4", "Sleep stage R", 
                                       "Sleep stage N2", "Sleep stage N3"]
        self.sleepStages            = self.sleepDeterminingStages
        self.minTimeNREM            = 15
        self.minTimeLastNREM        = 5
        self.minTimeREMExceptFirst  = 5
        self.minTimeLastREM         = 5
        self.maxPeriodWithoutREM    = 15
        self.cylesEndWithNextNREM   = True
        self.type                   = "Feinberg & Floyd"



     def setSchulz(self):
        self.minStage1ForSleep      = 0
        self.sleepDeterminingStages = ["Sleep stage 2", "Sleep stage 3", 
                                       "Sleep stage 4", "Sleep stage R", 
                                       "Sleep stage N2", "Sleep stage N3"]
        self.sleepStages            = self.sleepDeterminingStages 
        self.minTimeNREM            = 15
        self.minTimeLastNREM        = 5
        self.minTimeREMExceptFirst  = 0
        self.minTimeLastREM         = 5
        self.maxPeriodWithoutREM    = 15
        self.cylesEndWithNextNREM   = True
        self.type                   = "Schulz"
    

     def setAeschbach(self):
        self.minStage1ForSleep      = 0
        self.sleepDeterminingStages = ["Sleep stage 2", "Sleep stage 3", 
                                       "Sleep stage 4", "Sleep stage R", 
                                       "Sleep stage N2", "Sleep stage N3"]
        self.sleepStages            = self.sleepDeterminingStages 
        self.minTimeNREM            = 15
        self.minTimeLastNREM        = 5
        self.minTimeREMExceptFirst  = 5
        self.minTimeLastREM         = 0
        self.maxPeriodWithoutREM    = 15
        self.cylesEndWithNextNREM   = True
        self.type                   = "Aeschbach"
        

     def setMinimal(self):
        self.minStage1ForSleep      = 1
        self.sleepDeterminingStages = ["Sleep stage 1", "Sleep stage 2", 
                                       "Sleep stage 3", "Sleep stage 4", 
                                       "Sleep stage R", "Sleep stage N1",
                                       "Sleep stage N2", "Sleep stage N3"]
        self.sleepStages            = self.sleepDeterminingStages 
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
        self.timeStartNREM      = 0.0
        self.durationNREM   = 0.0
        self.NREMpages          = []
        
        self.timeStartREM       = 0.0
        self.durationREM    = 0.0
        self.REMpages           = []

        self.completeCycle      = True

    def duration(self):    
        return self.durationNREM + self.durationREM

    def timeStart(self):    
        return self.timeStartNREM
        
    def timeEnd(self):    
        return self.timeStart() + self.duration()


    def __str__(self):
        return ("NREM (start, duration) : (%f, %f) | REM (start, duration) : (%f, %f) " 
                %(self.timeStartNREM, self.durationNREM, self.timeStartREM, self.durationREM))


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
                cycles[i].durationREM    =  cycles[i+1].timeStartNREM  - cycles[i].timeStartREM 
    
        for cycle in cycles:
            cycle.REMpages = filter(lambda e: e.startTime >= cycle.timeStartREM and e.startTime < cycle.timeStartREM + cycle.durationREM, events) 
            cycle.NREMpages = filter(lambda e: e.startTime >= cycle.timeStartNREM and e.startTime < cycle.timeStartNREM + cycle.durationNREM, events) 

    # If no cycle has been found, make a single cycle of the whole recording.
    # This is useful to avoid breaking functions relying in cycles when
    # using shorter recording which full recording time is not enough 
    # for defining a complete cycle. This can happen for example when using
    # a 30-minute recording during which the subject stays in stage sleep 2-4.
    else:
        cycle = DreamCycle()

        startTimes = [e.startTime for e in events]
        endEvent     = filter(lambda e: e.startTime == max(startTimes), events)[0] 

        cycle.timeStartNREM    = min(startTimes)
        cycle.durationNREM = endEvent.startTime + endEvent.timeLength - cycle.timeStartNREM 

        cycles.append(cycle)


    return cycles
    
    
    
def computeUneDreamCyle(events, cycleDefinition, sleeping, noCycle): 
    
    cycle = DreamCycle()
    cycle, lastingEvents, sleeping = computeNREM(cycle, events, cycleDefinition, sleeping, noCycle)
    cycle, lastingEvents, sleeping = computeREM(cycle, lastingEvents, cycleDefinition, sleeping, noCycle)
    
    return cycle, lastingEvents, sleeping    
    
    
    
def computeNREM(cycle, events, cycleDefinition, sleeping, noCycle):
    
    possibleStartTime = None
    startStage1       = None
    for i, event in zip(range(len(events)), events):
        #print sleeping, possibleStartTime, event.name, cycleDefinition.sleepDeterminingStages     
        if not sleeping:        
            if event.name in cycleDefinition.sleepDeterminingStages:
                if event.name == "Sleep stage 1" or  event.name == "Sleep stage N1":
                    if startStage1 is None :
                        startStage1 = event.startTime
                        possibleStartTime = event.startTime
                    elif event.startTime - startStage1 + event.timeLength >= cycleDefinition.minStage1ForSleep*60.0 :
                        sleeping = True
                elif  event.name == "Sleep stage R":
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
                        
                        cycle.timeStartNREM      = possibleStartTime
                        cycle.durationNREM   = event.startTime - possibleStartTime
                        cycle.completeCycle      = False

                        return cycle, events[i:len(events)], sleeping
                else:
                    sleeping = True
                    if possibleStartTime is None :
                        possibleStartTime = event.startTime    

        else:
            if event.name in cycleDefinition.sleepStages:    
                if possibleStartTime is None:
                    possibleStartTime = event.startTime
                
                if  event.name == "Sleep stage R":
                    if noCycle == 1:
                        enterREM = True
                    else:
                        enterREM = goInREM(events, cycleDefinition, i)                    

                    if enterREM:
                        if possibleStartTime is None:
                            possibleStartTime   = event.startTime                  
                        
                        if event.startTime - possibleStartTime >= cycleDefinition.minTimeNREM :                       
                            # On a notre période NREM!
                            cycle.timeStartNREM  = possibleStartTime
                            cycle.durationNREM   = event.startTime - possibleStartTime

                            return cycle, events[i:len(events)], sleeping
                        else:
                            # Période pas assez longue...
                            # Le cycle sera incomplet avec aucune période NREM
                            cycle.timeStartNREM      = possibleStartTime
                            cycle.durationNREM   = event.startTime - possibleStartTime
                            cycle.completeCycle      = False

                            return cycle, events[i:len(events)], sleeping

    
    return None, [], sleeping



def goInREM(events, cycleDefinition, i):
    startTime = events[i].startTime
    
    for event in events[(i+1):len(events)] :
        if event.name !=  "Sleep stage R":
            return event.startTime - startTime >= cycleDefinition.minTimeREMExceptFirst*60.0
    


def computeREM(cycle, events, cycleDefinition, sleeping, noCycle):
 
    if not len(events) : 
        return None, [], sleeping
        
    if events[0].name != "Sleep stage R":
        cycle.completeCycle      = False
        return cycle, events, sleeping
        
    startTime   = events[0].startTime       
        

    for event, i in zip(events, range(len(events))):
        if  event.name == "Sleep stage R":
            finPossibleTime     = None
            finI                = None
        else:
            if finPossibleTime is None:
                finPossibleTime   = event.startTime
                finI              = i
                
            if event.startTime - finPossibleTime + event.timeLength >= cycleDefinition.maxPeriodWithoutREM*60.0:
                # Fin de la période REM
                cycle.timeStartREM       = startTime
                cycle.durationREM    = finPossibleTime - startTime             
                return cycle, events[finI:len(events)], events[finI].name in cycleDefinition.sleepStages  
     
     
     
    # Fin de la nuit    
    if finPossibleTime is None:
        finPossibleTime   = events[-1].startTime + events[-1].timeLength  
        
    cycle.timeStartREM       = startTime
    cycle.durationREM    = finPossibleTime - startTime             
    return cycle, [], False    



