# -*- coding: utf-8 -*-

###############################################################################
# License: For personnal, educationnal, and research purpose, this software is 
#          provided under the Gnu GPL (V.3) license. To use this software in
#          commercial application, please contact the author. 
#
#
# Author: Christian O'Reilly (christian.oreilly@umontreal.ca)
# Date  : June 11, 2013
#
###############################################################################

from scipy import array
import numpy as np
import re
from numpy import sqrt

try:
    from spyndle import Line, Point
except ImportError:
    pass

from spyndle.EEG.electrodesSVG import getElectrodeCoordinates
from spyndle.EEG.mapping import getTransformedCoord


electrodes = ['Fp1', 'Fp2', 'F7', 'F8', 'F3', 'F4', 
              'T3', 'T4', 'C3', 'C4', 'T5', 'T6', 'P3', 
               'P4', 'O1', 'O2', 'Fz', 'Cz', 'Pz', 'Oz']


"""
These positions are taken from Okamoto, M., Dan, H., Sakamoto, K., Takeo, K.,
Shimizu, K., Kohno, S., . . . Dan, I. (2004). Three-dimensional probabilistic
anatomical cranio-cerebral correlation via the international 10-20 system 
oriented for transcranial functional brain mapping. NeuroImage, 21(1), 99-111. 
"""
mni_location_head_surface = {}
mni_location_head_surface["Fp1"] =	[-26.1,	83.5,	-0.1]
mni_location_head_surface["Fp2"] =	[32.1,	81.3,	-0.2]
mni_location_head_surface["Fz"]  =	[-0.1,	53.3,	70.6]
mni_location_head_surface["F3"]  =	[-42.6,	58.0,	39.6]
mni_location_head_surface["F4"]  =	[46.9,	56.7,	40.2]
mni_location_head_surface["F7"]  =	[-68.5,	37.9,	-5.9]
mni_location_head_surface["F8"]  =	[71.4,	35.6,	-7.8]
mni_location_head_surface["Cz"]  =	[0.6,	-12.7,	101.4]
mni_location_head_surface["C3"]  =	[-62.7,	-12.7,	69.8]
mni_location_head_surface["C4"]  =	[64.2,	-15.2,	69.9]
mni_location_head_surface["T3"]  =	[-84.6,	-20.7,	-10.8]
mni_location_head_surface["T4"]  =	[86.0,	-25.5,	-9.4]
mni_location_head_surface["Pz"]  =	[-0.2,	-76.6,	88.5]
mni_location_head_surface["P3"]  =	[-46.8,	-88.2,	58.8]
mni_location_head_surface["P4"]  =	[44.8,	-87.9,	59.6]
mni_location_head_surface["T5"]  =	[-72.2,	-72.4,	0.6]
mni_location_head_surface["T6"]  =	[70.8,	-75.6,	4.1]
mni_location_head_surface["O1"]  =	[-31.9,	-112.6,	17.3]
mni_location_head_surface["O2"]  =	[28.1,	-112.9,	19.3]


def get10_20AdjacentElectrodes():
    electSource = array(["F4","F4","F4","F4","F4","F4",
                        "Fz","Fz","Fz","Fz","Fz","Fz","Fz",
                        "Pz","Pz","Pz","Pz","Pz","Pz","Pz",
                        "P3","P3","P3","P3","P3","P3",
                        "T5","T5","T5","T5",
                        "O1","O1","O1","O1",
                        "F3","F3","F3","F3","F3","F3",
                        "C3","C3","C3","C3","C3","C3","C3","C3",
                        "C4","C4","C4","C4","C4","C4","C4","C4", 
                        "P4","P4","P4","P4","P4","P4",
                        "Cz","Cz","Cz","Cz","Cz","Cz","Cz","Cz",
                        "F7","F7","F7","F7", 
                        "T6","T6","T6","T6",
                        "Fp2","Fp2","Fp2","Fp2",
                        "O2","O2","O2","O2",
                        "F8","F8","F8","F8",
                        "Fp1","Fp1","Fp1","Fp1",
                        "T4","T4","T4","T4","T4",
                        "T3","T3","T3","T3","T3"])
    
    electTest = array(["Fz", "Fp2", "F8", "T4", "C4", "Cz",
                  "Fp1", "Fp2", "F4", "C4", "Cz", "C3", "F3",
                  "C3", "Cz", "C4", "P4", "O2", "O1", "P3",
                  "T3", "C3", "Cz", "Pz", "O1", "T5",
                  "T3", "C3", "P3", "O1",
                  "T5", "P3", "Pz", "O2",
                  "Fz", "Fp1", "F7", "T3", "C3", "Cz",
                  "F7", "F3", "Fz", "Cz", "Pz", "P3", "T5", "T3",
                  "Fz", "F4", "F8", "T4", "T6", "P4", "Pz", "Cz",      
                  "Cz", "C4", "T4", "T6", "O2", "Pz",
                  "F3", "Fz", "F4", "C4", "P4", "Pz", "P3", "C3",
                  "Fp1", "F3", "C3", "T3",
                  "O2", "P4", "C4", "T4",
                  "Fp1", "Fz", "F4", "F8",
                  "O1", "Pz", "P4", "T6",
                  "Fp2", "F4", "C4", "T4",
                  "F7", "F3", "Fz", "Fp2",
                  "T6", "P4", "C4", "F4", "F8",
                  "F7", "F3", "C3", "P3", "T5"])
    
    return electSource, electTest





def get_10_20_electrode_distance(sagital=0.35, transverse=0.35, circonf= 0.55, anterior=0.22, posterior=0.24, parasagitleft=0.24, parasagitright=None):

    if parasagitright is None:
        parasagitright = parasagitleft
        
    distance = {}
    
    for ref, test in zip(("Cz","Cz","Fz", "Pz"), ("Fz","Pz", "Fpz", "Oz")):
        distance[(ref, test)] = sagital*0.2
    for ref, test in zip(("Nz","Iz"), ("Fpz", "Oz")):
        distance[(ref, test)] = sagital*0.1

    for ref, test in zip(("Fp1","Fp1", "T3","T3","O1","O1","T6","T6","F8","F8"), 
                         ("Fp2","F7","F7","T5","T5","O2","O2","T4","T4","Fp2")):
        distance[(ref, test)] = circonf*0.1
    for ref, test in zip(("Fp1","Fp2","O1","O2"), ("Fpz","Fpz","Oz","Oz")):
        distance[(ref, test)] = circonf*0.05    
  

    for ref, test in zip(("Cz", "Cz", "C3", "C4"), ("C3", "C4", "T3", "T4")):
        distance[(ref, test)] = transverse*0.2    
    for ref, test in zip(("Fz", "Fz", "F3", "F4"), ("F3", "F4", "F7", "F8")):
        distance[(ref, test)] = anterior*0.25 
    for ref, test in zip(("Pz", "Pz", "P3", "P4"), ("P3", "P4", "T5", "T6")):
        distance[(ref, test)] = posterior*0.25     
    for ref, test in zip(("C3", "C3", "F3", "P3"), ("P3", "F3", "Fp1", "O1")):
        distance[(ref, test)] = parasagitleft*0.25    
    for ref, test in zip(("C4", "C4", "F4", "P4"), ("P4", "F4", "Fp2", "O2")):
        distance[(ref, test)] = parasagitright*0.25      


    distance[("Cz", "F4")] = 0.5*sqrt((parasagitright*0.25)**2+(transverse*0.2)**2)+\
                             0.5*sqrt((sagital*0.2)**2+(anterior*0.25)**2)
    distance[("C4", "Fz")] = 0.5*sqrt((parasagitright*0.25)**2+(anterior*0.25)**2)+\
                             0.5*sqrt((sagital*0.2)**2+(transverse*0.2)**2)    
    distance[("Cz", "F3")] = 0.5*sqrt((parasagitleft*0.25)**2+(transverse*0.2)**2)+\
                             0.5*sqrt((sagital*0.2)**2+(anterior*0.25)**2)
    distance[("C3", "Fz")] = 0.5*sqrt((parasagitleft*0.25)**2+(anterior*0.25)**2)+\
                             0.5*sqrt((sagital*0.2)**2+(transverse*0.2)**2)    

    distance[("T3", "F3")] = 0.5*sqrt((parasagitleft*0.25)**2+(transverse*0.2)**2)+\
                             0.5*sqrt((circonf*0.1)**2+(anterior*0.25)**2)

    distance[("C3", "F7")] = 0.5*sqrt((parasagitleft*0.25)**2+(anterior*0.25)**2)+\
                             0.5*sqrt((circonf*0.1)**2+(transverse*0.2)**2)

    distance[("T4", "F4")] = 0.5*sqrt((parasagitright*0.25)**2+(transverse*0.2)**2)+\
                             0.5*sqrt((circonf*0.1)**2+(anterior*0.25)**2)

    distance[("C4", "F8")] = 0.5*sqrt((parasagitright*0.25)**2+(anterior*0.25)**2)+\
                             0.5*sqrt((circonf*0.1)**2+(transverse*0.2)**2) 

    
    
    
    
    distance[("Cz", "P4")] = 0.5*sqrt((parasagitright*0.25)**2+(transverse*0.2)**2)+\
                             0.5*sqrt((sagital*0.2)**2+(posterior*0.25)**2)
    distance[("C4", "Pz")] = 0.5*sqrt((parasagitright*0.25)**2+(posterior*0.25)**2)+\
                             0.5*sqrt((sagital*0.2)**2+(transverse*0.2)**2)

    distance[("Cz", "P3")] = 0.5*sqrt((parasagitleft*0.25)**2+(transverse*0.2)**2)+\
                             0.5*sqrt((sagital*0.2)**2+(posterior*0.25)**2)
    distance[("C3", "Pz")] = 0.5*sqrt((parasagitleft*0.25)**2+(posterior*0.25)**2)+\
                             0.5*sqrt((sagital*0.2)**2+(transverse*0.2)**2)
    
    distance[("T3", "P3")] = 0.5*sqrt((parasagitleft*0.25)**2+(transverse*0.2)**2)+\
                             0.5*sqrt((circonf*0.1)**2+(posterior*0.25)**2)
    distance[("C3", "T5")] = 0.5*sqrt((parasagitleft*0.25)**2+(posterior*0.25)**2)+\
                             0.5*sqrt((circonf*0.1)**2+(transverse*0.2)**2)
    distance[("T4", "P4")] = 0.5*sqrt((parasagitright*0.25)**2+(transverse*0.2)**2)+\
                             0.5*sqrt((circonf*0.1)**2+(posterior*0.25)**2)
    distance[("C4", "T6")] = 0.5*sqrt((parasagitright*0.25)**2+(posterior*0.25)**2)+\
                             0.5*sqrt((circonf*0.1)**2+(transverse*0.2)**2)
    
    
    
    
    
    
    distance[("Fz", "Fp2")] = 0.5*sqrt((parasagitright*0.25)**2+(anterior*0.25)**2)+\
                              0.5*sqrt((sagital*0.2)**2+(circonf*0.05)**2)
    distance[("F4", "Fpz")] = 0.5*sqrt((parasagitright*0.25)**2+(circonf*0.05)**2)+\
                              0.5*sqrt((sagital*0.2)**2+(anterior*0.25)**2)
    
    distance[("Fz", "Fp1")] = 0.5*sqrt((parasagitleft*0.25)**2+(anterior*0.25)**2)+\
                              0.5*sqrt((sagital*0.2)**2+(circonf*0.05)**2)
    distance[("F3", "Fpz")] = 0.5*sqrt((parasagitleft*0.25)**2+(circonf*0.05)**2)+\
                              0.5*sqrt((sagital*0.2)**2+(anterior*0.25)**2)
    
    
    
    
    distance[("Pz", "O2")] = 0.5*sqrt((parasagitright*0.25)**2+(posterior*0.25)**2)+\
                             0.5*sqrt((sagital*0.2)**2+(circonf*0.05)**2)
    distance[("P4", "Oz")] = 0.5*sqrt((parasagitright*0.25)**2+(circonf*0.05)**2)+\
                             0.5*sqrt((sagital*0.2)**2+(posterior*0.25)**2)
    distance[("Pz", "O1")] = 0.5*sqrt((parasagitleft*0.25)**2+(posterior*0.25)**2)+\
                             0.5*sqrt((sagital*0.2)**2+(circonf*0.05)**2)
    distance[("P3", "Oz")] = 0.5*sqrt((parasagitleft*0.25)**2+(circonf*0.05)**2)+\
                             0.5*sqrt((sagital*0.2)**2+(posterior*0.25)**2)
    
    
    for ref in ("Fp1", "Fp2"):
        distance[(ref, "Nz")] = sqrt((sagital*0.1)**2+(circonf*0.05)**2)

    for ref in ("O1", "O2"):
        distance[(ref, "Iz")] = sqrt((sagital*0.1)**2+(circonf*0.05)**2)
       
    keys = list(distance.keys())
    for ref, test in keys:
        distance[(test, ref)] = distance[(ref, test)]
    
    return distance



"""
 This function tries to determinate pairs of adjacents electrodes from a subset
 of the electrodes of the 10-20 system (see the eeg_electrodes_10-20_small.svg
 for the schema of the electrodes and of their placement)
"""
def getAdjacentsPairs_10_20(electrodes):
    x, y, transform = getElectrodeCoordinates(electrodes)
    for elect in electrodes:
        x[elect], y[elect] = tuple(getTransformedCoord(x[elect], y[elect], transform[elect]))
    
    NbIntersectMax = 0
    electRef  = []
    electTest = []
    for key1 in electrodes:
        for key2 in [key for key in electrodes if key not in [key1]]:
            l1 = Line(Point(x[key1], y[key1]), Point(x[key2], y[key2]))        
            
            nbIntersect = 0
            for key3 in [key for key in electrodes if not key in [key1, key2]]:
                for key4 in [key for key in electrodes if not key in [key1, key2, key3]]:
                    l2 = Line(Point(x[key3], y[key3]), Point(x[key4], y[key4])) 
                      
                    if l1.isIntersecting(l2):
                        if l1.length() > 1.2*l2.length():
                            nbIntersect += 1
                    if nbIntersect > NbIntersectMax:
                        break
                if nbIntersect > NbIntersectMax:
                    break
                
            if nbIntersect <= NbIntersectMax:
                electRef.append(key1)
                electTest.append(key2)
        
    
    return list(zip(electRef, electTest))
    
    
    
    


def get_10_20_1electrodeDistant_electrodes(adjElect):

    results = []
    
    refs  = np.array([ref  for (ref, test) in adjElect])
    tests = np.array([test for (ref, test) in adjElect]) 
    for ref in np.unique(refs) :
        nextAdjs = []
        for adj in tests[refs == ref] :
            nextAdjs.extend(tests[refs == adj])

        nextAdjs = np.unique(nextAdjs)
        nextAdjs = nextAdjs[np.logical_not(np.in1d(nextAdjs, tests[refs == ref]))]
        nextAdjs = [nextAdj for nextAdj in  nextAdjs if nextAdj != ref]
        for test in nextAdjs:
            results.append((ref, test))

    return results



def get_10_20_2electrodeDistant_electrodes(adjElect):

    results = []
    nextAdjElect = get_10_20_1electrodeDistant_electrodes(adjElect)
        
    refs  = np.array([ref  for (ref, test) in adjElect])
    tests = np.array([test for (ref, test) in adjElect])
    nextRefs  = np.array([ref  for (ref, test) in nextAdjElect])
    nextTests = np.array([test for (ref, test) in nextAdjElect])
    
    for ref in np.unique(nextRefs) :

        secNextAdj = []
        for nextAdj in nextTests[nextRefs == ref]:
            secNextAdj.extend(tests[refs == nextAdj])

        secNextAdj = np.unique(secNextAdj)
        secNextAdj = secNextAdj[np.logical_not(np.in1d(secNextAdj, nextTests[nextRefs == ref]))]
        secNextAdj = secNextAdj[np.logical_not(np.in1d(secNextAdj, tests[refs == ref]))]
        secNextAdj = [secNext for secNext in secNextAdj if secNext != ref]
        
        if len(secNextAdj) :
            for secNext in secNextAdj :
                results.append((ref, secNext))

    return results



def get_10_20_3electrodeDistant_electrodes(adjElect):

    results = []
    secNextAdjElect = get_10_20_2electrodeDistant_electrodes(adjElect)
    nextAdjElect    = get_10_20_1electrodeDistant_electrodes(adjElect)
        
    refs         = np.array([ref  for (ref, test) in adjElect])
    tests        = np.array([test for (ref, test) in adjElect])
    nextRefs     = np.array([ref  for (ref, test) in nextAdjElect])
    nextTests    = np.array([test for (ref, test) in nextAdjElect])
    secNextRefs  = np.array([ref  for (ref, test) in secNextAdjElect])
    secNextTests = np.array([test for (ref, test) in secNextAdjElect])    

    for ref in np.unique(secNextRefs) :

        thirdNextAdj = []
        for secNextAdj in secNextTests[secNextRefs == ref] :
            thirdNextAdj.extend(nextTests[nextRefs == secNextAdj])
 
        thirdNextAdj = np.unique(thirdNextAdj)
        
        thirdNextAdj = thirdNextAdj[np.logical_not(np.in1d(thirdNextAdj, secNextTests[secNextRefs == ref]))]
        thirdNextAdj = thirdNextAdj[np.logical_not(np.in1d(thirdNextAdj, nextTests[nextRefs == ref]))]
        thirdNextAdj = thirdNextAdj[np.logical_not(np.in1d(thirdNextAdj, tests[refs == ref]))]
        thirdNextAdj = [thirdNext for thirdNext in thirdNextAdj if thirdNext != ref]
        
        if len(thirdNextAdj):
            for thirdNext in thirdNextAdj :            
                results.append((ref, thirdNext))

    return results
   
    
    







"""
 Considering an EEG channel as constituted of two electrodes, a passive
 (the reference) and an active (the other), this function return the 
 active electrode name from the channel label. This is used, for example
 to get 'F3' from 'F3-ref'. The function tries to 
 find in channelLabel the presence of the strings electrodeNames
"""
def getActiveElectrode(channelLabel, electrodeNames=electrodes, excludePatterns=[], caseSensitive=False):
    
    retName = ""
    for electrodeName in electrodeNames:
        if caseSensitive:
            if electrodeName in channelLabel :
                if retName == "":
                    retName = electrodeName
                else:
                    raise "Conficting electrode names."

        else :
            if electrodeName.lower() in channelLabel.lower():
                if retName == "":
                    retName = electrodeName
                else:
                    raise "Conficting electrode names."


    for pattern in excludePatterns:
        if pattern in channelLabel :
            return ""
            
    return retName

def getActiveElectrodes(channelLabels, electrodeNames=electrodes, excludePatterns=[], caseSensitive=False):
    return [getActiveElectrode(channelLabel, electrodeNames, excludePatterns, caseSensitive) 
                                                                for channelLabel in channelLabels]




"""
 Return the subset of channels which are EEG channels of the 10-20 system.
"""
def getEEGChannels(channelLabels, electrodeNames=electrodes, excludePatterns=[]):
    
    OK = len(channelLabels)*[False]
    for ind, channel in enumerate(channelLabels):
        if len([include for include in electrodeNames if include.lower() in channel.lower()]):
            OK[ind] = True
        if len([include for include in excludePatterns if include in channel]):
            OK[ind] = False                
            
    return [channel for ind, channel in enumerate(channelLabels) if OK[ind]]
    





"""
 Return the subset of channels which are EEG channels of the 10-20 system.
"""
def getPropagationDirection(pairs):
    
    def renameElectrode(electrode):    
        if electrode == "Fp1":
            electrode = "Fp3"
        elif electrode == "FT7":
            electrode = "FC7"
        elif electrode == "FT9":
            electrode = "FC9"
        elif electrode == "T7":
            electrode = "C7"
        elif electrode == "T9":
            electrode = "C9"
        elif electrode == "A1":
            electrode = "C11"
        elif electrode == "TP7":
            electrode = "CP7"
        elif electrode == "TP9":
            electrode = "CP9"
        elif electrode == "O1":
            electrode = "O3"
        elif electrode == "T3":
            electrode = "C7"
        elif electrode == "T5":
            electrode = "P7"
                
        elif electrode == "Fp2":
            electrode = "Fp4"
        elif electrode == "FT8":
            electrode = "FC8"
        elif electrode == "FT10":
            electrode = "FC10"
        elif electrode == "T8":
            electrode = "C8"
        elif electrode == "T10":
            electrode = "C10"
        elif electrode == "A2":
            electrode = "C12"
        elif electrode == "TP8":
            electrode = "CP8"
        elif electrode == "TP10":
            electrode = "CP10"
        elif electrode == "O2":
            electrode = "O4"
        elif electrode == "T4":
            electrode = "C8"
        elif electrode == "T6":
            electrode = "P8"

        return electrode


    def splitElectrodeName(electrode):    
            
            if electrode[-1] == "z" :
                return electrode[:-1], 0
            else :
                reMatch = re.match(r'(?P<first>[A-z]*)(?P<second>[0-9]*)' , electrode)
                return reMatch.group("first"), int(reMatch.group("second"))



    def getDirections(pair):
        ref1, ref2   = splitElectrodeName(renameElectrode(pair[0]))
        test1, test2 = splitElectrodeName(renameElectrode(pair[1]))
        
        rowOrder    = ["I", "O", "PO", "P", "CP", "C", "FC", "F", "AF", "Fp", "N"]
        colOrder    = list(range(11, 0, -2)) + list(range(0, 13, 2))
        towardFront = rowOrder.index(test1) - rowOrder.index(ref1)
        towardRight = colOrder.index(test2) - colOrder.index(ref2)
            
        transEmisphere = False
        towardMedial   = np.nan                      
        if ref2 == 0:
            towardMedial   = -abs(towardRight)       
        elif test2 == 0:           
            towardMedial   = abs(towardRight)  
        elif np.mod(ref2, 2) == np.mod(test2, 2):
            transEmisphere = True            
        elif np.mod(ref2, 2) :
            towardMedial   = towardRight
        else:
            towardMedial   = -towardRight                    
                
        return towardFront, towardRight, towardMedial, transEmisphere

    return [pair + getDirections(pair) for pair in pairs]
