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


def get10_20AdjacentElectrodes():
    electSource = ["F4","F4","F4","F4","F4","F4",
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
                "T3","T3","T3","T3","T3"]
    
    electTest = ["Fz", "Fp2", "F8", "T4", "C4", "Cz",
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
                  "F7", "F3", "C3", "P3", "T5"]
    
    return electSource, electTest

