

# -*- coding: utf-8 -*-

import os,sys
parentdir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,parentdir + "\\..\\..\\..") 


from spyndle.propagation import computeSPF

import unittest

class TestAnalysis(unittest.TestCase) :
    
    
    def testAnalysisPropagation(self):
        electrodes = ["O1", "O2", "T5", "P3", "Pz", "P4", "T6", "T3", "C3", "Cz", "C4", 
                     "T4", "F7", "F3", "Fz", "F4", "F8", "Fp1", "Fp2"]
        offFiles = ["offset_" + e + ".txt" for e in electrodes]
        syncFiles = ["sync_" + e + ".txt" for e in electrodes]
        computeSPF(parentdir + "\\", "test", syncFiles, offFiles, 
                    electrodes, electrodes, alpha = 10.0, verbose = True)
                   
    
    
        #self.failIf(SNR < 50.0)
    

def main():
    unittest.main()

if __name__ == '__main__':
    main()






