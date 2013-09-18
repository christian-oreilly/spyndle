# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 16:50:26 2013

@author: oreichri
"""



"""
 NOTE: Border effects cause these test to fail. This situation might by
       resolved by using fully adaptative IFIR filters with removed border
       effect (Batista, E. L. O., Tobias, O. J., & Seara, R. (2008). A fully 
       adaptive IFIR filter with removed border effect. Paper presented at the 
       Acoustics, Speech and Signal Processing, 2008. ICASSP 2008. IEEE 
       International Conference on.).
"""


import numpy as np
import unittest
from spyndle import Filter

class FilterTests(unittest.TestCase) :
    
    def setUp(self):
        self.fs     = 200
        self.t      = np.arange(0, 30, 1.0/self.fs)
        self.middle = np.sin(2*np.pi*13.0*self.t) 
        self.low    = np.sin(2*np.pi*5.0*self.t)
        self.high   = np.sin(2*np.pi*25.0*self.t)
        self.sig = self.middle + self.low + self.high 


    def testLowPass(self):
        # Defining EEG filters
        lowPassFilter = Filter(self.fs)
        lowPassFilter.create(low_crit_freq=None, 
                              high_crit_freq=10.0, order=1001, 
                              btype="lowpass", ftype="FIR", useFiltFilt=True)           
                
        filtSig = lowPassFilter.applyFilter(self.sig)   

        import matplotlib.pyplot as plt 
        plt.plot(range(len(filtSig)), abs(filtSig - self.low))
        plt.plot(range(len(filtSig)), self.low)
        plt.plot(range(len(filtSig)), filtSig)
        plt.show()

        
        print abs(filtSig - self.low)
        self.assertTrue(np.allclose(filtSig, self.low))
        

    def testHighPass(self):
        # Defining EEG filters
        highPassFilter = Filter(self.fs)
        highPassFilter.create(low_crit_freq=16, 
                              high_crit_freq=None, order=1001, 
                              btype="highpass", ftype="FIR", useFiltFilt=True)               
        
        filtSig = highPassFilter.applyFilter(self.high)   
        
        self.assertTrue(np.allclose(filtSig, self.sig))


    def testBandPass(self):
        # Defining EEG filters
        bandPassFilter = Filter(self.fs)
        bandPassFilter.create(low_crit_freq=11.0, 
                              high_crit_freq=16.0, order=1001, 
                              btype="highpass", ftype="FIR", useFiltFilt=True)               
        
        filtSig = bandPassFilter.applyFilter(self.middle)   
        
        self.assertTrue(np.allclose(filtSig, self.sig))





def main():
    unittest.main()

if __name__ == '__main__':
    main()
