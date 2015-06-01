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
import time
from spyndle import Filter
from spyndle.miscellaneous import filters

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

        #import matplotlib.pyplot as plt 
        #plt.plot(range(len(filtSig)), abs(filtSig - self.low))
        #plt.plot(range(len(filtSig)), self.low)
        #plt.plot(range(len(filtSig)), filtSig)
        #plt.show()
        #print abs(filtSig - self.low)
        self.assertTrue(np.allclose(filtSig[1000:5000], self.low[1000:5000]))
        

    def testHighPass(self):
        # Defining EEG filters
        highPassFilter = Filter(self.fs)
        highPassFilter.create(low_crit_freq=16, 
                              high_crit_freq=None, order=1001, 
                              btype="highpass", ftype="FIR", useFiltFilt=True)               
        
        filtSig = highPassFilter.applyFilter(self.sig)   
        
        self.assertTrue(np.allclose(filtSig[1000:5000], self.high[1000:5000]))


    def testBandPass(self):
        # Defining EEG filters
        bandPassFilter = Filter(self.fs)
        bandPassFilter.create(low_crit_freq=11.0, 
                              high_crit_freq=16.0, order=1001, 
                              btype="bandpass", ftype="FIR", useFiltFilt=True)               
        
        filtSig = bandPassFilter.applyFilter(self.sig)   
        

        #import matplotlib.pyplot as plt 
        #plt.plot(range(len(filtSig[1000:5000])), abs(filtSig[1000:5000] - self.middle[1000:5000]))
        #plt.plot(range(len(filtSig[1000:5000])), self.middle[1000:5000])
        #plt.plot(range(len(filtSig[1000:5000])), filtSig[1000:5000])
        #plt.show()    
        #print max(np.abs(filtSig[1000:5000] - self.middle[1000:5000]))
        
        self.assertTrue(np.allclose(filtSig[1000:5000], self.middle[1000:5000], rtol=1e-04, atol=1e-05))




    def testNoWeave_filtfilt(self):
        # Defining EEG filters
        bck = filters.filtConf   
        filters.filtConf.useWeave = True 
        
        t1 = time.time()        
        bandPassFilter = Filter(self.fs)
        bandPassFilter.create(low_crit_freq=11.0, 
                              high_crit_freq=16.0, order=1001, 
                              btype="bandpass", ftype="FIR", useFiltFilt=True)               
        
        filtSigWeave = bandPassFilter.applyFilter(self.sig)   
        t2 = time.time()
        
        
        filters.filtConf.useWeave = False        

        t3 = time.time()
        bandPassFilter = Filter(self.fs)
        bandPassFilter.create(low_crit_freq=11.0, 
                              high_crit_freq=16.0, order=1001, 
                              btype="bandpass", ftype="FIR", useFiltFilt=True)               
        
        filtSigNoWeave = bandPassFilter.applyFilter(self.sig)   
        t4 = time.time()
        
        filters.filtConf = bck        


        print(("Time sweave: %f" % (t2 - t1)))   
        print(("Time no sweave: %f" % (t4 - t3)))   

        
        self.assertTrue(np.allclose(filtSigWeave, filtSigNoWeave))






    def testNoWeave(self):
        # Defining EEG filters
        bck = filters.filtConf   
        filters.filtConf.useWeave = True 
        
        t1 = time.time()        
        bandPassFilter = Filter(self.fs)
        bandPassFilter.create(low_crit_freq=11.0, 
                              high_crit_freq=16.0, order=1001, 
                              btype="bandpass", ftype="FIR", useFiltFilt=False)               
        
        filtSigWeave = bandPassFilter.applyFilter(self.sig)   
        t2 = time.time()
        
        
        filters.filtConf.useWeave = False        

        t3 = time.time()
        bandPassFilter = Filter(self.fs)
        bandPassFilter.create(low_crit_freq=11.0, 
                              high_crit_freq=16.0, order=1001, 
                              btype="bandpass", ftype="FIR", useFiltFilt=False)               
        
        filtSigNoWeave = bandPassFilter.applyFilter(self.sig)   
        t4 = time.time()
        
        filters.filtConf = bck        


        print(("Time sweave: %f" % (t2 - t1)))   
        print(("Time no sweave: %f" % (t4 - t3)))   

        
        self.assertTrue(np.allclose(filtSigWeave, filtSigNoWeave))




def main():
    unittest.main()

if __name__ == '__main__':
    main()
