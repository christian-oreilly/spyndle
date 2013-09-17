# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 16:50:26 2013

@author: oreichri
"""


import numpy as np
from spyndle.miscellaneous import Filter
t   = np.arange(0, 30, 1.0/200.0)
sig = np.sin(2*np.pi*13.0*t) + np.sin(2*np.pi*5.0*t) + np.sin(2*np.pi*25.0*t)
# Defining EEG filters
bandPassFilter = Filter(200.0)
bandPassFilter.create(low_crit_freq=11.0, 
                      high_crit_freq=16.0, order=1001, 
                      btype="bandpass", ftype="FIR", useFiltFilt=True)          

sig2 = bandPassFilter.applyFilter(sig)     

from matplotlib import pyplot as plt
plt.figure()
plt.plot(range(len(sig)), sig)
plt.plot(range(len(sig)), sig2)
plt.show()   
   
   
print np.sqrt(np.mean(sig2**2))