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
 Date  : July 12, 2012

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
 [4] Assous S, Boashash B. Evaluation of the modified S-transform for 
     time-frequency synchrony analysis and source localisation. EURASIP 
     Journal on Advances in Signal Processing, 2012; 2012: 49.     

"""



from math import floor, pi
from scipy import fft, var, exp, ifft, where, transpose, roll

from numpy.fft import fftfreq


###############################################################################
# Computing the Modified S Transform. Using m=0 and k=1 compute the
# S Transform. Translated, adapted and modified from a Matlab code available 
# in the paper [4].
#
def computeMST(sig, fs, m=None, k=None, fmin=None, fmax=None):
    
    N  = len(sig)               # get the length of the signal
    N2 = int(floor(N/2.0)) 


    f = fftfreq(N)     
    MST = []  
    
    SIG = fft(sig, N);          # compute the signal spectrum
    
    if m is None: 
        m   = (1.0/N)
    if k is None:
        k   = 4.0*var(sig) #4.0*var(sig) 
        
    fOut = f[0:N2]*fs        
        
    if fmin is None:
        iMin=1
    else:        
        iMin = where(fOut >= fmin)[0][0]

    if fmax is None:        
        iMax=N2
    else:        
        iMax = where(fOut <= fmax)[0][-1]        
        
    g   = m*f+k                             # parameter gamma
    for i in range(iMin, iMax):
        SIGs = roll(SIG, -(i-1))            # circshift the spectrum SIG
        W = (g[i]/f[i])*2.0*pi*f            # Scale Gaussian
        G = exp((-W**2)/2.0)                # W in Fourier domain
    
        MST.append(ifft(SIGs*G))            # Compute the complex values of MST
    
    return transpose(MST), fOut[iMin:iMax]  

