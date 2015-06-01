# -*- coding: utf-8 -*-

"""
 Pure python implementation of the fast ST transformed as proposed in [1].This
 code inspired from the PyGFT package proposed by these authors. It is 
 most probably signigicantly slower, but it dont need the installation of the
 FFTw library, it is all in pure python, it is well integrated with numpy/scipy,
 and it is much faster than the discrete ST transformed.
 
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
    If used for research purpose, the reference [1] should be cited in the 
    whenever appropriate.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


 Author: Christian O'Reilly (christian.oreilly@umontreal.ca)
 Date  : June 7, 2013

 [1]  Brown, R.A., Lauzon, M.L., and Frayne, R., “A General Description of 
      Linear Time-Frequency Transforms and Formulation of a Fast, Invertible 
      Transform That Samples the Continuous S-Transform Spectrum Nonredundantly,” 
      IEEE Trans. Signal Process. vol. 58, pp. 281-290, Jan 2010    

"""



import numpy
from numpy import abs, round, meshgrid, vectorize
from scipy import sqrt, pi, exp, arange, complex128
from scipy import log2, zeros, concatenate, fft, ifft, transpose
from scipy import where
from scipy.interpolate import griddata



from numpy.fft import fftfreq



def windows(N, windowFct):

    win = zeros(N, dtype=complex128)
    	
    # For each of the GFT frequency bands.  Using a dyadic scale.

    # Frequency 0 and -1 are special cases
    win[0]  = 1.0;
    win[-1] = 1.0;
	
    fstart = 1
    while fstart < N/2 :
        if fstart < 2:
             fwidth = 1
             fcentre = fstart +1
        else: 
             fwidth = fstart
             fcentre = fstart + fwidth/2 +1
		
        # Construct the window (in the Fourier domain)
        temp = roll(windowFct(N, fcentre), fwidth/2)
  		
        # Put the window in the proper +- partitions
        for i in range(fwidth):
            win[fstart+i]       = temp[i];
            win[N-fstart-1-i]   = temp[i];

        fstart *= 2 

    return win


"""
 Cannot use the scipy roll as the "take" function seem to complain about
 complex numbers.
"""
def roll(a, shift):
    shift = -int(round(shift))
    return concatenate((a[shift:],a[0:shift]))



def gaussian(N, freq):
    
    # Make a gaussian between 0 and 1 centered at 0.5
    x = arange(0,1,1.0/N).astype(complex128)
    win = abs(freq)/sqrt(2*pi)*exp(-(x-0.5)**2*abs(freq)**2/2)
     
    # Normalize its area
    win = win/sum(win) #trapz(win, dx=1.0/N)

    return fft(roll(win, -N/2.0), N)





def diadicPartitions(N):
    sf = 1
    cf = 1
    width = 1
    pcount = 0
    partitions = zeros(round(log2(N))*2, dtype=numpy.int)
    pOff       = round(log2(N))*2-1
	
    while sf < N/2 :
        ep = cf+width/2-1
        sn = N-cf-width/2+1; en = N-cf+width/2+1
        if ep > N :
            ep = N
        if sn < 0:
            sn = 0
        if width/2 == 0 :
            ep+=1; sn-=1
        partitions[pcount] = ep
        partitions[pOff-pcount] = en
        pcount += 1
		
        sf = sf+width;
        if (sf > 2) :
            width *= 2;
        
        cf = sf+width/2;

    return partitions





def gft1d(signal, windowType):

    N = len(signal) 
    
    # Compute the windows signal
    if (windowType == 'gaussian'):
        windowFct = gaussian
    else:
        raise NotImplemented
        #windowFct = &box    
    
    win = windows(N, windowFct)    
  
    # Do the initial FFT of the signal   
    signal = fft(signal, N)

    # Apply the windows	
    signal = signal*win

    # For each of the GFT frequency bands
    fstart = 0;
    for fend in diadicPartitions(N) : 

        # frequency band that we're working with : signal[fstart*2:(fstart+fend)]
        # inverse FFT to transform to S-space
        signal[fstart:fend] = ifft(signal[fstart:fend], fend-fstart)

        fstart = fend
 
    return roll(signal,N/2)  




def dyadic2D(N):   
    # For each of the GFT frequency bands
    fstart = 0;
    F = zeros(N)
    T = zeros(N)
    for fend in diadicPartitions(N) : 
        F[fstart:fend] = (fstart+fend)/2.0
        T[fstart:fend] = ((arange(0, 1, 1.0/(fend-fstart))+0.5/(fend-fstart))*N).astype(int)
        fstart = fend
    return roll(F,N/2), roll(T,N/2)





def computeFastST(sig, fs, fmin=0, fmax=None):

    if not fmax:
        fmax = fs/2.0

    N  = len(sig)               # get the length of the signal

    f    = fftfreq(N)     
    fOut = f*fs    
    indF = where((fOut >= fmin)*(fOut <= fmax))[0]
    
    # Compute the one-to-one fast ST
    SIGNAL = gft1d(sig,'gaussian');          
        
    # Compute the frequency and time associated with each sample of the one-to-one
    # fast ST
    F, T = dyadic2D(N)

    # Compute the extrapolation at the corner of the spectrogram using the nearest
    # neighboor approximation
    cornerT = [0, 0, N-1, N-1]
    cornerF = [indF[0], indF[-1], indF[0], indF[-1]]
    grid_real = griddata((T, F), SIGNAL.real, (cornerT, cornerF), method='nearest')
    grid_imag = griddata((T, F), SIGNAL.imag, (cornerT, cornerF), method='nearest')


    Tnew, Fnew = meshgrid(list(range(N)), indF)
    # Compute the linear interpolation, adding the values at the four corners to 
    # avoid extrapolation problems
    grid_real = griddata((concatenate((T, cornerT)), concatenate((F, cornerF))), 
                       concatenate((SIGNAL.real, grid_real.flatten())), (Tnew, Fnew), method='linear')    
    grid_imag = griddata((concatenate((T, cornerT)), concatenate((F, cornerF))), 
                       concatenate((SIGNAL.imag, grid_imag.flatten())), (Tnew, Fnew), method='linear')            
    return transpose(vectorize(complex)(grid_real, grid_imag)), fOut[indF] 



"""
    As computeFastST but compute only the real part, for a faster execution.
"""
def computeFastST_real(sig, fs, fmin=0, fmax=None):

    if not fmax:
        fmax = fs/2.0

    N  = len(sig)               # get the length of the signal

    f    = fftfreq(N)     
    fOut = f*fs    
    indF = where((fOut >= fmin)*(fOut <= fmax))[0]
    
    # Compute the one-to-one fast ST
    SIGNAL = gft1d(sig,'gaussian');          
        
    # Compute the frequency and time associated with each sample of the one-to-one
    # fast ST
    F, T = dyadic2D(N)

    # Compute the extrapolation at the corner of the spectrogram using the nearest
    # neighboor approximation
    cornerT = [0, 0, N-1, N-1]
    cornerF = [indF[0], indF[-1], indF[0], indF[-1]]
    grid_real = griddata((T, F), SIGNAL.real, (cornerT, cornerF), method='nearest')



    Tnew, Fnew = meshgrid(list(range(N)), indF)
    # Compute the linear interpolation, adding the values at the four corners to 
    # avoid extrapolation problems
    grid_real = griddata((concatenate((T, cornerT)), concatenate((F, cornerF))), 
                       concatenate((SIGNAL.real, grid_real.flatten())), (Tnew, Fnew), method='linear')    
            
    return transpose(grid_real), fOut[indF] 





