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

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


 Author: Christian O'Reilly (christian.oreilly@umontreal.ca)
 Date  : June 21, 2012

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



from scipy.signal import firwin
from scipy.signal._arraytools import odd_ext, even_ext, const_ext      

from scipy import signal, zeros, angle
import numpy, random

from scipy.weave import converters
import scipy.weave as weave
import numpy as np


channelType = {"EEG": 1, "EMG": 3, "EOG": 7, "ECG": 8, "MIC": 10, "RSP": 23}


class Filter:
    def __init__(self, samplingRate):
        self.samplingRate = samplingRate
        

    def create(self, low_crit_freq, high_crit_freq, order, btype="lowpass", ftype="butter",
               rp = None, rs = None, useFiltFilt=True, channelType = channelType["EEG"]):        


        self.low_crit_freq  = low_crit_freq        
        self.high_crit_freq = high_crit_freq   
        self.order          = order 
        self.btype          = btype 
        self.ftype          = ftype 
        self.rp             = rp 
        self.rs             = rs    
        self.useFiltFilt    = useFiltFilt
        
        self.channelType = channelType    # EEG, OEG, EMG...    

        self.update()
            

    # À appeller après la modification des propriétés de l'objet pour que le filtre y soit conforme.
    def update(self):
        Nyq = self.samplingRate/2.0        
        
        
        ############################# FIR filters #############################
        if self.ftype == "FIR":
            if  self.btype=='lowpass':
                pass_zero = True
                self.low_crit_freq  = 0.0
            elif self.btype=='highpass':
                pass_zero = False
                self.high_crit_freq = Nyq
            elif self.btype=='bandpass':
                pass_zero = False
                
            
            self.b = firwin(self.order, [self.low_crit_freq/Nyq, self.high_crit_freq/Nyq], window="hanning", pass_zero=pass_zero)
            self.a = [1.0]
                
       ########################### other filters ##############################     
        else:
                
            if self.btype=='bandpass' or  self.btype=='lowpass' :
                Wnh = self.high_crit_freq/Nyq
                self.bl, self.al = signal.iirfilter(self.order, Wnh, rp=self.rp, rs=self.rs, btype="lowpass", ftype=self.ftype)
                
            if self.btype=='bandpass' or  self.btype=='highpass' :
                Wnl = self.low_crit_freq/Nyq
                self.bh, self.ah = signal.iirfilter(self.order, Wnl, rp=self.rp, rs=self.rs, btype="highpass", ftype=self.ftype)        
            
            

    def applyFilter(self, sig2filt):
    
        if self.ftype == "FIR":
            
            try:
                if self.useFiltFilt:
                    sig2filt = filtfilt_FFTWEAVE(self.b, sig2filt, padtype='even')
                else:
                    sig2filt = filt_FFTWEAVE(self.b, sig2filt, padtype='even')       
            except ValueError:
                print "order:", self.order
                raise

        else:
            Nd = min(3*self.samplingRate, len(sig2filt))   
            
            if self.useFiltFilt:
                if self.btype=='bandpass': 
                    sig2filt = signal.filtfilt(self.bl, self.al, signal.filtfilt(self.bh, self.ah, sig2filt, padtype="even", padlen=Nd), padtype="even", padlen=Nd)
                elif self.btype=='lowpass' :
                    sig2filt = signal.filtfilt(self.bl, self.al, sig2filt, padtype="even", padlen=Nd)
                elif self.btype=='highpass' : 
                    sig2filt = signal.filtfilt(self.bh, self.ah, sig2filt, padtype="even", padlen=Nd)
            else:
                if self.btype=='bandpass': 
                    sig2filt = signal.lfilter(self.bl, self.al, signal.lfilter(self.bh, self.ah, sig2filt, padtype="even", padlen=Nd), padtype="even", padlen=Nd)
                elif self.btype=='lowpass' :
                    sig2filt = signal.lfilter(self.bl, self.al, sig2filt, padtype="even", padlen=Nd)
                elif self.btype=='highpass' : 
                    sig2filt = signal.lfilter(self.bh, self.ah, sig2filt, padtype="even", padlen=Nd)            
            
        return sig2filt
    


    def computeFilterResponse(self):        
        N = int(100*self.samplingRate)

        noise = numpy.zeros(N)
        freq = numpy.fft.fftfreq(N, 1.0/self.samplingRate)            

        for i in range(N) :
            noise[i] = random.random()-0.5
            
        fftN = numpy.fft.fft(noise)[range(N/2)]          
        fftS = numpy.fft.fft(self.applyFilter(noise))[range(N/2)]    
        fftNA  = abs(fftN)            
        fftSA  = abs(fftS)    
        fftP  = angle(fftS) - angle(fftN)    

        Niter = 100
        for i in range(1, Niter):
            for i in range(N) :
                noise[i] = random.random()-0.5

            fftN = numpy.fft.fft(noise)[range(N/2)]          
            fftS = numpy.fft.fft(self.applyFilter(noise))[range(N/2)]    
            fftNA  += abs(fftN)            
            fftSA  += abs(fftS)    
            fftP  += (angle(fftS) - angle(fftN))                

        return freq[range(N/2)], fftSA/Niter, fftNA/Niter, fftP/Niter        

        
        
        
        
        
        

# https://gist.github.com/astrofrog/837209

# b must be odd, lena > lenb
def convolve_weave(a,b):

    lena = len(a)
    lenb = len(b)
    convol = zeros(a.shape)               

    code = """
                int pad = (lenb-1)/2;
                int i, j;
             
                for (i=pad; i<pad+lena; i++)
                {
                    int kmin, kmax, k;
                    j = i-pad;

                    convol(j)  = 0;
                
                    kmin = (i >= lenb - 1) ? i - (lenb - 1) : 0;
                    kmax = (i <  lena - 1) ? i : lena - 1;
                
                    for (k = kmin; k <= kmax; k++)
                    {
                      convol(j)  += a(k)*b(i - k);
                    }                       
                }
                return_val = 1;
           """

    weave.inline(code, [ 'a', 'b', 'lena', 'lenb', 'convol'],
    type_converters=converters.blitz, compiler = 'gcc')

    return convol








# b must be odd, lena > lenb
def convolve_weaveFiltFilt(a,b):

    lena = len(a)
    lenb = len(b)
    convol = zeros(a.shape)

    code = """
                // Forward convolution
                int pad = (lenb-1)/2;
                int i, j;
             
                for (i=pad; i<pad+lena; i++)
                {
                    int kmin, kmax, k;
                    // Reverse indexing for the next pass
                    j = lena-1-i+pad;

                    convol(j)  = 0;
                
                    kmin = (i >= lenb - 1) ? i - (lenb - 1) : 0;
                    kmax = (i <  lena - 1) ? i : lena - 1;
                
                    for (k = kmin; k <= kmax; k++)
                    {
                      convol(j)  += a(k)*b(i - k);
                    }                       
                }
                
                
                // Backward convolution (the signal in convol has been
                // reversed using reversed indexes)            
                for (i=pad; i<pad+lena; i++)
                {
                    int kmin, kmax, k;
                    // Reverse indexing for reordering the output vector
                    j = lena-1-i+pad;

                    a(j)  = 0;
                
                    kmin = (i >= lenb - 1) ? i - (lenb - 1) : 0;
                    kmax = (i <  lena - 1) ? i : lena - 1;
                
                    for (k = kmin; k <= kmax; k++)
                    {
                      a(j)  += convol(k)*b(i - k);
                    }                       
                }                
                
                return_val = 1;
         """  

    weave.inline(code, [ 'a', 'b', 'lena', 'lenb', 'convol'],
    type_converters=converters.blitz, compiler = 'gcc')


            

## COPIED FROM scipy filtfilt but changing the call to lfilter to
## fftconvolve because it is much faster.

def filt_FFTWEAVE(b, x, padtype='odd', padlen=None):
    """A forward-backward filter.

This function applies a linear filter twice, once forward
and once backwards. The combined filter has linear phase.

Before applying the filter, the function can pad the data along the
given axis in one of three ways: odd, even or constant. The odd
and even extensions have the corresponding symmetry about the end point
of the data. The constant extension extends the data with the values
at end points. On both the forward and backwards passes, the
initial condition of the filter is found by using lfilter_zi and
scaling it by the end point of the extended data.

Parameters
----------
b : array_like, 1-D
The numerator coefficient vector of the filter.
a : array_like, 1-D
The denominator coefficient vector of the filter. If a[0]
is not 1, then both a and b are normalized by a[0].
x : array_like
The array of data to be filtered.
padtype : str or None, optional
Must be 'odd', 'even', 'constant', or None. This determines the
type of extension to use for the padded signal to which the filter
is applied. If `padtype` is None, no padding is used. The default
is 'odd'.
padlen : int or None, optional
The number of elements by which to extend `x` at both ends of
`axis` before applying the filter. This value must be less than
`x.shape[axis]-1`. `padlen=0` implies no padding.
The default value is 3*max(len(a),len(b)).

Returns
-------
y : ndarray
The filtered output, an array of type numpy.float64 with the same
shape as `x`.

See Also
--------
lfilter_zi
lfilter

Examples
--------
First we create a one second signal that is the sum of two pure sine
waves, with frequencies 5 Hz and 250 Hz, sampled at 2000 Hz.

>>> t = np.linspace(0, 1.0, 2001)
>>> xlow = np.sin(2 * np.pi * 5 * t)
>>> xhigh = np.sin(2 * np.pi * 250 * t)
>>> x = xlow + xhigh

Now create a lowpass Butterworth filter with a cutoff of 0.125 times
the Nyquist rate, or 125 Hz, and apply it to x with filtfilt. The
result should be approximately xlow, with no phase shift.

>>> from scipy.signal import butter
>>> b, a = butter(8, 0.125)
>>> y = filtfilt(b, a, x, padlen=150)
>>> np.abs(y - xlow).max()
9.1086182074789912e-06

We get a fairly clean result for this artificial example because
the odd extension is exact, and with the moderately long padding,
the filter's transients have dissipated by the time the actual data
is reached. In general, transient effects at the edges are
unavoidable.
"""

    if padtype not in ['even', 'odd', 'constant', None]:
        raise ValueError(("Unknown value '%s' given to padtype. padtype must "
                         "be 'even', 'odd', 'constant', or None.") %
                            padtype)

    b = np.asarray(b)
    x = np.asarray(x)

    ntaps = len(b)

    if padtype is None:
        padlen = 0

    if padlen is None:
        # Original padding; preserved for backwards compatibility.
        edge = ntaps * 3
    else:
        edge = padlen

    # x's 'axis' dimension must be bigger than edge.
    #if x.shape[axis] <= edge:
    if len(x) <= edge:
        raise ValueError("The length of the input vector x must be larger than "
                         "padlen, which is %d." % edge)

    if padtype is not None and edge > 0:
        # Make an extension of length `edge` at each
        # end of the input array.
        if padtype == 'even':
            ext = even_ext(x, edge)#, axis=axis)
        elif padtype == 'odd':
            ext = odd_ext(x, edge)#, axis=axis)
        else:
            ext = const_ext(x, edge)#, axis=axis)
    else:
        ext = x

    # Get the steady state of the filter's step response.
    #zi = lfilter_zi(b, a)

    # Reshape zi and create x0 so that zi*x0 broadcasts
    # to the correct value for the 'zi' keyword argument
    # to lfilter.
    #zi_shape = [1] * x.ndim
    #zi_shape[axis] = zi.size
    #zi = np.reshape(zi, zi_shape)
    #x0 = axis_slice(ext, stop=1, axis=axis)

    # Forward filter.
    ext = convolve_weave(ext, b)


    if edge > 0:
        # Slice the actual signal from the extended signal. Reverse and return y.
        return ext[edge:-edge]
    else:
        # Reverse and return y.
        return ext
            
            
            

def filtfilt_FFTWEAVE(b, x, #axis=-1, 
                         padtype='odd', padlen=None):
    """A forward-backward filter.

This function applies a linear filter twice, once forward
and once backwards. The combined filter has linear phase.

Before applying the filter, the function can pad the data along the
given axis in one of three ways: odd, even or constant. The odd
and even extensions have the corresponding symmetry about the end point
of the data. The constant extension extends the data with the values
at end points. On both the forward and backwards passes, the
initial condition of the filter is found by using lfilter_zi and
scaling it by the end point of the extended data.

Parameters
----------
b : array_like, 1-D
The numerator coefficient vector of the filter.
a : array_like, 1-D
The denominator coefficient vector of the filter. If a[0]
is not 1, then both a and b are normalized by a[0].
x : array_like
The array of data to be filtered.
axis : int, optional
The axis of `x` to which the filter is applied.
Default is -1.
padtype : str or None, optional
Must be 'odd', 'even', 'constant', or None. This determines the
type of extension to use for the padded signal to which the filter
is applied. If `padtype` is None, no padding is used. The default
is 'odd'.
padlen : int or None, optional
The number of elements by which to extend `x` at both ends of
`axis` before applying the filter. This value must be less than
`x.shape[axis]-1`. `padlen=0` implies no padding.
The default value is 3*max(len(a),len(b)).

Returns
-------
y : ndarray
The filtered output, an array of type numpy.float64 with the same
shape as `x`.

See Also
--------
lfilter_zi
lfilter

Examples
--------
First we create a one second signal that is the sum of two pure sine
waves, with frequencies 5 Hz and 250 Hz, sampled at 2000 Hz.

>>> t = np.linspace(0, 1.0, 2001)
>>> xlow = np.sin(2 * np.pi * 5 * t)
>>> xhigh = np.sin(2 * np.pi * 250 * t)
>>> x = xlow + xhigh

Now create a lowpass Butterworth filter with a cutoff of 0.125 times
the Nyquist rate, or 125 Hz, and apply it to x with filtfilt. The
result should be approximately xlow, with no phase shift.

>>> from scipy.signal import butter
>>> b, a = butter(8, 0.125)
>>> y = filtfilt(b, a, x, padlen=150)
>>> np.abs(y - xlow).max()
9.1086182074789912e-06

We get a fairly clean result for this artificial example because
the odd extension is exact, and with the moderately long padding,
the filter's transients have dissipated by the time the actual data
is reached. In general, transient effects at the edges are
unavoidable.
"""

    if padtype not in ['even', 'odd', 'constant', None]:
        raise ValueError(("Unknown value '%s' given to padtype. padtype must "
                         "be 'even', 'odd', 'constant', or None.") %
                            padtype)

    b = np.asarray(b)
    x = np.asarray(x)

    ntaps = len(b)

    if padtype is None:
        padlen = 0

    if padlen is None:
        # Original padding; preserved for backwards compatibility.
        edge = ntaps * 3
    else:
        edge = padlen

    # x's 'axis' dimension must be bigger than edge.
    #if x.shape[axis] <= edge:
    if len(x) <= edge:
        raise ValueError("The length of the input vector x must be larger than "
                         "padlen, which is %d." % edge)

    if padtype is not None and edge > 0:
        # Make an extension of length `edge` at each
        # end of the input array.
        if padtype == 'even':
            ext = even_ext(x, edge)#, axis=axis)
        elif padtype == 'odd':
            ext = odd_ext(x, edge)#, axis=axis)
        else:
            ext = const_ext(x, edge)#, axis=axis)
    else:
        ext = x

    convolve_weaveFiltFilt(ext, b)

    if edge > 0:
        # Slice the actual signal from the extended signal. Reverse and return y.
        return ext[edge:-edge] 
    else:
        # Reverse and return y.
        return ext
                                    