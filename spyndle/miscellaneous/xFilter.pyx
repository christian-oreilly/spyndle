# -*- coding: utf-8 -*-


"""
    Copyright (C) 2012-2015  Christian O'Reilly

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
 Date  : February 17, 2015

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
        

# https://gist.github.com/astrofrog/837209

# b must be odd, lena > lenb

import numpy as np
cimport numpy as np
cimport cython

@cython.boundscheck(False) # turn of bounds-checking for entire function
def convolve_cython(np.ndarray[double] a, np.ndarray[double] b):

    lena = len(a)
    lenb = len(b)
    cdef double[:] convol = np.zeros(lena)    

    cdef int pad = (lenb-1)/2
    cdef int i, j
    cdef int kmin, kmax, k

 
    for i in xrange(pad, pad+lena):
        j = i-pad

        convol[j]  = 0
    
        kmin = i - (lenb - 1) if (i >= lenb - 1) else 0
        kmax = i if (i <  lena - 1) else lena - 1
    
        for k in xrange(kmin, kmax+1):
          convol[j]  += a[k]*b[i - k]

    return convol



# b must be odd, lena > lenb

@cython.boundscheck(False) # turn of bounds-checking for entire function
def convolve_cythonFiltFilt(np.ndarray[double] a, np.ndarray[double] b):
    
    lena = len(a)
    lenb = len(b)
    cdef double[:] convol = np.zeros(lena)

    # Forward convolution
    cdef int pad = (lenb-1)/2
    cdef int i, j
    cdef int kmin, kmax, k
 
    for i in xrange(pad, pad+lena):
        # Reverse indexing for the next pass
        j = lena-1-i+pad

        convol[j]  = 0
    
        kmin = i - (lenb - 1) if (i >= lenb - 1) else 0
        kmax = i if (i <  lena - 1) else lena - 1
    
        for k in xrange(kmin, kmax+1):
          convol[j]  += a[k]*b[i - k]
    
    
    # Backward convolution (the signal in convol has been
    # reversed using reversed indexes)            
    for i in xrange(pad, pad+lena):

        # Reverse indexing for reordering the output vector
        j = lena-1-i+pad

        a[j]  = 0
    
        kmin = i - (lenb - 1) if (i >= lenb - 1) else 0
        kmax = i if (i <  lena - 1) else lena - 1
    
        for k in xrange(kmin, kmax+1): 
          a[j]  += convol[k]*b[i - k]

    return a
            
