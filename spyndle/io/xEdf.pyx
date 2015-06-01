# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 13:53:47 2015

@author: oreichri
"""


"""
    int i, i2;
    for (i=0, i2=0; i<N; i+=3, i2++)
        ret(i2) = ((samplesU1(i+2) < 128 ? 0 : 255) << 24) | (samplesU1(i+2) << 16) | (samplesU1(i+1) << 8) | samplesU1(i);
    return_val = 1;
"""

import numpy as np
cimport numpy as np
cimport cython

#def encodeBDF(input):
@cython.boundscheck(False) # turn of bounds-checking for entire function
def encodeBDF(np.ndarray[unsigned char] input):
    cdef int i, i2, N
    N = len(input)
    cdef np.ndarray[np.int32_t] ret = np.zeros(N/3, dtype=np.dtype('int32'))
    for i, i2 in zip(range(0, N, 3), range(N/3)):
        ret[i2] = ((0 if input[i+2] < 128 else 255) << 24)  | \
                       (input[i+2] << 16) | (input[i+1] << 8)   | \
                       input[i]
    return ret





"""
@cython.boundscheck(False) # turn of bounds-checking for entire function
def encodeBDF(np.ndarray[unsigned char] input):
    cdef int i, i2, N
    N = len(input)
    
    cdef np.ndarray[np.int32_t] ret = np.zeros(N/3, dtype=np.dtype('uint32'))
    for i, i2 in zip(range(0, N, 3), range(N/3)):
        ret[i2] = ((0 if input[i+2] < 128 else 255) << 24)  | \
                       (input[i+2] << 16) | (input[i+1] << 8)   | \
                       input[i]
    return ret
"""
