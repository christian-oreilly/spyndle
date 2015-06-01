# -*- coding: utf-8 -*-


"""
    Various utilitary functions used by the Spyndle package.

    Copyright (C) 2008-2013  Christian O'Reilly

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
 Date  : August 28, 2008

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



from scipy import zeros, arange
import sys, os
import numpy as np

def setUnbufferedPrint():
    sys.stdout=Unbuffered(sys.stdout)

class Unbuffered:
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)



def terminate_process(pid):
    if sys.platform == 'win32':
        import ctypes
        PROCESS_TERMINATE = 1
        handle = ctypes.windll.kernel32.OpenProcess(PROCESS_TERMINATE, False, pid)
        ctypes.windll.kernel32.TerminateProcess(handle, -1)
        ctypes.windll.kernel32.CloseHandle(handle)
    else:
        os.kill(pid, signal.SIGTERM)




def findMaxima(y) :
    #  Find location of local maxima

    # Identify whether signal is rising or falling
    upordown = np.sign(np.diff(y))
    if len(upordown) == 0 :
        return []
    
    # Find points where signal is rising before, falling after
    maxflags = np.concatenate(([upordown[0]<0], np.diff(upordown)<0, [upordown[-1]>0]))
    
    return np.where(maxflags)[0]
    
    
def findMinima(y) :
    #  Find location of local minima

    return findMaxima(-1.0*np.array(y))    
    
    
    
def findZeroCrossings(y):
    return np.where(y[:-1] * y[1:] < 0)[0]


def diff2(x, y) :
    # Second order numerical derivative. 
    #
    # First point is approximated using a forward derivative, the last point
    # using a backward derivative and the other points using a centered derivative.
    #
    #
    # Dérivé numérique d'ordre 2
    #
    # Le permier point est approximé par une dérivée avant, le dernier points
    # par une dérivée arrière et les autres points par une dérivée centrée. 
    #
    # Voir Analyse numérique pour ingénieurs, 2eme ed. de André Fortin, p. 318
    #
    # Christian O'Reilly
    # 26 août 2008
    #
    # Modifications
    #    15 septembre 2008 : Traduction de Matlab à Python.

    if len(x) != len(y) :
         raise 'Les vecteurs x et y doivent être de même taille. (' + str(len(x)) + ', ' + str(len(y)) + ')'
    
    if len(x) < 3 :
        raise 'Les vecteurs x et y doivent avoir au moins 3 éléments.'
    
    dy = zeros(len(x))
    
    # Dérivée numérique avant
    dy[0]   = ( -y[2] + 4*y[1] - 3*y[0] )/( x[2] - x[0] )
    
    # Dérivée numérique arrière
    dy[-1] = ( 3*y[-1] - 4*y[-2] + y[-3] )/( x[-1] - x[-3] )
    
    # Dérivée numérique centrée
    i = arange(1, len(x)-1)
    dy[i] = ( y[i+1] - y[i-1] )/( x[i+1] - x[i-1] )
    
    return dy
     