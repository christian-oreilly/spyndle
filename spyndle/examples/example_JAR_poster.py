# -*- coding: utf-8 -*-
"""
    Code used to produce figures included in a poster presention. The poster 
    can be found at https://www.researchgate.net/profile/Christian_OReilly/
    publication/256903596_Spyndle_Un_paquetage_Python__code_ouvert_pour_
    l%27analyse_des_transitoires_EEG/file/5046352408141c648b.pdf

    Copyright (C) 2013  Christian O'Reilly

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
 Date  : September 16, 2013


"""

import os
import numpy as np
from scipy.stats.mstats import mquantiles

from spyndle import DevuystDB
from spyndle.io import EDFReader
from spyndle.detector import SpindleDetectorRMS
from spyndle.plot import stPlot
from spyndle import Filter


path = os.path.join(os.path.dirname(__file__), "imagesJARexample")

if not os.path.exists(path):
    os.makedirs(path)

# Loading example signal from Devuyst database
filename = os.path.join(DevuystDB._devuystDBPath_, "excerpt1.BDF")
if not os.path.isfile(filename):
    from spyndle.DevuystDB import convertDevuystToBDF
    convertDevuystToBDF()
    
reader  = EDFReader(filename)
channel = reader.getChannelLabels()[0]

recSignals = reader.read([channel], 0, 35.0)

signal = recSignals[channel].signal
fs     = recSignals[channel].samplingRate 
t      = np.arange(0, len(signal))/fs



###############################################################################
# Making and saving an example of ST plot
###############################################################################
stPlot(t[:300], signal[:300], fmin=0, fmax=20, 
       figname=os.path.join(path, "stExample.png"))


###############################################################################
# Making and saving an example of filtering plots
###############################################################################

# Defining EEG filters
lowPassFilter = Filter(fs)
lowPassFilter.create(low_crit_freq=None, 
                      high_crit_freq=10.0, order=1001, 
                      btype="lowpass", ftype="FIR", useFiltFilt=True)          
                      
highPassFilter = Filter(fs)
highPassFilter.create(low_crit_freq=8.0, 
                      high_crit_freq=None, order=1001, 
                      btype="highpass", ftype="FIR", useFiltFilt=True)          


bandPassFilter = Filter(fs)
bandPassFilter.create(low_crit_freq=8.0, 
                      high_crit_freq=10.0, order=1001, 
                      btype="bandpass", ftype="FIR", useFiltFilt=True)          



lowSignal  = lowPassFilter.applyFilter(signal)     
highSignal = highPassFilter.applyFilter(signal)     
bandSignal = bandPassFilter.applyFilter(signal)     

from matplotlib import pyplot as plt

f, axarr = plt.subplots(3, sharex=True)

axarr[0].plot(t[:100], signal[:100])
axarr[0].plot(t[:100], lowSignal[:100])
axarr[0].set_title(u"passe-bas à 10 Hz")
axarr[0].set_ylabel(u"amplitude (uV)")
axarr[0].legend((u"original", u"filtré"), loc="lower center", ncol =2, frameon=False)

axarr[1].plot(t[:100], signal[:100])
axarr[1].plot(t[:100], highSignal[:100])
axarr[1].set_title(u"passe-haut à 8 Hz")
axarr[1].set_ylabel(u"amplitude (uV)")
axarr[1].legend((u"original", u"filtré"), loc="lower center", ncol =2, frameon=False) 
 
axarr[2].plot(t[:100], signal[:100])
axarr[2].plot(t[:100], bandSignal[:100])
axarr[2].set_title(u"passe-bande 8-10 Hz")
axarr[2].set_xlabel(u"temps (s)")
axarr[2].set_ylabel(u"amplitude (uV)")
axarr[2].legend((u"original", u"filtré"), loc="lower center", ncol =2, frameon=False)

plt.savefig(os.path.join(path, "filterExample.png"))
print "Saving " + os.path.join(path, "filterExample.png")
plt.show()



###############################################################################
# Making and saving a figure showing examples of spindles.
###############################################################################

spindleEvents = filter(lambda e: e.name == "spindleV1" and e.channel == channel, reader.events)    

startPad    = 4.0                               # Start zero padding 

bandPassFilter = Filter(fs)
bandPassFilter.create(low_crit_freq=2.0, 
                      high_crit_freq=30.0, order=251, 
                      btype="bandpass", ftype="FIR", useFiltFilt=True)              


f, axarr = plt.subplots(2, 2, sharex='col', sharey='row')
f.set_size_inches(12, 3)
f.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)
for event, ax, i in zip(np.array(spindleEvents)[[0, 9, 18, 49]], axarr.flat, range(4)) :
    data = reader.read([channel], event.startTime -startPad, 8.0)  
    t    = np.arange(len(data[channel].signal))/fs    # sampling time

    ax.plot(t, bandPassFilter.applyFilter(data[channel].signal))     
    ax.set_ylim([-60, 60])
    ax.set_frame_on(False)
    if i == 0 or i == 2:
        ax.get_yaxis().tick_left()
        ax.set_ylabel("uV")
        ax.set_yticks([-60, 0, 60])
    else:
        ax.get_yaxis().set_visible(False)
    if i >= 2:
        ax.get_xaxis().tick_bottom()
        ax.set_xlabel("secondes")        
        ax.set_xticks([0, 2, 4, 6, 8])
    else:
        ax.get_xaxis().set_visible(False)

    ax.plot([startPad, startPad], [-70, 70], color="red", linewidth=2, linestyle="--")     
    ax.plot([startPad+event.duration(), startPad+event.duration()], [-70, 70],
                                      color="red", linewidth=2, linestyle="--")     

plt.savefig(os.path.join(path, "spindleExample.png"), dpi=250)
print "Saving " + os.path.join(path, "spindleExample.png")
plt.show()







###############################################################################
# Making and saving a figure showing the spindle detection principle.
###############################################################################

bandPassFilter = Filter(fs)
bandPassFilter.create(low_crit_freq=11.0, 
                      high_crit_freq=16.0, order=251, 
                      btype="bandpass", ftype="FIR", useFiltFilt=True)              

detector = SpindleDetectorRMS()

event = spindleEvents[49]
data = reader.read([channel], event.startTime -2.0, 20.0)  

raw = data[channel].signal
signal = bandPassFilter.applyFilter(data[channel].signal)     
                   
windowNbSample = int(round(detector.averagingWindowSize*fs))
if np.mod(windowNbSample, 2) == 0 : # We need an odd number.
    windowNbSample += 1

rectSig = detector.averaging(np.abs(signal), windowNbSample)

treshold = mquantiles(rectSig, detector.threshold)
t = np.arange(len(rectSig[0:500]))/fs

f, axarr = plt.subplots(1, 2)
f.set_size_inches(12, 2)
f.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.2)

axarr[0].plot(t, raw[0:500])
axarr[0].set_ylim([-60, 60])
axarr[0].set_xlim([0, 5])
axarr[0].set_xlabel("secondes")   
axarr[0].set_ylabel("EEG (uV)")       

axarr[1].plot(t, rectSig[0:500])
axarr[1].plot([-1, 9], [treshold, treshold], color="green", linewidth=2, linestyle=":")
axarr[1].set_xlim([0, 5])
axarr[1].set_ylim([0, 16])
inds = np.where(rectSig[0:500] >= treshold)[0]
axarr[1].plot(t[[inds[0], inds[0]]], [-1, 17], color="red", linewidth=2, linestyle="--")
axarr[1].plot(t[[inds[-1], inds[-1]]], [-1, 17], color="red", linewidth=2, linestyle="--")
axarr[1].set_xlabel("secondes")   
axarr[1].set_ylabel("Amplitude f(t)")     

plt.savefig(os.path.join(path, "detectionExample.png"), dpi=250)
print "Saving " + os.path.join(path, "detectionExample.png")
plt.show()

