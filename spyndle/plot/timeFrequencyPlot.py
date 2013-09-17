# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 16:07:19 2012

@author: REVESTECH
"""

import sys

from spyndle.miscellaneous import computeST

from scipy import where, array, transpose, fft
from scipy.fftpack import fftfreq

import numpy as np

from PyQt4.QtCore import *
from PyQt4.QtGui import *
from PyQt4.QtWebKit import *
from PyQt4.QtSvg import *


import matplotlib.gridspec as gridspec
#from matplotlib import colorbar
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.colors import Normalize
from matplotlib import cm




class TimeFreqPlot(QDialog):
    def __init__(self, parent=None): #, currentTime):
        QDialog.__init__(self, parent)

        self.setWindowTitle('Time-Frequency')
        self.mainLayout = QGridLayout()        
        self.setLayout(self.mainLayout)  
        self.resize(700,600)

        ####################################################################
        # Single graph
        gs = gridspec.GridSpec(2, 3,width_ratios=[2,8,1], height_ratios=[4,1])
        gs.update(left=0.075, right=0.95, wspace=0.025, hspace=0.025, top=0.95, bottom=0.075)

        self.figSingle = Figure(figsize=(600,600), dpi=72, facecolor=(1,1,1), edgecolor=(0,0,0))

        self.axContour  = self.figSingle.add_subplot(gs[1])
        self.axTime     = self.figSingle.add_subplot(gs[4])
        self.axColorbar = self.figSingle.add_subplot(gs[2])
        self.axFreq     = self.figSingle.add_subplot(gs[0])
      
        self.axContour.get_xaxis().set_visible(False)
        self.axContour.get_yaxis().set_visible(False)

        #self.axFreq.get_xaxis().set_visible(False)
        #self.axTime.get_yaxis().set_visible(False)
        
        self.axFreq.get_yaxis().set_label_text("Frequency (Hz)")
        self.axTime.get_xaxis().set_label_text("time (s)")
        self.axFreq.get_xaxis().set_label_text("ampltude (uV)")
        self.axTime.get_yaxis().set_label_text("ampltude (uV)")

        self.spectogramPlot = FigureCanvas(self.figSingle) #contourPlot() # Qwt.QwtPlot(self)            

        self.mainLayout.addWidget(self.spectogramPlot, 1, 1, 1, 1)   


    def plot(self, signal, fs, t, fmin=None, fmax=None, figname=None, 
                                             cm=cm.get_cmap("jet")):
        X, fX = computeST(signal, fs, fmin, fmax)    
        #return computeMST(signal, fs, m=0.0, k=1.0)
        
        
        tX = t 
        

        self.axContour.cla()
        self.axTime.cla()
        self.axFreq.cla()
        
        
        f = fftfreq(signal.size, 1.0/fs)[0:((signal.size+1)/2)]
        

        tmin = min(tX)
        tmax = max(tX)
        fmin = min(fX)
        fmax = max(fX)
        
        self.zmax = np.max(np.max(abs(X)))
        
        self.axContour.set_xlim(tmin, tmax)
        self.axContour.set_ylim(fmin, fmax)
        self.axTime.set_xlim(tmin, tmax)
        self.axFreq.set_ylim(fmin, fmax)
        
        self.axFreq.get_yaxis().set_label_text("Frequency (Hz)")
        self.axTime.get_xaxis().set_label_text("time (s)")
        self.axFreq.get_xaxis().set_label_text("ampltude\n(uV*s)")
        self.axTime.get_yaxis().set_label_text("ampltude (uV)")
        
        

        norm = Normalize(vmin=0, vmax=self.zmax)
        
        
        indT = where((array(tX) >= tmin)*(array(tX) <= tmax))[0]
        indF = where((array(fX) >= fmin)*(array(fX) <= fmax))[0]
        CS = self.axContour.imshow(abs(transpose(X[indT[0]:(indT[-1]+1), indF[0]:(indF[-1]+1)])), 
                                  interpolation = "bilinear",  aspect='auto', 
                                extent=[tX[indT[0]], tX[indT[-1]], fX[indF[0]], fX[indF[-1]]], 
                                    cmap=cm, norm=norm, origin="lower") # 
        
        
        
        self.axTime.plot(t, signal)
        self.axFreq.plot(abs(fft(signal)[0:((signal.size+1)/2)]), f)
        
        self.axTime.locator_params(axis ="y", nbins = 3)
        self.axTime.get_yaxis().set_ticks_position("right")
        self.axTime.get_yaxis().set_label_position("right")
        self.axFreq.locator_params(axis ="x", nbins = 2)    
        self.axFreq.invert_xaxis() 
         
         
        #cax, kw = self.axColorbar.make_axes_gridspec(CS)
         
        self.axContour.figure.colorbar(CS, use_gridspec=True, cax=self.axColorbar)    
        self.axContour.get_yaxis().set_label_text("ampltude (uV)")
        self.axColorbar.set_label("ampltitude (uV*s*Hz)")

        self.spectogramPlot.draw()

        if not figname is None:
            self.figSingle.savefig(figname, dpi=250, facecolor='w', edgecolor='w',
                    orientation='portrait', papertype=None, format=None,
                    transparent=False, bbox_inches=None, pad_inches=0.1)




def stPlot(t, signal, fmin=None, fmax=None, figname=None, cm=cm.get_cmap("jet")) :
    app = QApplication([])
    win = TimeFreqPlot()
    win.show()  
             
    fs = float(len(t)-1)/(t[-1]-t[0])
    win.plot(signal , fs, t, fmin=fmin, fmax=fmax, figname=figname, cm=cm)     
        
    sys.exit(app.exec_())

