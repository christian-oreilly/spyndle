# -*- coding: utf-8 -*-


"""
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
 Date  : July 27, 2012

"""
import  os
import numpy as np
import re
import math


from matplotlib.colors import hsv_to_rgb    
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
import matplotlib            

from scipy.interpolate import griddata
from scipy import array, reshape, concatenate, dot, arange, meshgrid, nan, unique, \
    arctan2, pi, cos, sin, zeros, where, isnan, maximum, minimum, ones

from spyndle.EEG import electrodesSVG


def getTransformedCoord(cx, cy, transform):
    coord = array(map(float, [cx, cy, 1.0]))
    
    mat = transforToMatrix(transform)
    
    return dot(mat, coord)[:-1]


def transforToMatrix(transform):
    mat = map(float, re.split(r'\(|\)|,', transform)[1:-1])
    mat = reshape(mat, (2, 3), order='F')
    mat = reshape(mat, (6))
    mat = concatenate((mat, [0, 0, 1]))
    mat = reshape(mat, (3, 3))
    return mat




# For an ellipse with center (cx, cy) and rays rx and ry, transformed following
# the affine transform "transform", compute the Cartesian coordinates arround the
# ellipse for angulare value r in (rads).
def transformEllipticalCoordinates(r, cx, cy, rx, ry, transform):
    xl = float(rx)*cos(r) + float(cx)
    yl = float(ry)*sin(r) + float(cy)

    mat = transforToMatrix(transform)
    
    res = zeros((len(xl), 2))
    for x, y, i in zip(xl, yl, range(len(xl))):
         #print mat.shape, array([x, y, 1.0]).shape
         res[i, :] = dot(mat, array([x, y, 1.0]))[:-1]

    return res




def getHeadCircle():
    cx, cy, rx, ry, transform = electrodesSVG.getHeadCircleProp()

    r = arange(0, 2.0*pi+0.1, 0.1)
    return transformEllipticalCoordinates(r, cx, cy, rx, ry, transform)



def transformCircle((cx, cy, rx, ry, transform)):

    r = array([0.0, 0.5])*pi
    res = transformEllipticalCoordinates(r, cx, cy, rx, ry, transform)
         
    cx = res[1, 0]
    rx = abs(res[0, 0] - cx)
    cy = res[0, 1]
    ry = abs(res[1, 1] - cy)
    res[:, 1]

    return cx, cy, rx, ry




def showColorMap(listeElectrodes, zVal):
    import pylab    

    zz = getZZ(listeElectrodes, zVal)

    circle = getHeadCircle()
    pylab.plot(circle[:, 0], circle[:, 1], 'k')    
            
    img=mpimg.imread(os.path.dirname(os.path.realpath(__file__))  + "\\eeg_electrodes_10-20_small.png")  
    
    #extent: [ None | scalars (left, right, bottom, top) ]
    pylab.imshow(zz, origin="upper", extent=[0.025, 0.97, -0.015, 0.935])    
    pylab.colorbar()
    pylab.imshow(img, extent=[0.0, 1.0, 0.0, 1.0])    
    pylab.show()    



def getZZ(listeElectrodes, zVal):

    cx, cy, transform = electrodesSVG.getElectrodeCoordinates(listeElectrodes)
    
    # Computing the coordinate of the electrodes.
    x_electrode = []
    y_electrode = []
    z_electrode = []
    for electrode in listeElectrodes:
        transformedCoord = getTransformedCoord(cx[electrode], cy[electrode], transform[electrode])
        x_electrode.append(transformedCoord[0])
        y_electrode.append(transformedCoord[1])
        z_electrode.append(zVal[electrode])
       

    # Computing inter/extrapolated amplitude of the zVal variable at coordinates
    # along the head circle. 
    circle = getHeadCircle()
    z_head = griddata((x_electrode, y_electrode), array(z_electrode), 
                     (circle[:, 0], circle[:, 1]), method='nearest')        
        
    x = concatenate((x_electrode, circle[:, 0]))
    y = concatenate((y_electrode, circle[:, 1]))
    z = concatenate((z_electrode, z_head))
        
    xx, yy = meshgrid(arange(min(x), max(x)), 
                      arange(min(y), max(y)))
    
    return griddata((x, y), z, (xx, yy), method='cubic')
            




def computeColorMapPval(listeElectrodes, zVal, pVal, minZZ=-0.02, maxZZ=0.02):

    import pylab        
    
    zz = getZZ(listeElectrodes, zVal)
    pp = getZZ(listeElectrodes, pVal)
    
    zz[where(isnan(zz))] = 0.0
    pp[where(isnan(pp))] = 0.0

    circle = getHeadCircle()
    pylab.plot(circle[:, 0], circle[:, 1], 'k')    
            
    
    s = (pp -0.2)/(0.8-0.2)
    h = (zz-minZZ)/(maxZZ-minZZ)

    Mat = ones((zz.shape[0], zz.shape[1], 3))    
    Mat[:, :, 0] = minimum(maximum(h, 0.0))   # H
    Mat[:, :, 1] = minimum(maximum(s, 0.0), 1.0)     # S
    # V = 1
    
    Mat = hsv_to_rgb(Mat)

    #extent: [ None | scalars (left, right, bottom, top) ]
    pylab.imshow(Mat, origin="upper", extent=[0.025, 0.97, -0.015, 0.935])       

    img=mpimg.imread(os.path.dirname(os.path.realpath(__file__))  + "\\eeg_electrodes_10-20_small.png")  
    pylab.imshow(img, extent=[0.0, 1.0, 0.0, 1.0])    






def computeColorMap(listeElectrodes, zVal, minZZ=0.00, maxZZ=0.02):
        
    import pylab    
        
    zz = getZZ(listeElectrodes, zVal)
    
    zz[where(isnan(zz))] = 0.0
    #s = zz*0.0+1.0
    #s[where(isnan(zz))] = 0.0
    #print s

    circle = getHeadCircle()
    pylab.plot(circle[:, 0], circle[:, 1], 'k')    
            
    

    M = zz.shape[0]
    N = zz.shape[1]
    #Mat = zeros((M, N, 3))    
    Mat = zeros((M, N))    

    deltaZZ = maxZZ-minZZ

    h = (zz-minZZ)/deltaZZ
 
    #v = 1.0
    for i1 in range(M):
        for i2 in range(N):
            hij = min(max(h[i1, i2], 0.0), 1.0)   
                
            #r, g, b = hsv2rgb(hij*240.0, s[i1, i2], v)
            #Mat[i1, i2, 0] = float(r)/255.0
            #Mat[i1, i2, 1] = float(g)/255.0
            #Mat[i1, i2, 2] = float(b)/255.0
            Mat[i1, i2] = hij
 
    pylab.imshow(Mat, origin="upper", extent=[0.025, 0.97, -0.015, 0.935])       


    img=mpimg.imread(os.path.dirname(os.path.realpath(__file__))  + "\\eeg_electrodes_10-20_small_whiteContour.png")  
    pylab.imshow(img, extent=[0.0, 1.0, 0.0, 1.0])    



def showColorMapPval(listeElectrodes, zVal, pVal):
    from matplotlib import pyplot
    computeColorMapPval(listeElectrodes, zVal, pVal)
    pyplot.show()    

def saveColorMapPval(listeElectrodes, zVal, pVal, filename, minZZ=-0.02, maxZZ=0.02):
    from matplotlib import pyplot
    computeColorMapPval(listeElectrodes, zVal, pVal, minZZ, maxZZ)
    pyplot.axis('off')
    pyplot.savefig(filename, transparent=False, bbox_inches='tight', pad_inches=0)




def saveColorMap(listeElectrodes, zVal, pVal, filename, minZZ=0.00, maxZZ=0.02):
    from matplotlib import pyplot
    computeColorMap(listeElectrodes, zVal, minZZ, maxZZ)
    pyplot.axis('off')
    pyplot.savefig(filename, transparent=False, bbox_inches='tight', pad_inches=0)







class HeadDrawing:
    def __init__(self, electrodeList = None):
        from matplotlib import patches
        
        if electrodeList is None:
            cx, cy, transform = electrodesSVG.getElectrodeCoordinates()
            self.__electrodeList = cx.keys()  
        else:
            self.__electrodeList = electrodeList  
            
        self.cx = {}
        self.cy = {}
        self.rx = {} 
        self.ry = {}
        self.electrodePatch = {}
        for elect in self.__electrodeList :
            self.cx[elect], self.cy[elect], self.rx[elect], self.ry[elect] = transformCircle(electrodesSVG.getElectCircleProp(elect))  
            self.electrodePatch[elect] = patches.Ellipse((self.cx[elect], self.cy[elect]), 2.0*self.rx[elect], 
                                                                              2.0*self.ry[elect], angle=0.0, fc="white")
    
        self.arrowPatch     = []
        
        self.ColorMap       = None
        self.s              = None
        self.colorbarDict   = None
    
    
    def setElectrodeList(self, electrodeList):
        self.__electrodeList = electrodeList
    
    
    def setElectrodeColor(self, elect, color):
        self.electrodePatch[elect].set_facecolor(color)
    
    def addArrow(self, source, dest, color="#000000", alpha=1.0):
        from matplotlib import patches
        angle = arctan2(self.cy[dest]-self.cy[source], self.cx[dest]-self.cx[source])

        x1 = self.rx[source]*cos(angle)  + self.cx[source] + 5*cos(angle+pi/2.0)
        y1 = self.rx[source]*sin(angle)  + self.cy[source] + 5*sin(angle+pi/2.0)
        x2 = -self.rx[dest]*cos(angle) + self.cx[dest] + 5*cos(angle+pi/2.0)
        y2 = -self.rx[dest]*sin(angle) + self.cy[dest] + 5*sin(angle+pi/2.0)
        
        self.arrowPatch.append(patches.FancyArrowPatch((x1, y1), (x2, y2), arrowstyle='-|>',mutation_scale=20, color=color, alpha=alpha))


    def addColorMap(self, listeElectrodes, zVal, minZZ=None, maxZZ=None, 
                                           sVal=None, smin=0.0, smax=1.0):
        
        valDict = {}      
        for elect, val in zip(listeElectrodes, zVal):
            valDict[elect] = val
        
        if minZZ is None:
            minZZ = min(zVal)
        if maxZZ is None:
            maxZZ = max(zVal)
        
                        
        h = (getZZ(listeElectrodes, valDict)-minZZ)/(maxZZ-minZZ)            
        
        s = None            
        if not sVal is None:
            sDict = {}      
            for elect, val in zip(listeElectrodes, sVal):
                sDict[elect] = val
            
            ss = getZZ(listeElectrodes, sDict)              

            deltaSS = smax-smin     
            if deltaSS == 0 :
                s = (ss >= smin).astype(np.float32)
            else:            
                s = (ss-smin)/deltaSS
                    
            self.s = ss


        self.ColorMap = ones((h.shape[0], h.shape[1], 3))                
        self.ColorMap[:, :, 0] = -(minimum(maximum(h, 0.0), 1.0)-1)*2.0/3.0      # H
        if not s is None:
            self.ColorMap[:, :, 1] = minimum(maximum(s, 0.0), 1.0)          # S
        #else S = 1
        # V = 1
        self.ColorMap[isnan(h), :] = [0.0, 0.0, 1.0]  #use white for nans       
        self.ColorMap = hsv_to_rgb(self.ColorMap)
     

        if s is None:    
            steps = arange(0, 1.1, 0.1)
            hsv = ones((1, 11, 3))
            hsv[:, :, 0] = -(steps-1)*2.0/3.0
            rgb = hsv_to_rgb(hsv)
            R = [((step, r, r)) for step, r in zip(steps, rgb[0, :, 0])]
            G = [((step, g, g)) for step, g in zip(steps, rgb[0, :, 1])]
            B = [((step, b, b)) for step, b in zip(steps, rgb[0, :, 2])]
                
            self.colorbarDict = {'red'  : tuple(R),
                                 'green': tuple(G),
                                 'blue' : tuple(B),
                                 'min':minZZ,
                                 'max':maxZZ}       





    def addColorMap_grey(self, listeElectrodes, zVal, minZZ=None, maxZZ=None, 
                                           sVal=None, smin=0.0, smax=1.0):
        
        valDict = {}      
        for elect, val in zip(listeElectrodes, zVal):
            valDict[elect] = val
        
        if minZZ is None:
            minZZ = min(zVal)
        if maxZZ is None:
            maxZZ = max(zVal)
        
                        
        h = (getZZ(listeElectrodes, valDict)-minZZ)/(maxZZ-minZZ)            
        
        s = None            
        if not sVal is None:
            sDict = {}      
            for elect, val in zip(listeElectrodes, sVal):
                sDict[elect] = val
            
            ss = getZZ(listeElectrodes, sDict)              

            deltaSS = smax-smin     
            if deltaSS == 0 :
                s = (ss >= smin).astype(np.float32)
            else:            
                s = (ss-smin)/deltaSS
                    
            self.s = ss
        
        
        self.ColorMap = ones((h.shape[0], h.shape[1], 3))         
        self.ColorMap[:, :, 0] = minimum(maximum(h, 0.0), 1.0)
        self.ColorMap[:, :, 1] = self.ColorMap[:, :, 0]
        self.ColorMap[:, :, 2] = self.ColorMap[:, :, 0]
        #else S = 1
        # V = 1
        self.ColorMap[isnan(h), :] = [1.0, 1.0, 1.0]  #use white for nans       
     
     

        steps = arange(0, 1.1, 0.1)
        hsv = ones((1, 11, 3))
        hsv[:, :, 0] = -(steps-1)*2.0/3.0
        rgb = hsv_to_rgb(hsv)
        R = [((step, r, r)) for step, r in zip(steps, rgb[0, :, 0])]
        G = [((step, g, g)) for step, g in zip(steps, rgb[0, :, 1])]
        B = [((step, b, b)) for step, b in zip(steps, rgb[0, :, 2])]
            
        self.colorbarDict = {'red'  : tuple(R),
                             'green': tuple(G),
                             'blue' : tuple(B),
                             'min':minZZ,
                             'max':maxZZ}       


    
    
    def setColorBarDict(self, cdict):
        self.colorbarDict = cdict
    
    
    def getColorBarDict(self):
        return self.colorbarDict


    def plot(self, save=False, filename="fig.svg", showColorBar=False):
            
        if save:
            matplotlib.use('Agg')
        import pylab    
    
    
        cx, cy, rx, ry = transformCircle(electrodesSVG.getHeadCircleProp())
        
        fig = pylab.figure(figsize=(10, 10.5))   
        
        if showColorBar and not self.colorbarDict is None:    
            gs = gridspec.GridSpec(1, 2, width_ratios=[10, 0.5]) 
            ax = pylab.subplot(gs[0])            
        else:
            ax = pylab.gca()        
        
        if not self.ColorMap is None:
            pylab.imshow(self.ColorMap, origin="upper", extent=[cx-rx, cx+rx, cy+ry, cy-ry])  
            if not self.s is None:
                pylab.contour(self.s, [-1, 0.95], colors='k',  origin="upper", linestyles="dashed", 
                              linewidths=3.0, extent=[cx-rx, cx+rx, cy+ry, cy-ry])
             
    
        img=mpimg.imread(os.path.dirname(os.path.realpath(__file__))  + "\\eeg_electrodes_10-20_small_clear.png")  
        pylab.imshow(img, extent=[cx-rx-59, cx+rx+59, cy+ry+17, cy-ry-72])  
    
    
        for elect in self.__electrodeList : 
            ax.add_patch(self.electrodePatch[elect] )
    

            ax.text(self.cx[elect], self.cy[elect], elect,
                        horizontalalignment='center',
                        verticalalignment='center',
                        fontsize=14, color='black')

    
    
    
        for arrow in self.arrowPatch:
            ax.add_patch(arrow)
                
        pylab.axis('off')
        
        if showColorBar and not self.colorbarDict is None:
            ax = pylab.subplot(gs[1])
            my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',self.colorbarDict, 256)
            data = np.outer(np.arange(101),np.ones(5)) 
            pylab.imshow(data, cmap=my_cmap, origin="lower")
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().tick_right()
            ticks = np.arange(0, 101, 25)
            ax.set_yticks(ticks)
            ax.set_yticklabels(np.round(ticks/100.0*(self.colorbarDict["max"]-
                        self.colorbarDict["min"])+self.colorbarDict["min"], 4))
            gs.tight_layout(fig)



        if save:
            pylab.savefig(filename,dpi=100, transparent=False, bbox_inches='tight', pad_inches=0)
        else:
            pylab.show()









def plotArrowsBidirect(meanVal, pairs, minDelay=None, maxDelay=None, 
                                      filename=None, addColorbar=True):

    listVal = []
    for ref, test in pairs :
        try:        
            listVal.append(meanVal[ref][test])
        except KeyError:
            listVal.append(nan)
    
    listVal = array(listVal)        
    
    plotArrowDiagram(listVal, pairs, 
                     minDelay, maxDelay, filename, addColorbar)    



def plotArrowsDifferential(meanVal, pairs, minDelay=None, maxDelay=None, 
                          filename=None, addColorbar=True):
    
    listVal = []
    for ref, test in pairs :
        try:        
            listVal.extend([meanVal[ref][test] - meanVal[test][ref]])
        except KeyError:
            listVal.append(nan)
    
    
    
    listVal = array(listVal)
    
    IND       = where(listVal > 0.0)[0]
    listVal   = listVal[IND]
    pairs     = [pairs[i] for i in IND]
    
    plotArrowDiagram(listVal, pairs, 
                     minDelay, maxDelay, filename, addColorbar)
    
    
    
def plotArrowDiagram(listVal, pairs, minDelay=None, maxDelay=None, 
                          filename=None, addColorbar=True):    
    
    if isinstance(listVal, list):
        listVal = array(listVal)
    
    if minDelay is None :
        minDelay = min(listVal[np.logical_not(np.isnan(listVal))])
    if maxDelay is None :    
        maxDelay = max(listVal[np.logical_not(np.isnan(listVal))])
    
    normDelay = (listVal - minDelay)/(maxDelay - minDelay)
    normDelay = np.maximum(normDelay, 0.0)
    normDelay = np.minimum(normDelay, 1.0)

    drawing = HeadDrawing()
    refs  = unique([pair[0] for pair in pairs])
    tests = unique([pair[1] for pair in pairs])
    drawing.setElectrodeList(unique(concatenate((refs, tests))))

    for pair, delay in zip(pairs, normDelay) :
        if not isnan(delay):
            r  = min(1.0, delay*2.0)
            gb = max(0.0, (delay*2.0 - 1.0)*200/255)                     
            drawing.addArrow(pair[0], pair[1], (r, gb, gb), 1.0)        
    
    if addColorbar:
        cdict = {'red': ((0.0, 0.0, 0.0),
                         (0.5, 1.0, 1.0),
                         (1.0, 1.0, 1.0)),
                 'green': ((0.0, 0.0, 0.0),
                           (0.5, 0.0, 0.0),
                           (1.0, 200./255., 200./255.)),
                 'blue': ((0.0, 0.0, 0.0),
                          (0.5, 0.0, 0.0),
                          (1.0, 200./255., 200./255.)),
                 'min':minDelay,
                 'max':maxDelay}   
        drawing.setColorBarDict(cdict)
    
    if filename is None:
        drawing.plot(False, showColorBar=addColorbar)
    else:
        drawing.plot(True, filename, showColorBar=addColorbar)



def plotColorMap(electrodes, var, filename=None, addColorbar=True, minVal=None, maxVal=None):
    if not filename is None:
        matplotlib.use('Agg')
        
    drawing = HeadDrawing(electrodeList = np.unique(electrodes))
    
    drawing.addColorMap(electrodes, var, minZZ=minVal, maxZZ=maxVal)

    if filename is None:
        drawing.plot(False, showColorBar=addColorbar)
    else:
        drawing.plot(True, filename, showColorBar=addColorbar)

