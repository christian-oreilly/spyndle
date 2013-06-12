# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 16:21:15 2012

@author: REVESTECH
"""
import sys, os
sys.path.append('\\\\HERMES\\DNL-net\\Members\\Christian\\Python\\src\\')


from scipy.interpolate import interp2d, bisplrep, bisplev, NearestNDInterpolator, griddata
from scipy import array, reshape, concatenate, dot, arange, meshgrid, transpose
from scipy import arctan2, pi, cos, sin, zeros, where, isnan
from EEG import colorElectrodes

import numpy as np

import re


import math

def hsv2rgb(h, s, v):
    h = float(h)
    s = float(s)
    v = float(v)
    h60 = h / 60.0
    h60f = math.floor(h60)
    hi = int(h60f) % 6
    f = h60 - h60f
    p = v * (1 - s)
    q = v * (1 - f * s)
    t = v * (1 - (1 - f) * s)
    r, g, b = 0, 0, 0
    if hi == 0: r, g, b = v, t, p
    elif hi == 1: r, g, b = q, v, p
    elif hi == 2: r, g, b = p, v, t
    elif hi == 3: r, g, b = p, q, v
    elif hi == 4: r, g, b = t, p, v
    elif hi == 5: r, g, b = v, p, q
    r, g, b = int(r * 255), int(g * 255), int(b * 255)
    return r, g, b

def getTransformedCoord(cx, cy, transform):
    coord = array(map(float, [cx, cy, 1.0]))
    
    mat = map(float, re.split(r'\(|\)|,', transform)[1:-1])
    mat = reshape(mat, (2, 3), order='F')
    mat = reshape(mat, (6))
    mat = concatenate((mat, [0, 0, 1]))
    mat = reshape(mat, (3, 3))
    
    return dot(mat, coord)[:-1]



def getHeadCircle():
    cx, cy, rx, ry, transform = colorElectrodes.getHeadCircleProp()

    r = arange(0, 2.0*pi+0.1, 0.1)
    xl = float(rx)*cos(r) + float(cx)
    yl = float(ry)*sin(r) + float(cy)

    mat = map(float, re.split(r'\(|\)|,', transform)[1:-1])
    mat = reshape(mat, (2, 3), order='F')
    mat = reshape(mat, (6))
    mat = concatenate((mat, [0, 0, 1]))
    mat = reshape(mat, (3, 3))

    res = zeros((len(xl), 2))
    for x, y, i in zip(xl, yl, range(len(xl))):
         #print mat.shape, array([x, y, 1.0]).shape
         res[i, :] = dot(mat, array([x, y, 1.0]))[:-1]

    return res

#listeElectrodes = ['Fp1', 'Fp2', 'F3', 'F4', 'F7', 'F8', 'C3', 'C4', 'P3', 
 #                           'P4', 'O1', 'O2', 'T3', 'T4', 'T5', 'T6', 'Fz', 'Cz', 'Pz']


def showColorMap(listeElectrodes, zVal):


    zz = getZZ(listeElectrodes, zVal)

    import matplotlib.image as mpimg
    import pylab
    
    circle = getHeadCircle()
    pylab.plot(circle[:, 0], circle[:, 1], 'k')    
            
    img=mpimg.imread(os.path.dirname(os.path.realpath(__file__))  + "\\eeg_electrodes_10-20_small.png")  
    
    #extent: [ None | scalars (left, right, bottom, top) ]
    pylab.imshow(zz, origin="upper", extent=[0.025, 0.97, -0.015, 0.935])    
    pylab.colorbar()
    pylab.imshow(img, extent=[0.0, 1.0, 0.0, 1.0])    
    pylab.show()    



def getZZ(listeElectrodes, zVal):

    cx, cy, transform = colorElectrodes.getElectrodeCoordinates(listeElectrodes)
    
    x = []
    y = []
    z = []
    for electrode in listeElectrodes:
        transformedCoord = getTransformedCoord(cx[electrode], cy[electrode], transform[electrode])
        x.append(transformedCoord[0])
        y.append(transformedCoord[1])
        z.append(zVal[electrode])
       

    points  = array([x,y]).T
        
    circle = getHeadCircle()
    headZ = griddata(points, array(z), (circle[:, 0], circle[:, 1]), method='nearest')        
        
    points  = array([x,y]).T

    points2 = array([concatenate((x, circle[:, 0])), concatenate((y, circle[:, 1]))]).T
        
    xnew = arange(min(x)-100, max(x)+100, 1)
    ynew = arange(min(y)-100, max(y)+100, 1)
    xx, yy = meshgrid(xnew, ynew)
    
    return griddata(points2, concatenate((z, headZ)), (xx, yy), method='cubic')
            




def computeColorMapPval(listeElectrodes, zVal, pVal, minZZ=-0.02, maxZZ=0.02):
    zz = getZZ(listeElectrodes, zVal)
    pp = getZZ(listeElectrodes, pVal)
    
    zz[where(isnan(zz))] = 0.0
    pp[where(isnan(pp))] = 0.0

    import matplotlib.image as mpimg
    import pylab
    
    circle = getHeadCircle()
    pylab.plot(circle[:, 0], circle[:, 1], 'k')    
            
    img=mpimg.imread(os.path.dirname(os.path.realpath(__file__))  + "\\eeg_electrodes_10-20_small.png")  
    

    M = zz.shape[0]
    N = zz.shape[1]
    Mat = zeros((M, N, 3))    

    deltaZZ = maxZZ-minZZ

    print np.min(zz), np.max(zz)

    #print minZZ, maxZZ    
    
    #s = -(pp -0.001)/(0.2-0.001) + 1.0
    s = (pp -0.2)/(0.8-0.2)
    h = (zz-minZZ)/deltaZZ*240.0
    v = 1.0
    for i1 in range(M):
        for i2 in range(N):
            r, g, b = hsv2rgb(min(max(h[i1, i2], 0.0), 240.0), min(max(s[i1, i2], 0.0), 1.0), v)
            Mat[i1, i2, 0] = float(r)/255.0
            Mat[i1, i2, 1] = float(g)/255.0
            Mat[i1, i2, 2] = float(b)/255.0
 

    #extent: [ None | scalars (left, right, bottom, top) ]
    pylab.imshow(Mat, origin="upper", extent=[0.025, 0.97, -0.015, 0.935])       
    #pylab.colorbar()
    pylab.imshow(img, extent=[0.0, 1.0, 0.0, 1.0])    






def computeColorMap(listeElectrodes, zVal, minZZ=0.00, maxZZ=0.02):
    zz = getZZ(listeElectrodes, zVal)
    
    zz[where(isnan(zz))] = 0.0
    #s = zz*0.0+1.0
    #s[where(isnan(zz))] = 0.0
    #print s

    import matplotlib.image as mpimg
    import pylab
    
    circle = getHeadCircle()
    pylab.plot(circle[:, 0], circle[:, 1], 'k')    
            
    img=mpimg.imread(os.path.dirname(os.path.realpath(__file__))  + "\\eeg_electrodes_10-20_small_whiteContour.png")  
    

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
 

    #extent: [ None | scalars (left, right, bottom, top) ]
    pylab.imshow(Mat, origin="upper", extent=[0.025, 0.97, -0.015, 0.935])       
    #pylab.colorbar()
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






def transformCircle((cx, cy, rx, ry, transform)):

    r = array([0.0, 0.5])*pi
    xl = float(rx)*cos(r) + float(cx)
    yl = float(ry)*sin(r) + float(cy)

    mat = map(float, re.split(r'\(|\)|,', transform)[1:-1])
    mat = reshape(mat, (2, 3), order='F')
    mat = reshape(mat, (6))
    mat = concatenate((mat, [0, 0, 1]))
    mat = reshape(mat, (3, 3))

    res = zeros((len(xl), 2))
    for x, y, i in zip(xl, yl, range(len(xl))):
         #print mat.shape, array([x, y, 1.0]).shape
         res[i, :] = dot(mat, array([x, y, 1.0]))[:-1]
         
    cx = res[1, 0]
    rx = abs(res[0, 0] - cx)
    cy = res[0, 1]
    ry = abs(res[1, 1] - cy)
    res[:, 1]

    return cx, cy, rx, ry




class headDrawing:
    def __init__(self):
        from matplotlib import patches
        
        self.listeElectrodes = ['Fp1', 'Fp2', 'F3', 'F4', 'F7', 'F8', 'C3', 'C4', 'P3', 
                            'P4', 'O1', 'O2', 'T3', 'T4', 'T5', 'T6', 'Fz', 'Cz', 'Pz']    
  
        self.cx = {}
        self.cy = {}
        self.rx = {} 
        self.ry = {}
        self.electrodePatch = {}
        for elect in self.listeElectrodes :
            self.cx[elect], self.cy[elect], self.rx[elect], self.ry[elect] = transformCircle(colorElectrodes.getElectCircleProp(elect))  
            self.electrodePatch[elect] = patches.Ellipse((self.cx[elect], self.cy[elect]), 2.0*self.rx[elect], 
                                                                              2.0*self.ry[elect], angle=0.0, fc="white")
    
        self.arrowPatch = []
        
        self.ColorMap = None
        self.s        = None
    
    
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


    def addColorMap(self, listeElectrodes, zVal, minZZ=0.00, maxZZ=0.02, sVal=None, smin=0.0, smax=1.0):
        
        valDict = {}      
        for elect, val in zip(listeElectrodes, zVal):
            valDict[elect] = val
        
        #print valDict
        zz = getZZ(listeElectrodes, valDict)        
            
        M = zz.shape[0]
        N = zz.shape[1]

        self.ColorMap = zeros((M, N, 3))    
    
        deltaZZ = maxZZ-minZZ
    
        h = (zz-minZZ)/deltaZZ            
            
            
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

     
        v = 1.0
        for i1 in range(M):
            for i2 in range(N):
                if isnan(zz[i1, i2]):
                    self.ColorMap[i1, i2, 0] = 1.0
                    self.ColorMap[i1, i2, 1] = 1.0
                    self.ColorMap[i1, i2, 2] = 1.0                    
                else:             
                    if s is None:
                        sij = 1.0
                    else:
                        sij = min(max(s[i1, i2], 0.0), 1.0)   
                    hij = min(max(h[i1, i2], 0.0), 1.0)   
                    r, g, b = hsv2rgb(hij*240.0, sij, v)
                    self.ColorMap[i1, i2, 0] = float(r)/255.0
                    self.ColorMap[i1, i2, 1] = float(g)/255.0
                    self.ColorMap[i1, i2, 2] = float(b)/255.0
     


    
    def plot(self, save=False, filename="fig.svg"):
            
        if save:
            import matplotlib
            matplotlib.use('Agg')
            
        import matplotlib.image as mpimg
            
        import pylab
    
    
        cx, cy, rx, ry = transformCircle(colorElectrodes.getHeadCircleProp())
        
        pylab.figure(figsize=(10, 10))         
        #thismanager = pylab.get_current_fig_manager()
        #thismanager.window.wm_geometry("500x500+0+0")     
    
        if not self.ColorMap is None:
            pylab.imshow(self.ColorMap, origin="upper", extent=[cx-rx-52, cx+rx+50, cy+ry+44, cy-ry-43])  
            if not self.s is None:
                pylab.contour(self.s, [-1, 0.95], colors='k',  origin="upper", linestyles="dashed", 
                              linewidths=3.0, extent=[cx-rx-52, cx+rx+50, cy+ry+44, cy-ry-43])
             
    
        img=mpimg.imread(os.path.dirname(os.path.realpath(__file__))  + "\\eeg_electrodes_10-20_small_clear.png")  
        pylab.imshow(img, extent=[cx-rx-65, cx+rx+65, cy+ry+27, cy-ry-80.5])  
    
        ax=pylab.gca()
    
        for elect in self.listeElectrodes : 
            ax.add_patch(self.electrodePatch[elect] )
    

            ax.text(self.cx[elect], self.cy[elect], elect,
                        horizontalalignment='center',
                        verticalalignment='center',
                        fontsize=14, color='black')

    
    
    
        for arrow in self.arrowPatch:
            ax.add_patch(arrow)
                
        pylab.axis('off')
        pylab.show()

        if save:
            pylab.savefig(filename,dpi=100, transparent=False, bbox_inches='tight', pad_inches=0)


