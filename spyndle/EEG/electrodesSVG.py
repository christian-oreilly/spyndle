 # coding=ISO-8859-1

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
 Date  : June 10, 2012

"""

import warnings
import os
from xml.dom import minidom


def readElectrodesFile(fileName):

    with open(fileName, "r") as f:
        line =  f.readline()    
        electNames = line.split()
        line =  f.readline()    
        electColors = line.split()
        line =  f.readline()    
        strokeColors = line.split()
        
    electrodes = {}
    stroke     = {}
    k = 0

    for names in electNames:
        electrodes[names] = electColors[k]
        stroke[names]     = strokeColors [k]
        k += 1

    return electrodes, stroke





def generateSVGColoredMap_fromFile(fileIn, fileOut):
    electrodeColor, strokeColor  = readElectrodesFile(fileIn)
    generateSVGColoredMap(fileOut, electrodeColor, strokeColor)
    
    



def getElectrodeCoordinates(electrodeNames = None):

    if not electrodeNames is None:
        warnings.warn("The use of the electrodeNames parameters to call electrodeSVG.getElectrodeCoordinates() is now deprecated.", DeprecationWarning, stacklevel=2)

    dom = minidom.parse(os.path.dirname(os.path.realpath(__file__))  + "\\eeg_electrodes_10-20.svg")

    SVGnode = dom.getElementsByTagName("svg")[0]   

    cx = {}
    cy = {}
    transform = {}

        
    for pathNode in SVGnode.getElementsByTagName("path"):
        if pathNode.getAttribute("id")[0:6] == "circle":
            
            nameElectrode = pathNode.getAttribute("id")[6:]
            #if nameElectrode in electrodeNames:
            
            cx[nameElectrode] = pathNode.getAttribute("sodipodi:cx")
            cy[nameElectrode] = pathNode.getAttribute("sodipodi:cy")
            transform[nameElectrode] = pathNode.getAttribute("transform")


    for small, large in zip(["T3", "T5", "T4", "T6"], ["T7", "P7", "T8", "P8"]):
        cx[small]        = cx[large]
        cy[small]        = cy[large]
        transform[small] = transform[large]

    return cx, cy, transform



def getHeadCircleProp():

    dom = minidom.parse(os.path.dirname(os.path.realpath(__file__))  + "\\eeg_electrodes_10-20_small.svg")

    SVGnode = dom.getElementsByTagName("svg")[0]   

        
    for pathNode in SVGnode.getElementsByTagName("path"):
        if pathNode.getAttribute("id") == "path3097":
            cx = pathNode.getAttribute("sodipodi:cx")
            cy = pathNode.getAttribute("sodipodi:cy")
            rx = pathNode.getAttribute("sodipodi:ry")
            ry = pathNode.getAttribute("sodipodi:ry")
            transform = pathNode.getAttribute("transform")


    return cx, cy, rx, ry, transform



def getElectCircleProp(nameElect):

    dom = minidom.parse(os.path.dirname(os.path.realpath(__file__))  + "\\eeg_electrodes_10-20.svg")

    SVGnode = dom.getElementsByTagName("svg")[0]   

    if nameElect == "T3":
        nameElect = "T7"
    elif nameElect == "T5":
        nameElect = "P7"
    elif nameElect == "T4":
        nameElect = "T8"
    elif nameElect == "T6":
        nameElect = "P8"
    
    for pathNode in SVGnode.getElementsByTagName("path"):
        if pathNode.getAttribute("id") == "circle" + nameElect:
            cx = pathNode.getAttribute("sodipodi:cx")
            cy = pathNode.getAttribute("sodipodi:cy")
            rx = pathNode.getAttribute("sodipodi:ry")
            ry = pathNode.getAttribute("sodipodi:ry")
            transform = pathNode.getAttribute("transform")
            return cx, cy, rx, ry, transform

    print(("Electrode " + nameElect + " not found in eeg_electrodes_10-20.svg."))
    raise KeyError
    

def generateSVGColoredMap(fileName, electrodeColor, strokeColor, template=None):
    
    if template is None:
        template = os.path.dirname(os.path.realpath(__file__))  + "\\eeg_electrodes_10-20_small.svg"
   
    dom = minidom.parse(template)    

    SVGnode = dom.getElementsByTagName("svg")[0]   
        
    for pathNode in SVGnode.getElementsByTagName("path"):
        if pathNode.getAttribute("id")[0:6] == "circle":
            
            nameElectrode = pathNode.getAttribute("id")[6:]
            attribut = pathNode.getAttribute("style")
            subAttributesLst = attribut.split(';')
            
            newAttritut = ""
            for subattribut in subAttributesLst:
                
                if subattribut.split(':')[0] == 'fill':     
                    if nameElectrode in electrodeColor:
                        subattribut = 'fill:' + electrodeColor[nameElectrode][0:7]
                if subattribut.split(':')[0] == 'stroke':     
                    if nameElectrode in electrodeColor:
                        subattribut = 'stroke:' + strokeColor[nameElectrode][0:7]
                if subattribut.split(':')[0] == 'stroke-width':     
                    if nameElectrode in electrodeColor:
                        subattribut = 'stroke-width:' + "5.0"
                 
                newAttritut = newAttritut + subattribut + ';'

            
            pathNode.setAttribute("style", newAttritut)
    
    
    file = open(fileName, "w")
    file.write(dom.toxml())
    file.close()
