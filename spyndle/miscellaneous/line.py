# -*- coding: utf-8 -*-

"""
    Copyright (C) 2013  Christian O'Reilly

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
 Date  : July 1, 2013

"""

from numpy import sqrt

class Point(object):
    def __init__(self, x=0, y=0):
        self.setPoint(x, y)
    
    def setPoint(self, x, y):
        self.__x = x
        self.__y = y

    def setx(self, x):
        self.__x = x
        
    def sety(self, y):
        self.__y = y

    def x(self) : return self.__x
    def y(self) : return self.__y


    def __sub__(self, p2):
        return Point(self.x()-p2.x(), self.y()-p2.y())

    def __add__(self, p2):
        return Point(self.x()+p2.x(), self.y()+p2.y())

    def __mul__(self, p2):
        return Point(self.x()*p2.x(), self.y()*p2.y())

    def __div__(self, p2):
        return Point(self.x()/p2.x(), self.y()/p2.y())




class Line(object):
    def __init__(self, p1=Point(), p2=Point()):
        self.setP1(p1)
        self.setP2(p2)
    
    #def setLine(self, x1, y1, x2, y2):
      #  self.__p1 = Point(x1, y1)
     #   self.__p2 = Point(x2, y2)        


    def setP1(self, p1):
        self.__p1 = p1
        
    def setP2(self, p2):
        self.__p2 = p2
        
    def x1(self) : return self.__p1.x()
    def y1(self) : return self.__p1.y()
    def x2(self) : return self.__p2.x()
    def y2(self) : return self.__p2.y()

    def p1(self) : return self.__p1
    def p2(self) : return self.__p2


    def length(self):
        return sqrt((self.p1().x() - self.p2().x())**2 + (self.p1().y() - self.p2().y())**2)

    # Based on Antonio, F. "Faster Line Segment Intersection. Ch. IV.6 in 
    # Graphics Gems III (Ed. D. Kirk). San Diego: Academic Press, pp. 199-202 
    # and 500-501, 1992.     
    def isIntersecting(self, l2):
            
        a = self.p2() - self.p1()
        b = l2.p1()   - l2.p2()
        c = self.p1() - l2.p1()
    
        den = a.y() *b.x() - a.x() * b.y()
    
        if den == 0:
            return False

        num1 = b.y()*c.x() - b.x()*c.y()
        
        if den > 0 :
            if num1 < 0 or num1 > den:
                return False
        elif num1 > 0 or num1 < den:
            return False

        num2 = a.x()*c.y() - a.y()*c.x()

        if den > 0 :
            if num2 <0 or num2 > den:
                return False
        elif num2 > 0 or num2 < den:
            return False

        return True    
