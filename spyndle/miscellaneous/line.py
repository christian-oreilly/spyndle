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

try:
    from PySide import QtCore, QtGui    
except ImportError:
    from PyQt4 import QtCore, QtGui   


from numpy import sqrt

class Polygon(QtGui.QPolygon):
    def __init__(self, points):

        QtGui.QPolygon.__init__(self)        
        for point in points:
            if isinstance(point, Point):
                point = QtCore.QPoint(point.x, point.y)
            elif not isinstance(point, QtCore.QPoint):
                raise TypeError("The points argument contain an item which is equal to " \
                                + str(point) + " and of type " + str(type(point)))
        
            self.append(point)


class Point(object):
    def __init__(self, x=0, y=0):
        if isinstance(x, QtCore.QPoint):
            self.setPoint(x.x(), x.y())
        else:
            self.setPoint(x, y)
    
    def setPoint(self, x, y):
        self._x = x
        self._y = y


    def QPoint(self):
        return QtCore.QPoint(self.x, self.y)


    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, x):
        self._x = x

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y):
        self._y = y
        
    

    def __sub__(self, p2):
        if isinstance(p2, QtCore.QPoint):
            p2 = Point(p2)      
        elif not isinstance(p2, Point):
            raise TypeError("The p2 argument is equal to " \
                            + str(p2) + " and of type " + str(type(p2)))
        return Point(self.x-p2.x, self.y-p2.y)

    def __add__(self, p2):
        if isinstance(p2, QtCore.QPoint):
            p2 = Point(p2)      
        elif not isinstance(p2, Point):
            raise TypeError("The p2 argument is equal to " \
                            + str(p2) + " and of type " + str(type(p2)))
        return Point(self.x+p2.x, self.y+p2.y)

    def __mul__(self, p2):
        if isinstance(p2, QtCore.QPoint):
            p2 = Point(p2)      
        elif not isinstance(p2, Point):
            raise TypeError("The p2 argument is equal to " \
                            + str(p2) + " and of type " + str(type(p2)))
        return Point(self.x*p2.x, self.y*p2.y)

    def __div__(self, p2):
        if isinstance(p2, QtCore.QPoint):
            p2 = Point(p2)      
        elif not isinstance(p2, Point):
            raise TypeError("The p2 argument is equal to " \
                            + str(p2) + " and of type " + str(type(p2)))
        return Point(self.x/p2.x, self.y/p2.y)




class Line(object):
    def __init__(self, p1=Point(), p2=Point()):
        
        self._p1 = p1
        self._p2 = p2
    
    #def setLine(self, x1, y1, x2, y2):
      #  self._p1 = Point(x1, y1)
     #   self._p2 = Point(x2, y2)        

    @property
    def p1(self):
        return self._p1

    @p1.setter
    def p1(self, p1):
        if isinstance(p1, QtCore.QPoint):
            self._p1 = Point(p1)
        elif isinstance(p1, Point):
            self._p1 = p1
        else:
            raise TypeError
        
        
        
        
        
    @property
    def p2(self):
        return self._p2

    @p2.setter
    def p2(self, p2):
        if isinstance(p2, QtCore.QPoint):
            self._p2 = Point(p2)
        elif isinstance(p2, Point):
            self._p2 = p2
        else:
            raise TypeError
        
        
        
        
    @property
    def x1(self):
        if self._p1 is None:
            return None
        return self._p1.x      
        
    @property
    def x2(self):
        if self._p2 is None:
            return None        
        return self._p2.x       
        
        

    @property
    def y1(self):
        if self._p1 is None:
            return None        
        return self._p1.y      
        
    @property
    def y2(self):
        if self._p2 is None:
            return None        
        return self._p2.y       
        
                



    def length(self):
        return sqrt((self.x1 - self.x2)**2 + (self.y1 - self.y2)**2)

    # Based on Antonio, F. "Faster Line Segment Intersection. Ch. IV.6 in 
    # Graphics Gems III (Ed. D. Kirk). San Diego: Academic Press, pp. 199-202 
    # and 500-501, 1992.     
    def isIntersecting(self, l2):
            
        a = self.p2 - self.p1
        b = l2.p1   - l2.p2
        c = self.p1 - l2.p1
    
        den = a.y *b.x - a.x * b.y
    
        if den == 0:
            return False

        num1 = b.y*c.x - b.x*c.y
        
        if den > 0 :
            if num1 < 0 or num1 > den:
                return False
        elif num1 > 0 or num1 < den:
            return False

        num2 = a.x*c.y - a.y*c.x

        if den > 0 :
            if num2 <0 or num2 > den:
                return False
        elif num2 > 0 or num2 < den:
            return False

        return True    
