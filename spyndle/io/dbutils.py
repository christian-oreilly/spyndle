# -*- coding: utf-8 -*-

"""
    Module implementing database interaction utilitary functions.
    
    Copyright (C) 2013  Christian O'Reilly

    For personnal, educationnal, and research purpose, this software is 
    provided under the GNU GPL (V.3) license: you can redistribute it and/or
    modify it under the terms of the version 3 of the GNU General Public 
    License as published by the Free Software Foundation.
          
    To use this software in commercial application, please contact the author.
    
    If this code is used for research purpose, the reference [1] should be
    cited in the derived publication.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


 Author: Christian O'Reilly (christian.oreilly@umontreal.ca)
 Date  : November 1, 2013

"""


import pandas as pd


def rows2df(rows):
    """
    Convert a set of rows obtained with session.query.all() in a pandas DataFrame 
    object.
    """     
    if len(rows) == 0:
        return pd.DataFrame()
    
    d = {}
    for column in rows[0].__table__.columns:
        d[column.name] = [getattr(rows[0], column.name)]

    for row in rows[1:]:
        for column in row.__table__.columns:
            d[column.name].append(getattr(row, column.name))

    return pd.DataFrame(d)         
