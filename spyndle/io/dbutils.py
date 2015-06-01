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
import numpy as np
from sqlalchemy.util._collections import KeyedTuple
from warnings import warn

def rows2df(rows, columnsIn=[], columnsOut=[]):
    """
    Convert a set of rows obtained with session.query.all() in a pandas DataFrame 
    object. If a list a column names is passed as the columnsIn argument, only
    these columns will be included. If columnsOut contains column names, these
    will be excluded.
    """     
    if len(rows) == 0:
        return pd.DataFrame()

    if len(columnsIn) != 0 and len(columnsOut) != 0:
        warn("The columnsIn shadows the columnsOut argument in rows2df.")
    
    def setColumn(d, dKey, table, columnName):
        if len(columnsIn) != 0:
            if columnName in columnsIn:
                if dKey in d:
                    d[dKey].append(getattr(table, columnName))
                else:
                    d[dKey] = [getattr(table, columnName)] 
            return
        else:
            if not columnName in columnsOut:
                if dKey in d:
                    d[dKey].append(getattr(table, columnName))
                else:
                    d[dKey] = [getattr(table, columnName)]         
        
        
    d = {}
    if isinstance(rows[0], KeyedTuple):
        """
        KeyedTuple are received from query of type query(Class1, Class2, ..., ClassN)
        """
        
        """
        Finding the column that have the same name accross different tables.
        These names cannot be conveyed directly to the DataFrame because this
        column would have twice (if two columns has the same name) as many
        items as the other columns, half of it comming from the column of this
        name in each of the two tables. Thus we identify columns names that are
        the same accross more than one table such that we can prefix them with
        the table name.
        """
        columnNames = []
        for table in rows[0]:
            columnNames.extend([c.name for c in table.__table__.columns])
        columnNames = np.sort(columnNames)
        commonNames = np.unique(columnNames[columnNames[1:] == columnNames[:-1]])      


        for row in rows:
            for table in row:
                tblName = table.__table__.name
                for column in table.__table__.columns:
                    if column.name in commonNames: 
                        setColumn(d, tblName + "_" + column.name, table, column.name)
                    else:                     
                        setColumn(d, column.name, table, column.name)                        
                
                
    else:
        for row in rows:
            for column in row.__table__.columns:
                setColumn(d, column.name, row, column.name)  

    return pd.DataFrame(d)         
