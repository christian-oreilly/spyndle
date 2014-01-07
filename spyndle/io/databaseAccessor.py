# -*- coding: utf-8 -*-
"""
Created on Mon Jan 06 15:53:55 2014

@author: oreichri
"""

from dataAccessor import DataAccessor
from databaseMng import DatabaseMng
from dbutils import rows2df
from dataModel import TransientEvent

import sqlalchemy as sa    
import pandas as pd    

    
    
class DatabaseAccessor(DataAccessor):
    """
    Subclass of DataAccessor specific for accessing recording whitin a database.
    """
    
    def __init__(self):
        #TODO: The sharding must be implemented correctly before such a class  
        #      can be implemented appropriately.
        pass
    


    def getSpindleData(self, nights=[], eventNames=[],  
                              dbName = "sqlite:///:memory:", shards=[None]):
        
        if not isinstance(shards, list):
            raise TypeError("The argument shards must be a list.")
     
        DFs = []
        for shard in shards:            
            try:
                dbMng = DatabaseMng(dbName, shard=shard)
            except sa.exc.ArgumentError, error:
                 print "Error connecting to the specified database URL. "\
                       "The format of the database path is invalid."     \
                       "\n\nSQLAlchemy error message: " + str(error)
    
            query = dbMng.session.query(TransientEvent)
            if nights != []:
                query = query.filter(TransientEvent.psgNight.in_(nights))   
        
            if eventNames != []:
                query = query.filter(TransientEvent.eventName.in_(eventNames))  

            DFs.append(rows2df(query.all()))
        
        return pd.concat(DFs)

