# -*- coding: utf-8 -*-
"""
Created on Mon Jan 06 15:53:55 2014

@author: oreichri
"""

from .dataAccessor import DataAccessor
from .databaseMng import DatabaseMng
from .dbutils import rows2df
from .dataModel import TransientEvent, SpindleEvent, SlowWaveEvent, Propagation
from urllib.parse import urlparse, urlunparse

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


    def getShards(self, nights, dbName, shards):
        parsedUrl = urlparse(dbName)

        if parsedUrl.scheme == "mysql":
            if shards == [None]:
                if nights == []:
                    """
                    If no shard or nights have been passed in parameter and we
                    are using an mysql database, look if there is databases
                    with the with name staring with dbName. If so, list these
                    databases as shards of dbName.
                    """
                    # Remove the first caracter which is a "/"
                    shardPattern = parsedUrl.path[1:] + "___shard_"
    
                    msqlConnectString = list(urlparse(dbName))
                    msqlConnectString[2] = ""  # Removing the url path
                    msqlConnectString = urlunparse(msqlConnectString)
                    engine = sa.create_engine(msqlConnectString)
                    databases = sa.inspect(engine).get_schema_names()
                    shards = [db[len(shardPattern):] for db in databases if
                                 db[:len(shardPattern)] == shardPattern]
                    if len(shards) == 0:
                        if parsedUrl.path[1:] in databases:
                            shards = [None]   
                else:
                    shards = nights

        return shards
        
        

    def getSlowWaveData(self, nights=[], eventNames=[],
                        dbName="sqlite:///:memory:", shards=[None],
                        filterIn={}, filterOut={}, columnsIn=[], columnsOut=[]):

        return self.getTransientData(self.getSlowWaveQuery, [SlowWaveEvent], 
                                     nights, eventNames, dbName, shards, 
                                     filterIn,  filterOut, columnsIn, 
                                     columnsOut)



    def getSpindleData(self, nights=[], eventNames=[],
                        dbName="sqlite:///:memory:", shards=[None],
                        filterIn={}, filterOut={}, columnsIn=[], columnsOut=[]):

        return self.getTransientData(self.getSpindleQuery, [SpindleEvent], 
                                     nights, eventNames, dbName, shards, 
                                     filterIn,  filterOut, columnsIn, 
                                     columnsOut)



    def getPropagationData(self, nights=[], eventNames=[],
                        dbName="sqlite:///:memory:", shards=[None],
                        filterIn={}, filterOut={}, columnsIn=[], columnsOut=[]):

        return self.getTransientData(self.getPropagationQuery, [Propagation, SpindleEvent], 
                                     nights, eventNames, dbName, shards, 
                                     filterIn,  filterOut, columnsIn, 
                                     columnsOut)




    def getSlowWaveQuery(self, dbMng):
        return dbMng.session.query(TransientEvent, SlowWaveEvent)\
                            .join(SlowWaveEvent, TransientEvent.ID == SlowWaveEvent.ID)

    def getSpindleQuery(self, dbMng):
        return dbMng.session.query(TransientEvent, SpindleEvent)\
                            .join(SpindleEvent, TransientEvent.ID == SpindleEvent.ID)


    def getPropagationQuery(self, dbMng, EventClass=SpindleEvent):

        return dbMng.session.query(TransientEvent, EventClass, Propagation)\
                            .join(EventClass, TransientEvent.ID == EventClass.ID)\
                            .join(Propagation, TransientEvent.ID == Propagation.transientEventID)




    def getTransientData(self, queryFct, EventClasses, nights=[], eventNames=[],
                        dbName="sqlite:///:memory:", shards=[None],
                        filterIn={}, filterOut={}, columnsIn=[], columnsOut=[]):

        if not isinstance(shards, list):
            raise TypeError("The argument shards must be a list.")


        shards = self.getShards(nights, dbName, shards)
        
        DFs = []
        for shard in shards:     
            try: 
                with DatabaseMng(dbName, shard=shard) as dbMng:
            
                    query = queryFct(dbMng)
                    
                    if nights != []:
                        query = query.filter(TransientEvent.psgNight.in_(nights))   
             
                    if eventNames != []:
                        query = query.filter(TransientEvent.eventName.in_(eventNames))  
        
                    for key in filterIn:
                        for obj in EventClasses:
                            if key in dir(obj):
                                if isinstance(filterIn[key], list):
                                    query = query.filter(getattr(obj, key).in_(filterIn[key]))
                                else:
                                    query = query.filter(getattr(obj, key) == filterIn[key]) 
                                
                    for key in filterOut:
                        for obj in EventClasses:
                            if key in dir(obj):
                                if isinstance(filterOut[key], list):
                                    query = query.filter(sa.not_(getattr(obj, key).in_(filterOut[key])))
                                else:
                                    query = query.filter(getattr(obj, key) != filterOut[key]) 


                    DFs.append(rows2df(query.all(), columnsIn, columnsOut))
            except sa.exc.ArgumentError as error:
                 print(("Error connecting to the specified database URL. "\
                       "The format of the database path is invalid."     \
                       "\n\nSQLAlchemy error message: " + str(error)))
                 raise
        
        return pd.concat(DFs)




