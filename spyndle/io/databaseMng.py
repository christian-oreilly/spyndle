# -*- coding: utf-8 -*-

"""
    Code used for managing database interaction.
    
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


import sqlalchemy as sa
from spyndle.io import Session, Base, DataModelMng
from sqlalchemy.pool import NullPool, QueuePool


def clearDatabase(dbName):
    """
    Helper function to easily clear a given database without the needs to 
    create a DatabaseMng object. This function keeps the structure
    of the database but erase any data in it.
    """    
    dbMng = DatabaseMng()
    dbMng.connectDatabase(dbName)       
    dbMng.clearDatabase()
    dbMng.disconnectDatabase()





class DatabaseMng():
    """
    Class used to manage the interaction with the SQL database.
    """
    

    def __init__(self, dbName = "", shard=None):
        """
        DatabaseMng constructor.
         
        If a name dbName is passed, connect to the database and create the 
        tablesof the data model if those are not already existent. Else, only 
        set the session to none and do nothing.
        """ 
        self.session = None
        if dbName != "":
            if shard is None:
                self.connectDatabase(dbName)
            else:
                shard = shard.replace(":", "_").replace("/", "_").replace("\\", "_").replace(".", "_")
                self.connectDatabase(dbName + "___shard_" + shard)
            self.createTables()
            if not self.isConnected():
                raise IOError("Error connecting to the database.")   
                
        self.dmm = DataModelMng(self.session)  
        
            
            
    def __enter__(self):
        return self
    def __exit__(self, type, value, traceback):
        self.disconnectDatabase()   
  
      
     
    @property
    def dmm(self):
        """
         Property providing access to the data model manager object.
        """               
        return self.__dmm 

    @dmm.setter
    def dmm(self, dmm):
        self.__dmm = dmm
    
    
    

    def createTables(self):
        """
        Create the tables of the data model if those are not already existent.
        """            
        Base.metadata.create_all(self.dbEngine)   


     
    def connectDatabase(self, URL):        
        """
        Connect to the database URL. This creates the session object which can
        be used to performs operation on this database.
        
        URL is of the form dialect+driver://user:password@host/dbname[?key=value..]
        or is an URL.
        """         
        if URL.split(":")[0].split("+")[0] == "sqlite": 
            self.dbEngine = sa.create_engine(URL, poolclass=QueuePool, pool_size=20,  
                                             pool_timeout=180)#, isolation_level='READ UNCOMMITTED') #NullPool)
        elif URL.split(":")[0].split("+")[0] == "mysql":
            dbName = URL.split("@")[1].split("/")[1].split("[")[0]
            connectStr = URL.split("@")[0] + "@" + URL.split("@")[1].split("/")[0]
            self.dbEngine = sa.create_engine(connectStr, poolclass=QueuePool, pool_size=20,  
                                             pool_timeout=180)#, isolation_level='READ UNCOMMITTED') #NullPool)            
            self.dbEngine.execute("CREATE DATABASE IF NOT EXISTS " + dbName) #create db
            self.dbEngine.execute("USE " + dbName) # select new db
        else:
            raise ValueError("Only sqlite and mysql are acccepted as database type for now.")
            
            
        Session.configure(bind=self.dbEngine)   
        self.session = Session()            
    
    def clearDatabase(self):
        """
        Clear the connected database. That is, keep the structure
        of the database but erase any data in it.
        """            
        if not self.session is None:

            for name, table in Base.metadata.tables.items(): 
                print table.delete() 
                self.session.execute(table.delete()) 
    
            self.session.commit() 
        else:
            raise IOError("You must fist connect to a database before calling "\
                          "the DatabaseMng.clearDatabase(...) method. "\
                          "Alternatively, you can call the clearDatabase "\
                          "function, passing the database name as parameter.")
        
        
        
             
    def disconnectDatabase(self):
        """
        Disconnect from the currently connected database.
        """           
        if not self.session is None:        
            self.session.close()
            self.session = None
            
     
   
    def isConnected(self):
        """
        Check if a the database managed is currently connected to a database.
        """            
        return not self.session is None
      
        

    def __del__(self):
        self.disconnectDatabase()
