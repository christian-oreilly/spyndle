# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 19:54:09 2013

@author: oreichri
"""

import sqlalchemy as sa
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine

from spyndle.io import Session


# Getting an SQLAchemy base class
Base = declarative_base()

# Creating the ORM database classes 

# TODO: The definition should depend on the type of transient event analyzed.
class TransientEvent(Base): 
    __tablename__       = "transientEvent"
    
    no                  = sa.Column(sa.Integer, primary_key=True)
    psgNight            = sa.Column(sa.String, sa.ForeignKey("psgNight.fileName"))
    startTime           = sa.Column(sa.Float)
    duration            = sa.Column(sa.Float)
    channel             = sa.Column(sa.String)
    eventName           = sa.Column(sa.String)
    
    RMSamp              = sa.Column(sa.Float)
    meanFreq            = sa.Column(sa.Float)
    stage               = sa.Column(sa.String)
    cycle               = sa.Column(sa.Integer)
    slopeOrigin         = sa.Column(sa.Float)
    slope               = sa.Column(sa.Float)
    filteredRMSamp      = sa.Column(sa.Float)

    #propagagations      = sa.orm.relationship("Propagation")
              
def buildTransientEvent(event, psgNight):
             
    
    startTime           = event.startTime
    duration            = event.duration()
    channel             = event.channel
    eventName           = event.name
    
    RMSamp              = float(event.properties["RMSamp"]) if "RMSamp" in event.properties else None
    meanFreq            = float(event.properties["meanFreq"]) if "meanFreq" in event.properties else None
    stage               = event.properties["stage"] if "stage" in event.properties else None
    cycle               = int(event.properties["cycle"]) if "cycle" in event.properties else None
    slopeOrigin         = float(event.properties["slopeOrigin"]) if "slopeOrigin" in event.properties else None
    slope               = float(event.properties["slope"]) if "slope" in event.properties else None
    filteredRMSamp      = float(event.properties["filteredRMSamp"]) if "filteredRMSamp" in event.properties else None

    return TransientEvent(psgNight = psgNight, startTime = startTime, duration = duration, 
                          channel = channel, RMSamp = RMSamp, meanFreq = meanFreq,
                          stage = stage, cycle = cycle, slopeOrigin = slopeOrigin,
                          slope = slope, filteredRMSamp= filteredRMSamp, eventName=eventName)
         
             

class PSGNight(Base): 
    __tablename__       = "psgNight"
    
    fileName            = sa.Column(sa.String, primary_key=True)

class Propagation(Base): 
    __tablename__       = "propagation"

    no                  = sa.Column(sa.Integer, primary_key=True)
    spindleNo           = sa.Column(sa.Integer, sa.ForeignKey("transientEvent.no"))
    testChannel         = sa.Column(sa.String)
    similarity          = sa.Column(sa.Float)
    delay               = sa.Column(sa.Float)
    offset              = sa.Column(sa.Float)
    
    spindle             = sa.orm.relationship("TransientEvent", backref="propagations")





class DatabaseMng():

    def __init__(self):
        self.session = None
    
    
    def createDatabase(self):
        print "Creating database..."
        Base.metadata.create_all(self.dbEngine)   

    
    def connectDatabase(self, dbPath):        
        self.dbEngine = create_engine(dbPath)
        Session.configure(bind=self.dbEngine)   
        self.session = Session()
        
        
    def disconnectDatabase(self):
        self.session.close()
        self.session = None
        
        
    def isConnected(self):
        return not self.session is None
        
        
    
        
        
        