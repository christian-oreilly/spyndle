# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 19:54:09 2013

@author: oreichri
"""

import sqlalchemy as sa
import sys
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
    channelName         = sa.Column(sa.String, sa.ForeignKey("channel.name"))
    eventName           = sa.Column(sa.String)
    
    RMSamp              = sa.Column(sa.Float)
    meanFreq            = sa.Column(sa.Float)
    stage               = sa.Column(sa.String)
    cycle               = sa.Column(sa.Integer)
    slopeOrigin         = sa.Column(sa.Float)
    slope               = sa.Column(sa.Float)
    filteredRMSamp      = sa.Column(sa.Float)
    
    #channel             = sa.orm.relationship("Channel", backref="transientEvent")    

    #propagagations      = sa.orm.relationship("Propagation")
              
def buildTransientEvent(event, psgNight):
             
    
    startTime           = event.startTime
    duration            = event.duration()
    channelName         = event.channel
    
    eventName           = event.name
    
    RMSamp              = float(event.properties["RMSamp"]) if "RMSamp" in event.properties else None
    meanFreq            = float(event.properties["meanFreq"]) if "meanFreq" in event.properties else None
    stage               = event.properties["stage"] if "stage" in event.properties else None
    cycle               = int(event.properties["cycle"]) if "cycle" in event.properties else None
    slopeOrigin         = float(event.properties["slopeOrigin"]) if "slopeOrigin" in event.properties else None
    slope               = float(event.properties["slope"]) if "slope" in event.properties else None
    filteredRMSamp      = float(event.properties["filteredRMSamp"]) if "filteredRMSamp" in event.properties else None

    return TransientEvent(psgNight = psgNight, startTime = startTime, duration = duration, 
                          channelName = channelName, RMSamp = RMSamp, meanFreq = meanFreq,
                          stage = stage, cycle = cycle, slopeOrigin = slopeOrigin,
                          slope = slope, filteredRMSamp= filteredRMSamp, eventName=eventName)
         
             
class PSGNight(Base): 
    __tablename__       = "psgNight"
    
    fileName            = sa.Column(sa.String, primary_key=True)



class Channel(Base): 
    __tablename__       = "channel"
    
    name                = sa.Column(sa.String, primary_key=True)
    #activeElectrode = sa.Column(sa.String)
    #reference       = sa.Column(sa.String)




class DataManipulationProcess(Base): 
    __tablename__       = "dataManipulationProcess"
    
    no                  = sa.Column(sa.Integer, primary_key=True)
    name                = sa.Column(sa.String, primary_key=True)
    reprString          = sa.Column(sa.String)  
    datetime            = sa.Column(sa.DateTime)    
    
    #activeElectrode = sa.Column(sa.String)
    #reference       = sa.Column(sa.String)





             
             
             
             
             

class SPF(Base): 
    __tablename__       = "spf"
    no                  = sa.Column(sa.Integer, primary_key=True)
    psgNight            = sa.Column(sa.String, sa.ForeignKey("psgNight.fileName"))


"""
class SPF_propagation(Base): 
    __tablename__       = "spf_propagation"
    __table_args__      = (sa.UniqueConstraint("noSPF", "noPropagation"), )
    
    noSPF               = sa.Column(sa.Integer, sa.ForeignKey("spf.no"))
    noPropagation       = sa.Column(sa.Integer, sa.ForeignKey("propagation.no"))
"""


class Propagation(Base): 
    __tablename__       = "propagation"

    no                  = sa.Column(sa.Integer, primary_key=True)
    
    # ID of the spindle from which this propagation has been computed.
    spindleNo           = sa.Column(sa.Integer, sa.ForeignKey("transientEvent.no"))
    
    sourceChannelName   = sa.Column(sa.String, sa.ForeignKey("channel.name"))
    sinkChannelName     = sa.Column(sa.String, sa.ForeignKey("channel.name"))
    
    similarity          = sa.Column(sa.Float)
    delay               = sa.Column(sa.Float)
    offset              = sa.Column(sa.Float)
    
    transientEvent      = sa.orm.relationship("TransientEvent", #backref="propagations")
                                                primaryjoin='Propagation.spindleNo==TransientEvent.no',
                                                lazy='joined')    
    
    
    
    #startTime           = sa.orm.relationship("TransientEvent.startTime")
    #testChannel         = sa.orm.relationship("Channel", backref="propagations")

    ######################
    # Rejection criteria  

    # Is the propagation is a false positive according to a threhold on the FDR?
    # (criterion C1)  
    isFP                = sa.Column(sa.Boolean)

    # Is the propagation an outlier according to Tuckey's test?  
    # (criterion C2)  
    isOutlier           = sa.Column(sa.Boolean)


    # Used in SPF computation.
    bidirect            = sa.Column(sa.Integer)
    source              = sa.Column(sa.Integer, sa.ForeignKey("propagation.no"))
    noSPF               = sa.Column(sa.Integer, sa.ForeignKey("spf.no"))
    inverted            = sa.Column(sa.Boolean)


    def __init__(self, *args, **kwargs):

        super(Propagation, self).__init__(*args, **kwargs)  # call parent class init        
        
        # By default, propagations are considered unreliable. They must be proven
        # reliable.
        self.isFP       = True
        self.isOutlier  = True
        
        
        self.noSPF      = -1
        self.bidirect   = -1
        self.source     = -1
        self.inverted   = False
        






        
class DynamicProperty(Base): 
    __tablename__       = "dynamicProperty"

    name                = sa.Column(sa.String, primary_key=True)        
    type                = sa.Column(sa.String)    
        
        
class PropRelProperty(Base): 
    __tablename__       = "propagationRelationshipProperty"

    noPropRel           = sa.Column(sa.Integer, sa.ForeignKey("propagationRelationship.no"), primary_key=True)
    propertyName        = sa.Column(sa.String,  sa.ForeignKey("dynamicProperty.name"), primary_key=True)
    value               = sa.Column(sa.String)        

    def __repr__(self):
        return "{noPropRel:'%s', propertyName:'%s', value:'%s'}" % (self.noPropRel, self.propertyName, self.value)


class PropagationRelationship(Base): 
    __tablename__       = "propagationRelationship"
    __table_args__      = (sa.UniqueConstraint("sourceChannelName", 
                                               "sinkChannelName",
                                               "psgNight",
                                               "eventName"), )

    no                  = sa.Column(sa.Integer, primary_key=True)
    sourceChannelName   = sa.Column(sa.String, sa.ForeignKey("channel.name"))
    sinkChannelName     = sa.Column(sa.String, sa.ForeignKey("channel.name"))
    psgNight            = sa.Column(sa.String, sa.ForeignKey("psgNight.fileName"))
    eventName           = sa.Column(sa.String)
    
    cutoff              = sa.Column(sa.Float)
    FDR                 = sa.Column(sa.Float)
            

    def getDict(self, session):
        retDict = {"sourceChannelName":self.sourceChannelName,
                   "sinkChannelName":self.sinkChannelName,
                   "psgNight":self.psgNight,
                   "eventName":self.eventName,
                   "cutoff":self.cutoff,
                   "FDR":self.FDR}
                   
        properties = session.query(PropRelProperty, DynamicProperty.type)\
                            .join(DynamicProperty, DynamicProperty.name == PropRelProperty.propertyName)\
                            .filter(PropRelProperty.noPropRel == self.no)\
                            .all()
                            
        for (aProperty, aType) in properties:
            try:
                # Check if it's a builtin type
                cls = getattr(__builtins__, aType)
            except AttributeError:
                moduleName, typeName = aType.split(".")
                module = sys.modules[moduleName]
                cls = getattr(module, typeName)
                
            retDict[aProperty.propertyName] = cls(aProperty.value)
            
        return retDict
                   
                   
    def update(self, session, row, behavior="raise"):
        

        if behavior == "addMissingFields":
            for key in row:
                if session.query(DynamicProperty).filter_by(name = key).count() == 0:
                    session.add(DynamicProperty(name=key, type=str(type(row[key])).split("'")[1]))                    
                   
                        

        no = session.query(PropagationRelationship)\
                        .filter_by(sourceChannelName = self.sourceChannelName,
                                   sinkChannelName   = self.sinkChannelName,
                                   psgNight          = self.psgNight,
                                   eventName         = self.eventName )\
                        .one().no        
        
        
        session.add_all([PropRelProperty(noPropRel=no, propertyName=key, value=str(row[key])) for key in row]) 
        session.commit()


    def add(self, session, behavior = "raise"):
        
        try:
            session.add(self)
            session.commit()
            
        except sa.exc.IntegrityError:
            
            if behavior == "raise":
                raise
             
            session.rollback()
            # updateSilently keep the same primaryKey. We could implement
            # also replaceSilently which whould also change the primary key...
            if behavior == "updateSilently":
                
                session.query(PropagationRelationship)\
                            .filter_by(sourceChannelName = self.sourceChannelName,
                                       sinkChannelName   = self.sinkChannelName,
                                       psgNight          = self.psgNight,
                                       eventName         = self.eventName )\
                            .update({"cutoff": self.cutoff,
                                     "FDR"   : self.FDR})
                session.commit()


            elif behavior == "failSilently":
                pass
            
            else:
                raise ValueError("The value '" + behavior + "' is invalid for the behavior variable.")













"""
  Helper function to easily clear a given database.
"""
def clearDatabase(dbName):
    dbMng = DatabaseMng()
    dbMng.connectDatabase(dbName)       
    dbMng.clearDatabase()
    dbMng.disconnectDatabase()



class DatabaseMng():

    def __init__(self, dbName = ""):
        if dbName == "":
            self.session = None
        else:
            self.connectDatabase(dbName)
            self.createTables()
            if not self.isConnected():
                raise IOError("Error connecting to the database.")             
    
    
    def createTables(self):
        Base.metadata.create_all(self.dbEngine)   

    
    def connectDatabase(self, dbName):        
        self.dbEngine = create_engine(dbName)
        Session.configure(bind=self.dbEngine)   
        self.session = Session()
        
        
    def clearDatabase(self):
        if not self.session is None:

            for name, table in Base.metadata.tables.items(): 
                print table.delete() 
                self.session.execute(table.delete()) 
    
            self.session.commit() 
        else:
            print "You must fist connect to the database."
        
        
        
        
    def disconnectDatabase(self):
        self.session.close()
        self.session = None
        
        
    def isConnected(self):
        return not self.session is None
        
        
    
        



