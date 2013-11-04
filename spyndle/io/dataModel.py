# -*- coding: utf-8 -*-

"""
    Implementation of the ORM database classes for using SQLAlchemy to manage
    the recording of EEG analysis in an SQL database.
    
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


 [1] O'Reilly, C. & Nielsen, T. (2013) Assessing EEG sleep spindle propagation. 
     Part 1: Theory and proposed methodology, Journal of 
     Neuroscience Methods. doi: 10.1016/j.jneumeth.2013.08.013
 [2] O'Reilly, C., & Nielsen, T. (2013). The Spyndle toolbox: an open source 
     software package for the investigation of EEG transient events. 
     Revue internationale des technologies avancées, in press.     

"""

###############################################################################
# IMPORTS
###############################################################################
import sqlalchemy as sa
import numpy as np
from datetime import datetime

from spyndle.io import Base, rows2df





###############################################################################
# Creating the ORM database classes 
###############################################################################

class DataModelMng:
    """
    Class used to simplify the interaction with the data model. The code of this
    class has been separated from the DatabaseMng class to keep the latter class 
    independant of the data structure. 
    
    TODO: Complete this class to provide easy access to most data.
    """    

    def __init__(self, session):
        """
        DataModelMng constructor.
        
        Parameter:
            session : SQLAlchemy session used to interact with the data model.
        """    
        self.session = session
        

    def getTransientEvents(self, filteringDict={}, pandasFormat=False):
        """
        Return the transient events filtered using the filteringDict, either as
        a list of TransientEvent objects or a pandas.DataFrame object.
        """
        
        query = self.session.query(TransientEvent)
        
        for key in filteringDict:
            query = query.filter(getattr(TransientEvent, key) == filteringDict[key])

        if pandasFormat:
            return rows2df(query.all())
        
        return query.all()
    
        





###############################################################################
# Creating the ORM database classes 
###############################################################################


class TransientEvent(Base): 
    """
    Implement the ORM class for the representation of transient events such as 
    sleep spindles, K-complexes, slow waves, etc.
    
    TODO: The definition should depend on the type of transient event analyzed.
    """

    __tablename__       = "transientEvent"
    
    # The ID is a Universal Unique Identifier (UUID) that can be provided using
    # str(uuid.uuid1()) after the importation of the uuid module. 
    ID              = sa.Column(sa.String, primary_key=True)
    
    psgNight        = sa.Column(sa.String, sa.ForeignKey("psgNight.fileName"))
    channelName     = sa.Column(sa.String, sa.ForeignKey("channel.name"))
    eventName       = sa.Column(sa.String, sa.ForeignKey("eventClass.name"))
    startTime       = sa.Column(sa.Float)
    duration        = sa.Column(sa.Float)
    
    RMSamp          = sa.Column(sa.Float)
    meanFreq        = sa.Column(sa.Float)
    stage           = sa.Column(sa.String)
    cycle           = sa.Column(sa.Integer)
    slopeOrigin     = sa.Column(sa.Float)
    slope           = sa.Column(sa.Float)
    filteredRMSamp  = sa.Column(sa.Float)


    """
     Factory static method used to get a TransientEvent object from a
     spyndle.io.Event object and a psgNight string.
    """
    @staticmethod
    def fromEvent(event, psgNight):
        
        def getProperty(cls, label):
            return cls(event.properties[label]) if label in event.properties else None
                 
        startTime       = event.startTime
        duration        = event.duration()
        channelName     = event.channel
        eventName       = event.name
        
        RMSamp          = getProperty(float, "RMSamp")
        meanFreq        = getProperty(float, "meanFreq")
        stage           = getProperty(str,   "RMSamp")
        cycle           = getProperty(int,   "cycle")
        slopeOrigin     = getProperty(float, "slopeOrigin")
        slope           = getProperty(float, "slope")
        filteredRMSamp  = getProperty(float, "filteredRMSamp")   
 
        return TransientEvent(ID = event.ID, psgNight = psgNight, 
                              startTime = startTime, duration = duration, 
                              channelName = channelName, RMSamp = RMSamp, 
                              meanFreq = meanFreq, stage = stage, cycle = cycle, 
                              slopeOrigin = slopeOrigin, slope = slope, 
                              filteredRMSamp= filteredRMSamp, eventName=eventName)





class PSGNight(Base): 
    """
    Implement the ORM class for the representation of a polysomnographic (PSG) 
    recording night. For now, this class is essentially empty but it is implemented
    for future use (e.g., linking with subject record, recording montage, etc.).
    """    
    
    __tablename__       = "psgNight"
    
    fileName            = sa.Column(sa.String, primary_key=True)






class Channel(Base): 
    """
    Implement the ORM class for the representation of a recording channel. For now, 
    this class is essentially empty but it is implemented for future use 
    (e.g., recording gain, digitization frequency, etc.).
    """    
    __tablename__       = "channel"
    
    name                = sa.Column(sa.String, primary_key=True)
    #activeElectrode = sa.Column(sa.String)
    #reference       = sa.Column(sa.String)




     

class SPF(Base): 
    """
    Implement the ORM class for the representation of a spindle propagation field
    (SPF) recording channel (see [1] for the SPF definition). 
    """    
    __tablename__       = "spf"
    __table_args__      = {'sqlite_autoincrement': True}    
    
    no                  = sa.Column(sa.Integer, primary_key=True)
    psgNight            = sa.Column(sa.String, sa.ForeignKey("psgNight.fileName"))






class Propagation(Base): 
    """
    Implement the ORM class for the representation of the propagation of a 
    transient EEG event. 
    """
    
    __tablename__       = "propagation"
    __table_args__      = {'sqlite_autoincrement': True}

    no                  = sa.Column(sa.Integer, primary_key=True)
    
    # ID of the spindle from which this propagation has been computed.
    transientEventID    = sa.Column(sa.String, sa.ForeignKey("transientEvent.ID"))
    propRelNo           = sa.Column(sa.Integer, sa.ForeignKey("propagationRelationship.no"))
    
    sourceChannelName   = sa.Column(sa.String, sa.ForeignKey("channel.name"))
    sinkChannelName     = sa.Column(sa.String, sa.ForeignKey("channel.name"))
    
    similarity          = sa.Column(sa.Float)
    delay               = sa.Column(sa.Float)
    offset              = sa.Column(sa.Float)
    
    transientEvent      = sa.orm.relationship("TransientEvent", 
                                    primaryjoin='Propagation.transientEventID==TransientEvent.ID',
                                    lazy='joined')    
    
    
    # Used in SPF computation.
    noSPF               = sa.Column(sa.Integer, sa.ForeignKey("spf.no"))
    bidirect            = sa.Column(sa.Integer)
    source              = sa.Column(sa.Integer, sa.ForeignKey("propagation.no"))
    inverted            = sa.Column(sa.Boolean)    
    
    
    ######################
    # Rejection criteria  

    # Is the propagation is a false positive according to a threhold on the FDR?
    # (criterion C1)  
    isFP                = sa.Column(sa.Boolean)

    # Is the propagation an outlier according to Tuckey's test?  
    # (criterion C2)  
    isOutlier           = sa.Column(sa.Boolean)


    def __init__(self, *args, **kwargs):
        """
        Propagation constructor. 
        """
        super(Propagation, self).__init__(*args, **kwargs)  # call parent class init        
        
        # By default, propagations are considered unreliable. They must be proven
        # reliable.
        self.isFP       = True
        self.isOutlier  = True
        
        
        self.noSPF      = -1
        self.bidirect   = -1
        self.source     = -1
        self.inverted   = False
        
        if self.propRelNo is None:
            raise ValueError("The Propagation.propRelNo cannot be None.")




    def __repr__(self):
        """
        Propagation representative string. 
        """        
        return str(self.__dict__)







class EventClass(Base): 
    """
    Implement the ORM class for the representation of a class of events that have
    been annotated. They are specified by a name and the reference to a 
    DataManipulationProcess.
    """    
    
    __tablename__       = "eventClass"

    name                = sa.Column(sa.String, primary_key=True)
    dataManipNo         = sa.Column(sa.Integer, sa.ForeignKey("dataManipulationProcess.no"))





class DataManipulationProcess(Base): 
    """
    Implement the ORM class for the representation of a data manipulation process.
    Every manipulation of the raw data should be associated with a 
    DataManipulationProcess to provide the traceability information on how 
    computed data have been generated. For now, the information of the
    DataManipulationProcess are encoded in a reprString which can be the string
    representation of any data structure. This approach has been chosen because
    the information of a data manipulation process can vary a lot depending on the
    nature of the manipulation. 
    """    
    
    __tablename__       = "dataManipulationProcess"
    __table_args__      = {'sqlite_autoincrement': True}   
    
    no                  = sa.Column(sa.Integer, primary_key=True)
    reprStr             = sa.Column(sa.String)  
    datetime            = sa.Column(sa.DateTime)    
        
    
    #activeElectrode = sa.Column(sa.String)
    #reference       = sa.Column(sa.String)


    def __init__(self, reprStr):
        self.reprStr  = reprStr
        self.datetime = datetime.now()




class PropagationRelationship(Base): 
    """
    Implement the ORM class for the representation of a propagation relationshiop.
    Propagation relationships represents the average propagation behavior between
    a pair of source-sink channel for a given class of events on a given recoding
    night.
    """    
    
    __tablename__       = "propagationRelationship"
    __table_args__      = (
                            sa.UniqueConstraint("sourceChannelName", 
                                               "sinkChannelName",
                                               "psgNight",
                                               "eventName"), 
                            {'sqlite_autoincrement': True} 
                          )

    no                  = sa.Column(sa.Integer, primary_key=True)
    sourceChannelName   = sa.Column(sa.String, sa.ForeignKey("channel.name"))
    sinkChannelName     = sa.Column(sa.String, sa.ForeignKey("channel.name"))
    psgNight            = sa.Column(sa.String, sa.ForeignKey("psgNight.fileName"))
    eventName           = sa.Column(sa.String, sa.ForeignKey("eventClass.name"))
    
    cutoff              = sa.Column(sa.Float)
    FDR                 = sa.Column(sa.Float)
    nbValid             = sa.Column(sa.Integer)
    nbOut               = sa.Column(sa.Integer)
    delay_mean          = sa.Column(sa.Float)
    delay_sd            = sa.Column(sa.Float)
    
    isValidC3           = sa.Column(sa.Boolean)
    isValidC4           = sa.Column(sa.Boolean)
             
    

    def testRejectionC3(self, deltaWindow = 0.5, alphaSD = 0.2):
        """
        Method testing if the propagation relationship passes the C3 criterion 
        [1,2]. Stores the result in the isValidC3 attribute.
        """        
        # Applying rejection criterion c3 and saving the result
        sdThreshold = deltaWindow*alphaSD/np.sqrt(12)  
        self.isValidC3 = self.delay_sd < sdThreshold

           

    def testRejectionC4(self, minNbValid=40):
        """
        Method testing if the propagation relationship passes the C4 criterion 
        [1,2]. Stores the result in the isValidC4 attribute.
        """        
        # Applying rejection criterion c3 and saving the result
        self.isValidC4 = self.nbValid >= minNbValid
        



    def getValidPropagations(self, session):
        """
        Return the set of valid propagation (i.e. passing the C3 and C4 criteria)
        that are associated with the propagation relationship.
        """        
        return session.query(Propagation).filter_by(propRelNo  = self.no,
                                                    isFP       = False,
                                                    isOutlier  = False).all()
                    


    def add(self, session, behavior = "raise"):
        """
        This methods allow adding the propagation relationship to the database
        associated with session. It allows managing the different situations that
        can arise on adding the record such as integrity errors. Possible behaviors
        are:
            - raise          : simply raise the IntegrityError exeception.
            - updateSilently : if a conflicting record already exists in the 
                               database, replace it with this new record, taking
                               the propagation relationship numero of the old and
                               now deleted record. Do not raise any warning or error.
            - failSilently   : Fail to add the new record but do not raise any 
                               error or warning. In short, on error, do nothing.                            
        """
        
        try:
            session.add(self)
            session.commit()
            
        except sa.exc.IntegrityError:
            
            if behavior == "raise":
                raise
             
            session.rollback()
            if behavior == "updateSilently":
                
                old = session.query(PropagationRelationship)\
                            .filter_by(sourceChannelName = self.sourceChannelName,
                                       sinkChannelName   = self.sinkChannelName,
                                       psgNight          = self.psgNight,
                                       eventName         = self.eventName ).one()
                self.no = old.no
                session.delete(old)
                session.add(self)
                session.flush()

            elif behavior == "failSilently":
                pass
            
            else:
                raise ValueError("The value '" + behavior + "' is invalid for the behavior variable.")




    """
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
                   
    """
    """          
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
    """






"""        
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
"""