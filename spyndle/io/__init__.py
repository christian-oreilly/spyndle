# __init__.py
# -*- coding: utf-8 -*-

"""
    Initialization of the spyndle.io module.
    
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
import sqlalchemy.ext.declarative as sad

# Getting an SQLAchemy base class
Base    = sad.declarative_base()

# Getting an SQLAchemy session factory
Session = sa.orm.scoped_session(sa.orm.sessionmaker()) #sa.orm.scoped_session(sa.orm.sessionmaker())#expire_on_commit=False))


# Reader base
from spyndle.io.EEGDatabaseReader import Event, RecordedChannel, \
    EEGDBReaderBase, EventList
    
# Readers
from spyndle.io.edf import EDFReader, EDFMultiReader
from spyndle.io.harmonie import HarmonieReader
from spyndle.io.devuyst import DevuystReader

# Database management
from spyndle.io.dbutils import rows2df
from spyndle.io.dataModel import DataModelMng, TransientEvent, PSGNight, \
    Channel, SPF, Propagation, EventClass, DataManipulationProcess, \
    PropagationRelationship, SlowWaveEvent, SpindleEvent
from spyndle.io.databaseMng import DatabaseMng, clearDatabase

    
    
