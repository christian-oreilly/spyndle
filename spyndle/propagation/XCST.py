# -*- coding: utf-8 -*-

"""
    Code assessing transient event propogation through an array of sensors
    using cross-correlation of S-transform of the signal captured by
    the different sensors.

    Copyright (C) 2012-2013  Christian O'Reilly

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
 Date  : July 11, 2012

 [1] O’Reilly, C., & Nielsen, T. (2013). Assessing the propagation of EEG 
     transient activity, Proceedings of the 8th International Workshop on 
     Systems, Signal Processing and their Applications, Algiers, Algeria, 
     12-15 May 2013, 132-139.
 [2] O'Reilly, C. & Nielsen, T. (2013) Assessing EEG sleep spindle propagation. 
     Part 1: Theory and proposed methodology, Submitted to Journal of 
     Neuroscience Methods. doi: 10.1016/j.jneumeth.2013.08.013
 [3] O'Reilly, C. & Nielsen, T. (2013) Assessing EEG sleep spindle propagation. 
     Part 2: Experimental characterization, Submitted to Journal of 
     Neuroscience Methods. doi: 10.1016/j.jneumeth.2013.08.014    
 [4] Devuyst, S. (2013). The DREAMS Sleep Spindles Database, 2013, from 
     http://www.tcts.fpms.ac.be/~devuyst/Databases/DatabaseSpindles/.
"""


import os
from scipy import where, zeros, transpose, arange, shape, diff, sign, argmax, sqrt
import numpy as np
import sqlalchemy as sa
import serpent
from datetime import datetime

from spyndle import diff2, __version__
from spyndle import computeMST
from spyndle.io import Propagation, TransientEvent, \
    DatabaseMng, PropagationRelationship, DataManipulationProcess
from spyndle.EEG.system10_20 import getEEGChannels



def computeXCST(readerClass, fileName, eventName, dbSession=None,  dbPath=None,  shard=None,
                verbose=False, **kwargs):
    """
    Helper function for a one-function-call processing of the 2D cross-corellation 
    of S-Transform spectra described in [1, 2].
    
    Parameters:
        readerClass    : Class of EEG reader. As of now, the HarmonieReader 
                         (.sig/.sts), the EDFReader (.edf, .bdf), and the 
                         DevuystReader (data from the DREAMS database [4]) 
                         classes have been implemented. Users can implement their 
                         own reader for EEG file format they are using folowwing 
                         the example given by the existing reader classes.
        fileName       : File name of the EEG data.
        eventName      : Name used to label the transient event for which the 
                         propagation is XCST is to be computed
        dbSession      : [optional, default=None] A SQLAlchemy session used to 
                         read/write from the database. If equal to None, an 
                         in-memory sqlite database is used.
        dbPath         : [optional, default=None] Database path that will be 
                         used if no session are provided.                
        verbose        : [optional, default=False] If set to True, information 
                         messages about the execution will be printed.
        **kwargs       : Keyword arguments that will be passed to the 
                         spyndle.propagation.XCSTEvaluator.compute(...) method.
                         
    Note: This function returns nothing. Computation results are stored in an
    SQL database which can be chosen by setting the dbSession parameter or 
    that is an in-memory sqlite database.
    """                  
    evaluator = XCSTEvaluator(eventName, verbose)        
    evaluator.createEEGReader(fileName, readerClass)
    evaluator.prepareDatabase(dbSession, dbPath, shard)   
    evaluator.compute(**kwargs)
            
            
        



class PropagationEvaluator:
    """
    Base class for subclasses providing evaluator objects that can be used to 
    compute transient event propagation
    
    Default values of the parameters of the algorithm are specified in the 
    constructor. These are as follow: 
    
        delta [default = 0.5 (in seconds)] : Correspond to the delta_w 
        in [1]. The observation windows will be slid from -delta to +delta.
        
        beforePad [default = 0.5 (in seconds)] : Padding of the observation 
        windows before the beginning  of the spindle. Labeled delta_{sp} in [1].
        
        afterPad [default = 0.5 (in seconds)] : Padding of the observation 
        windows after the ending of the spindle. Labeled delta_{sp} in [1].
        
        artifactPad [default = 0.25 (in seconds)] : Used to pad the ST to 
        compute it on a larger window that is afterward trucated such as to 
        eliminate border artifacts.
        
        offset [default = 0.0 (in seconds)] : Should be equal to 0.0 for 
        computing the delays. A value of 5.0 can be used to compute false 
        positive propagation. See [1].
        
        similIndexType  = "inf", 
        
        intFct [default = numpy.sum] : Allows to choose the function used in 
        integration. By default numpy.sum is used but numpy.trapz could be used
        for more precision (at cost of a much slower execution).
        
    Different values for these parameters can be specified at computation time
    by specifying keyword arguments to the compute(...) method. 
    
    Other class attributes: 
    
        __spyndle_version : Version of the spyndle package. It is important
        that it is stored whitin the object for traceability when the __repr__
        of the object is associated with saved results.
        
        propRelNos : Filled during database preparation with a dictionnary 
        linking SLQ propagation relationship primary key. This is use for an
        easy access to these numbers for specifying the propagation foreign key
        at saving time. Keys are concatenation of reference and test channel 
        names.
        
        channelList: List of channels for which we want to compute the propagation.
        
        reader : EEG reader object used for assessing the EEG data.
        
        dbSession : SQLAlchemy session object for interfacing with the database.
        
        dbMng : DatabaseMng object used for opening and closing an SQL session 
        if no session are provided by external call. If a session is provided 
        by an external call, this attribute is left to a None value.
        
        dataManipObj : DataManipulationProcess ORM SQLAlchemy object used to 
        register the data manipulation process performed by this object to 
        provide traceability of data processing.
    
    """

    """
    Temporary patch because using simply np.sum or np.trapz for the value of 
    the intFct attribute does not work for an unknown reason.
    """    
    @staticmethod
    def intSum(X, axis=0): 
        return np.sum(X, axis=0)
        
    @staticmethod        
    def intTrapz(X, axis=0) : 
        return np.trapz(X, axis=0)
  

    
    """
     A flag specifying wheter an exception should be raised when an event is 
     encountered outside the recorded epoch. Note that such an "outside event" can
     also occurs if the offset applied to the XCST is such that an "inside event"
     is pushed outside of the recorded epoch by the offset. If the flag is set 
     to False, the "outside events" are simply dropped such that no propagation
     are computed for them.
    """
    __RAISE_ON_OUTSIDE_EVENT__ = False
  

    def __init__(self, eventName, verbose = False):
        """
        Constructor for the XCSTEvaluator objects.
    
        Parameters:
            eventName      : Name used to label the transient event for which the 
                             propagation is XCST is to be computed
            verbose        : [optional, default=False] If set to True, information 
                             messages about the execution will be printed.
        """                            
        self.eventName          = eventName    
        self.verbose            = verbose
    
        ######################################################################
        # Setting default values for parameters of the algorithm
        self.delta              = 0.5 
        self.beforePad          = 0.5 
        self.afterPad           = 0.5 
        self.artifactPad        = 0.25 
        self.offset             = 0.0 
        self.similIndexType     = "inf" 
        self.intFct             = XCSTEvaluator.intSum    
      
        # Setting default values for other attributes.
        self.__spyndle_version  = __version__       
        self.propRelNos         = {}
        self.channelList        = []     
        self.reader             = None
        self.dbSession          = None
        self.dbMng              = None
        self.dataManipObj       = None    




        
    def __repr(self):        
        """
        Use the serpent serialization algorithm to provide a string 
        representation of the object.
        """
        return serpent.dumps(self, indent=False, set_literals=False)        
        

    def createEEGReader(self, fileName, readerClass):
        """
        Create and EEG reader object to be used for accessing EEG data.
        Alternatively, the reader can be set using the reader property.
        
        Parameters:
            fileName    : Full name of that file containing the EEG data. The 
                          type of the file should agree with the readerClass.
            readerClass : Class of EEG reader. Should be one of the reader 
                          classes (e.g. EDFReader, HarmonieReader).
        """        
            
        try:
            if self.verbose:
                print('Loading data of the file ', end=' ')
                print(os.path.basename(fileName) + '...')
                
            self.reader = readerClass(fileName)    
            
        except IOError as e:     
            errMsg = "Error : computeXCST : " + fileName + "could not be opened.\n"
            IOError(errMsg + str(e))
            
        self.channelList = getEEGChannels(self.reader.getChannelLabels())            
   
   
   
   
    @property
    def reader(self):
        """
        Property to set or get the EEG reader object used for assessing the EEG
        data.
        """
        return self.__reader
    
    @reader.setter
    def reader(self, reader):
        self.__reader = reader
   



    @property
    def dbSession(self):
        """
        Property to set or get the SQLAlchemy session object for interfacing 
        with the database.
        """        
        return self.__dbSession
    
    @dbSession.setter
    def dbSession(self, dbSession):
        self.__dbSession = dbSession



    @property
    def fileName(self):
        """
        Property for accessing the name of the file containing the EEG data. 
        Only a getter is provided. Setting a new filename should be done by
        interacting directly the the EEG reader. Caution is advisable when
        changing the fileName to avoid breaking traceability of data 
        manipulations.
        """                
        return self.reader.fileName




    def prepareDatabase(self, dbSession=None, dbPath=None, shard=None): 
        """
        Set-up the SQL database for recording results from data processing. 
        Create the PropagationRelationship records corresponding to the channel
        contained in the EEG file, with an association to this specific file and
        the object eventName atribute.
        
        Parameters:
            dbSession : [optional, default=None] SQLAlchemy session object. 
                        If left to None, an in-memory sqlite database is used.      
        """    
        
        if dbSession is None:
            try:
                if dbPath is None:
                    if self.verbose:
                        print("No database session or path has been passed. "\
                              "Using in-memory database.")
                    dbPath = "sqlite://"
                        
                self.dbMng      = DatabaseMng(dbPath, shard=shard)
                self.dbSession  = self.dbMng.session

            except sa.exc.ArgumentError as error:
                errMsg = "Error connecting to the specified database URL. "\
                         "The format of the database path is invalid."     \
                         "\n\nSQLAlchemy error message: " + str(error)
                raise sa.exc.ArgumentError(errMsg)
        else:
            self.dbSession = dbSession
    
        # Create a new data manipulation process to keep track of 
        # the execution of this code.   
        self.dataManipObj = DataManipulationProcess(reprStr  = repr(self))
    
    
    
    
    def createPropRel(self):
        # Create the needed  propagation relationship objects and keep note of
        # their id.
        for ref in self.channelList:
            for test in self.channelList:
                if ref != test:            
                    propRel = PropagationRelationship(sinkChannelName    = test, 
                                                      sourceChannelName  = ref, 
                                                      psgNight           = self.fileName,
                                                      eventName          = self.eventName)
                                                      
                    # Add the propagation relationship if it does not already
                    # exist. If it already exist, keep its id number but update it
                    # with the new PropagationRelationship object.
                    propRel.add(self.dbSession, behavior = "updateSilently")
                    self.dbSession.flush()
                    self.propRelNos[ref + test] = propRel.no
    
        

    

    def compute(self, **kwargs):
        """
        Perform the 2D cross-corellation of S-Transform spectra described in 
        [1, 2] and save the results in the connected database.
        
        Parameters:
            **kwargs : [optional] Any keyword arguments corresponding to the  
                       classes attribute (see the module definition for the 
                       list of these attributes along with their default values).     
        """    

        print("compute...")

        for key in kwargs:
            if key in self.__dict__:
                setattr(self, key, kwargs[key])
            else:
                ValueError("The keyword " + str(key) + " is unknown.") 

        # Create the propagation relationship related to the selected channel
        # for this computation. 
        self.createPropRel()
  
        # Using the LIMIT approach has cons (see http://stackoverflow.com/
        # questions/7389759/memory-efficient-built-in-sqlalchemy-iterator-generator)
        # However, using the windows approach proposed in that link does not
        # works with all DB engine, include MySQL.  
  
  
        eventQuery = self.dbSession.query(TransientEvent)\
                        .filter(TransientEvent.channelName.in_(self.channelList),
                                TransientEvent.psgNight  == self.fileName,
                                TransientEvent.eventName == self.eventName)
  
        if self.verbose:
            print("Processing propagation for " + str(eventQuery.count()), end=' ')
            print(" transient events " + self.eventName + " present in file ", end=' ')
            print(os.path.basename(self.fileName) + " on channels ", end=' ')
            print(str(self.channelList))   
        
        #if self.verbose and len(events) > 100:
        #    
        #    t1 = datetime.now()     
        #    for event in events[:100]:        
        #        self._computeForEvent(event)   
        #    t2 = datetime.now()    
        #    
        #    print "One hundred events have been processed in " + str((t2-t1).total_seconds()),
        #    print " seconds. Estimated total processing time : ", 
        #    print str((t2-t1).total_seconds()/100.0*len(events)) + " seconds."
        #        
        #    for event in events[100:]:        
        #        self._computeForEvent(event)   
        #         
        #else:
  
  
        offset = 0
        while True:
            limitedQuery = eventQuery.offset(offset).limit(10000)
            if limitedQuery.count() == 0: 
                break
            for event in limitedQuery:
                # Doing this instead of "for event in eventQuery.all():" makes that the
                # ORM processing is made one event at a time instead of for all events
                # at once.                
                offset += limitedQuery.count()
                self._computeForEvent(event)              
            
        # Update the data manipulation process object.            
        self.dataManipObj.reprStr = repr(self)



            
            
    def getWindowedSignal(self, event, deltaSampleStart, deltaSampleEnd, 
                                  channelList, offset):

        refChannel  = event.channelName
        refFreq     = self.reader.getChannelFreq(refChannel)
        startTime   = event.startTime -  self.beforePad  - self.artifactPad + offset 
        duration    = event.duration + self.beforePad + self.afterPad + 2*self.artifactPad

        deltaStart       = deltaSampleStart/refFreq    
        deltaEnd         = deltaStart

        """
         The signal at startTime - delta might not be available and the 
         reader will start reading at the next available sample following 
         the startTime requested. In some cases, it can happen that events
         have been specified outside the time epoch covered by the recording 
         at hand. Given on the value of self.__RAISE_ON_OUTSIDE_EVENT__, this
         can be considered as an error or not.
        """
        refStartTime = self.reader.getNextSampleStartTime(refChannel, startTime - deltaStart)
        if refStartTime is None:
            if self.__RAISE_ON_OUTSIDE_EVENT__:
                raise ValueError("The event is not occurring during the recorded epoch.")
            else:
                return None, None, None        
        
        """         
         If startTime - delta was not occuring during recording period, the next
         available sample was taken as the beginning of the signal. When this 
         happens, we may miss some samples in this signal. We remove these 
         samples in both synchronous and asynchronous signals (for them to have 
         the same size) by reducing the delta accordingly.
        """
        deltaStartRef = int((refStartTime - (startTime - deltaStart ))*refFreq)
            
        assert(deltaStartRef >=0)
        deltaSampleStart -= deltaStartRef
        deltaStart        = deltaSampleStart/refFreq                
        
        refStart      = startTime - deltaStart
        refDuration   = duration + deltaStart + deltaEnd
        signalsDataRef = self.reader.read(channelList, refStart, refDuration)  
        
        """
         It can happen that the signal is not available for the whole refDuration
         period. In this case, we want to remove the unaivalable portion from
         the deltaEnd padding.
        """
        deltaEndRef = int(refDuration*refFreq) - len(signalsDataRef[refChannel].signal)

        assert(deltaEndRef >=-1)    
        
        deltaSampleEnd -= deltaEndRef       
            
            
        return signalsDataRef, deltaSampleStart, deltaSampleEnd
            
         
         
         
         
    def getSignals(self, event):

        refChannel        = event.channelName
        fs                = self.reader.getChannelFreq(refChannel)        
    
        deltaSampleStart  = int(self.delta*fs)          
        deltaSampleEnd    = deltaSampleStart    

        signalsDataRef, deltaSampleStartRef, deltaSampleEndRef = \
            self.getWindowedSignal(event, deltaSampleStart, deltaSampleEnd, 
                                   [refChannel], 0)  
        signalsDataCmp, deltaSampleStartCmp, deltaSampleEndCmp = \
            self.getWindowedSignal(event, deltaSampleStart, deltaSampleEnd, 
                                   self.channelList,  self.offset)  
        
        if signalsDataRef is None or signalsDataCmp is None:
            return 
        
        # In case  the reference or the test signal are not sufficiently long 
        # signal before and after the spindle to compute delays that are at 
        # least half of the limit asked by the parameter delta, we reject 
        # the spindle.
        if deltaSampleStartRef < deltaSampleStart/2  or \
           deltaSampleEndRef   < deltaSampleEnd/2    or \
           deltaSampleStartCmp < deltaSampleStart/2  or \
           deltaSampleEndCmp   < deltaSampleEnd/2    :      
            return


        
        #################################################################
        # We keep the intersection of the time duration of the two signals.
        dtStart = deltaSampleStartRef - deltaSampleStartCmp
        if dtStart < 0 :
            for channel in signalsDataCmp:
                signalsDataCmp[channel].signal = signalsDataCmp[channel].signal[-dtStart:]
                deltaSampleStart = deltaSampleStartRef                      
        elif dtStart > 0:
            signalsDataRef[refChannel].signal = signalsDataRef[refChannel].signal[dtStart:]
            deltaSampleStart = deltaSampleStartCmp                         
                
        dtEnd = deltaSampleEndRef - deltaSampleEndCmp
        if dtEnd < 0 :
            for channel in signalsDataCmp:
                signalsDataCmp[channel].signal = signalsDataCmp[channel].signal[:dtEnd]        
                deltaSampleEnd   = deltaSampleEndRef   
        elif dtEnd > 0:
            signalsDataRef[refChannel].signal = signalsDataRef[refChannel].signal[:-dtEnd]        
            deltaSampleEnd   = deltaSampleEndCmp   
                
        assert(len(signalsDataCmp[refChannel].signal) == len(signalsDataRef[refChannel].signal))
         
        return signalsDataRef, signalsDataCmp, deltaSampleStart, deltaSampleEnd


         
         
         
    def _computeForEvent(self, event):       
        """
        """  
        raise NotImplementedError("This method needs to be implemented in subclasses of PropagationEvaluator.")        
       






















class XCSTEvaluator(PropagationEvaluator):
    """
    Class providing evaluator objects that can be used to compute
    the 2D cross-corellation of S-Transform spectra described in [1, 2]. For 
    performing this computation with a one-liner, please refer to the
    spyndle.propagation.computeXCST(...) function.
    """

         
    def _computeForEvent(self, event):       
        """
        Perform the 2D cross-corellation of S-Transform spectra for a specific
        event.
        
        Parameters:
            event : Event for which the XCST will be computed.     
        """  

        refChannel        = event.channelName
        fs                = self.reader.getChannelFreq(refChannel)                
        artifactPadSample = int(self.artifactPad*fs)                        
        signalsDataRef, signalsDataCmp, deltaSampleStart, deltaSampleEnd = self.getSignals(event)        
        
        # Spectra computation
        X, fXRef = computeMST(signalsDataRef[refChannel].signal, fs, m=0.0, k=1.0, fmin=11.0, fmax=16.0)     
        X = abs(transpose(X))                
        try:
            #Y = zeros(X.shape)
            #Y[:, (deltaSampleStart+artifactPadSample):-(deltaSampleEnd+artifactPadSample)] = X[:, (deltaSampleStart+artifactPadSample):-(deltaSampleEnd+artifactPadSample)]
            X = {"Ref" :X[:, (deltaSampleStart+artifactPadSample):-(deltaSampleEnd+artifactPadSample)]}
            #X["Ref2"] = Y
        except IndexError:
            print(deltaSampleStart, deltaSampleEnd, X.shape, artifactPadSample)
            raise               
        
        for testChannel in self.channelList :      
            Xtmp, fXtmp = computeMST(signalsDataCmp[testChannel].signal, fs, m=0.0, k=1.0, fmin=11.0, fmax=16.0)         
            X[testChannel] = abs(transpose(Xtmp))[:, artifactPadSample:-artifactPadSample]              
        
        refShape = list(shape(X["Ref"]))      
        
        
        unnormCrossCorr = zeros(deltaSampleStart + deltaSampleEnd)
        time        = (arange(deltaSampleStart + deltaSampleEnd)-deltaSampleStart)/fs
        refSelfCor  = self.intFct(self.intFct(X["Ref"]*X["Ref"]))
        indOffsets = arange(deltaSampleStart + deltaSampleEnd)
        
        propagations = []
        for testChannel in self.channelList :
            if testChannel == refChannel:  
                continue                
                
            XCmpSquare = X[testChannel]*X[testChannel]
            XCmpSquare_colTrapz = self.intFct(XCmpSquare, axis=0)     
            
            
            """
            equivalent to:
            CmpSelfCor  = zeros(deltaSampleStart + deltaSampleEnd)                        
            for indOffset in indOffsets:
                CmpSelfCor[indOffset]       = sum(XCmpSquare_colTrapz[indOffset:(indOffset+refShape[1])])                         
            """
            cums = np.cumsum(XCmpSquare_colTrapz)
            CmpSelfCor = np.concatenate(( [cums[refShape[1]-1]], 
                                        cums[indOffsets[1:] + refShape[1]-1] - cums[indOffsets[:-1]]))
    
    
                              
            # This...                        
            #unnormCrossCorr = correlate2d(X[testChannel], X["Ref"], mode="valid")[0][indOffsets]   
    
            # and this...
            #K = X[testChannel].shape[1] 
            #H = X["Ref"].shape[0]
            #N = X["Ref"].shape[1]   
            #y = np.concatenate((X["Ref"], zeros((H, K-1))), axis=1)           
            #FFTy = np.conjugate(np.fft.fft(y, axis=1)) 
            #x = np.concatenate((zeros((H, N-1)), X[testChannel]), axis=1)                        
            #unnormCrossCorr = np.real(np.fft.ifft(np.fft.fft(x, axis=1)*FFTy, axis=1)[:, (N-1):(N-1 + deltaSampleStart + deltaSampleEnd)]) 
            #unnormCrossCorr = np.sum(unnormCrossCorr, axis=0)                        
    
    
            # are slower than that:
            for indOffset in indOffsets:
                XCmp                        = X[testChannel][:, indOffset:(indOffset+refShape[1])]
                unnormCrossCorr[indOffset]  = self.intFct(self.intFct(XCmp*X["Ref"]))            
    
    
            if self.similIndexType == "inf":
                den = np.maximum(CmpSelfCor, refSelfCor)
            elif self.similIndexType == "euclidean":
                den = sqrt(CmpSelfCor**2 + refSelfCor**2)                                
            else: 
                print("Error. The similarity index type " + self.similIndexType + "is unknown.")
                raise TypeError
                
            crosscor = unnormCrossCorr/den
            assert(np.all((crosscor >= 0.0)*(crosscor <= 1.0)))    
    
            ind          = argmax(crosscor)                    
            diffCrosscor = diff2(arange(deltaSampleStart + deltaSampleEnd), crosscor)                
    
            # Interpolation to find the position of the true maximum
            # which normally falls between samples
            maxInd = len(crosscor)-1
            if ind == 0 or ind == maxInd:
                maxDeltay   = time[ind]       
                maxCrosscor = crosscor[ind]                      
            else:
                if ind-2 > 0 :
                    indm = ind-2
                else:
                    indm = 0
                    
                if ind+2 < maxInd :
                    indp = ind+2
                else:
                    indp = maxInd
    
                pindinc = where(diff(sign(diffCrosscor[indm:(indp+1)])) != 0)[0][0]                
                indFloor = indm +pindinc 
                indCeil  = indFloor + 1
                indFrac  = diffCrosscor[indFloor]/(diffCrosscor[indFloor] - diffCrosscor[indCeil])
              
                maxDeltay   = time[indFloor]      + indFrac*(time[indCeil]-time[indFloor])       
                maxCrosscor = crosscor[indFloor]  + indFrac*(crosscor[indCeil]-crosscor[indFloor])       
                                
            assert(abs(maxDeltay) <= self.delta)
            
            propagations.append(Propagation(transientEventID=event.ID, sinkChannelName=testChannel, 
                                    sourceChannelName=refChannel, similarity=maxCrosscor,
                                     delay=maxDeltay, offset=self.offset, 
                                    propRelNo = self.propRelNos[refChannel + testChannel]))
        self.dbSession.add_all(propagations)
        self.dbSession.flush()
        
        
        

















        



class TDPropagationEvaluator(PropagationEvaluator):
    """
    Class providing evaluator objects that can be used to compute
    the best time alignement the time-domain signals using a cross-corellation.     
    """

       
    def _computeForEvent(self, event):       

        refChannel        = event.channelName
        fs                = self.reader.getChannelFreq(refChannel)                                 
        signalsDataRef, signalsDataCmp, deltaSampleStart, deltaSampleEnd = self.getSignals(event)        
        artifactPadSample = int(self.artifactPad*fs)                              
                      
        X = {"Ref" :signalsDataRef[refChannel].signal[(deltaSampleStart+artifactPadSample):-(deltaSampleEnd+artifactPadSample)]}
        for testChannel in self.channelList :      
            X[testChannel] = signalsDataCmp[testChannel].signal[artifactPadSample:-artifactPadSample]              
        
        sigLen = len(X["Ref"])      
        
        unnormCrossCorr = zeros(deltaSampleStart + deltaSampleEnd)
        time        = (arange(deltaSampleStart + deltaSampleEnd)-deltaSampleStart)/fs
        refSelfCor  = self.intFct(X["Ref"]*X["Ref"])
        indOffsets = arange(deltaSampleStart + deltaSampleEnd)
        
        propagations = []
        for testChannel in self.channelList :
            if testChannel == refChannel:  
                continue                
                
            XCmpSquare = X[testChannel]*X[testChannel]
            
            
            cums = np.cumsum(XCmpSquare)
            CmpSelfCor = np.concatenate(([cums[sigLen-1]], 
                                        cums[indOffsets[1:] + sigLen-1] - cums[indOffsets[:-1]]))
    
            # are slower than that:
            for indOffset in indOffsets:
                XCmp                        = X[testChannel][indOffset:(indOffset+sigLen)]
                unnormCrossCorr[indOffset]  = self.intFct(abs(XCmp)*abs(X["Ref"])) # This index can only be computer on positive real functions, thus the abs(...)            
    
    
            if self.similIndexType == "inf":
                den = np.maximum(CmpSelfCor, refSelfCor)
            elif self.similIndexType == "euclidean":
                den = sqrt(CmpSelfCor**2 + refSelfCor**2)                                
            else: 
                print("Error. The similarity index type " + self.similIndexType + "is unknown.")
                raise TypeError
                
            crosscor = unnormCrossCorr/den
            assert(np.all((crosscor >= 0.0)*(crosscor <= 1.0)))
    
            ind          = argmax(crosscor)                    
            diffCrosscor = diff2(arange(deltaSampleStart + deltaSampleEnd), crosscor)                
    
            # Interpolation to find the position of the true maximum
            # which normally falls between samples
            maxInd = len(crosscor)-1
            if ind == 0 or ind == maxInd:
                maxDeltay   = time[ind]       
                maxCrosscor = crosscor[ind]                      
            else:
                if ind-2 > 0 :
                    indm = ind-2
                else:
                    indm = 0
                    
                if ind+2 < maxInd :
                    indp = ind+2
                else:
                    indp = maxInd
    
                pindinc = where(diff(sign(diffCrosscor[indm:(indp+1)])) != 0)[0][0]                
                indFloor = indm +pindinc 
                indCeil  = indFloor + 1
                indFrac  = diffCrosscor[indFloor]/(diffCrosscor[indFloor] - diffCrosscor[indCeil])
              
                maxDeltay   = time[indFloor]      + indFrac*(time[indCeil]-time[indFloor])       
                maxCrosscor = crosscor[indFloor]  + indFrac*(crosscor[indCeil]-crosscor[indFloor])       
                                
            assert(abs(maxDeltay) <= self.delta)
            
            propagations.append(Propagation(transientEventID=event.ID, sinkChannelName=testChannel, 
                                    sourceChannelName=refChannel, similarity=maxCrosscor,
                                     delay=maxDeltay, offset=self.offset, 
                                    propRelNo = self.propRelNos[refChannel + testChannel]))
        self.dbSession.add_all(propagations)
        self.dbSession.flush()
        
        
        

