# -*- coding: utf-8 -*-

"""
    Module implementing the computation of the spindle propagation files (SPF)
    as defined in [1, 2].

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
    If used for research purpose, the reference [1, 2] should  be cited in the 
    derived publication to refere the reader to the description of the 
    methodology.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


 Author: Christian O'Reilly (christian.oreilly@umontreal.ca)
 Date  : Jun 10, 2013

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
"""


###############################################################################
# IMPORTS
###############################################################################
from scipy import unique, where, percentile, concatenate, array
from scipy import sqrt, nan, median, std
from sklearn.covariance import MinCovDet
from numpy import in1d
import sqlalchemy
from sqlalchemy import distinct, and_

from spyndle.propagation import computeXCST
from spyndle.io import Propagation, TransientEvent, \
    DatabaseMng, PropagationRelationship, Channel, rows2df



###############################################################################
# Internal functions
###############################################################################


"""
  Compute outliers according to the nonparametric test proposed by Tukey (1977).
  
  Return the indexes of the outliers. 
  
  Tukey JW. Exploratory data analysis. Addison-Wesley Pub. Co.: Reading, Mass. ; 
  Don Mills, Ont., 1977.
"""
def getOutliers(data, coef=1.5):
    t1, t2 = getOutlierThresholds(data, coef=1.5)
    return where((data > t2)*(data < t1))[0]

def getInliers(data, coef=1.5):
    t1, t2 = getOutlierThresholds(data, coef=1.5)
    return where((data <= t2)*(data >= t1))[0]

       
def getOutlierThresholds(data, coef=1.5):
    Q1   = percentile(data, 25.0)
    Q3   = percentile(data, 75.0)
    return Q1 - coef*(Q3-Q1), Q3 + coef*(Q3-Q1), 




class SPFEvaluator : 
    
    def __init__(self, night, eventName, dbSession=None, dbName = "sqlite:///:memory:", verbose = False):
        
        self.night      = night
        self.eventName  = eventName
        self.verbose    = verbose
        
        if not dbSession is None:
            self.session  = dbSession
            self.dbMng      = None
            
        else:
            if dbName == "sqlite:///:memory:":
                print "No database session passed. Using in-memory database."
            else:     
                print "No database session passed. Using the " + dbName + " database."           
            try:

                self.dbMng = DatabaseMng(dbName)
                self.session = self.dbMng.session
            except sqlalchemy.exc.ArgumentError, error:
                 print "Error connecting to the specified database URL. "\
                       "The format of the database path is invalid."     \
                       "\n\nSQLAlchemy error message: " + str(error)
            
        


    def __del__(self):
        if not self.dbMng is None:
            self.dbMng.disconnectDatabase()
    
    
    

    """
      Wrapper around the XCST.computeXCST(...) functions to simplify the computation
      of XCST when performed in for SPF computation.
    """
    def computeXCST(self, readerClass, **kwargs) :  
        computeXCST(readerClass, self.night, self.eventName, 
                    dbSession=self.session, verbose=self.verbose, **kwargs) 



    """
     Computing SPF for self.night.
    """
    def computeSPF(self, channelList = None, addBehavior='raise', alpha = 10.0):        
                        
        self.addBehavior = addBehavior
                
        if channelList is None:
            # Getting all the records for the night
            channelList = self.session.query(distinct(Channel.name))\
                            .join(TransientEvent, Channel.name == TransientEvent.channelName)\
                            .filter(and_(TransientEvent.psgNight == self.night,
                                         TransientEvent.eventName == self.eventName)).all()
            if len(channelList):
                channelList = list(zip(*channelList)[0])
            else:
                raise UserWarning("No channel availables for computing SPFs.")
    

        # Computing cutoff threshold from asychronized comparisons
        if self.verbose :
            print "Computing cutoff for night ", self.night
        self.computingCutOff(channelList, alpha)
    
        # FALSE DETECTION AND OUTLIERS REJECTION
        if self.verbose:
            print "Rejecting outliers"
        self.propagationRejection(channelList, alpha)
    
        # CORRECTING DELAYS TO HAVE ONLY PROPAGATION WITH POSITIVE TIME DELAYS 
        if self.verbose: 
            print "Correcting negative delays"   
        self.negativeDelayCorrection()    

        # COMPUTING SPINDLE PROPAGATION FIELD      
        if self.verbose: 
            print "Identifying SPF" 
        self.identifySPF()            
    
        # REMOVING DUPLICATE PROPAGATION DUE TO BIDIRECTIONNALITY
        if self.verbose: 
            print "Removing duplicates"        
        self.removeBidirectionnalDuplicates()  
                
    

    
    ###########################################################################    
    # Computing the cutoff
    ###########################################################################
    def computingCutOff(self, channelList, alpha = 10.0):        
    
        # Getting all the records for the night
        nightQuery = self.session.query(Propagation.similarity)\
                                        .join(TransientEvent, TransientEvent.ID == Propagation.transientEventID)\
                                        .filter(and_(TransientEvent.psgNight == self.night,
                                                     TransientEvent.eventName == self.eventName))         
        
        # Keeping only offset records        
        nightQuery = nightQuery.filter(Propagation.offset != 0.0)    

        for testChannel in channelList :  
            testQuery = nightQuery.filter(Propagation.sinkChannelName==testChannel)    
            for refChannel in channelList:
                if testChannel == refChannel:
                    continue

                similarities = testQuery.filter(Propagation.sourceChannelName==refChannel).all()       

                propRel = self.session.query(PropagationRelationship)\
                                                .filter_by(psgNight          = self.night,
                                                           eventName         = self.eventName,
                                                           sourceChannelName = refChannel,
                                                           sinkChannelName   = testChannel).one()     
                
                if len(similarities):                    
                    propRel.cutoff = percentile(similarities, 100-alpha)
                else:
                    propRel.cutoff = 1.0
          
            
    
    ###########################################################################    
    # FALSE DETECTION AND OUTLIERS REJECTION
    ###########################################################################
    
    def propagationRejection(self, channelList, alpha):        

    
        # Getting all the records for the night
        nightQueryProp   = self.session.query(Propagation)\
                                        .join(TransientEvent, TransientEvent.ID == Propagation.transientEventID)\
                                        .filter(TransientEvent.psgNight == self.night)\
                                        .filter(TransientEvent.eventName == self.eventName)
        
        nightQueryCutoff = self.session.query(PropagationRelationship)\
                                        .filter_by(psgNight  = self.night,
                                                   eventName = self.eventName)
        
        # Keeping only synchronous records        
        nightQueryProp = nightQueryProp.filter(Propagation.offset == 0.0)     

        for testChannel in channelList :  
            testQueryProp   = nightQueryProp.filter(Propagation.sinkChannelName==testChannel)    
            testQueryCutoff = nightQueryCutoff.filter(PropagationRelationship.sinkChannelName==testChannel)    
            
            if self.verbose :
                print "test:", testChannel                 
            
            for refChannel in channelList:
                if testChannel == refChannel:
                    continue
            
                propRel = testQueryCutoff.filter(PropagationRelationship.sourceChannelName==refChannel).all()                   

                # Cutoff should be either equal to a float or to []
                # The "if" clause is necessary to avoid using len() on
                # a fload, raising an error.
                if len(propRel) == 0: 
                    continue 
                elif  len(propRel) == 1:
                    propRel = propRel[0]
                    if not isinstance(propRel, PropagationRelationship): 
                        raise TypeError("propRel = " + str(propRel) + " with type" + str(type(propRel)))                 
                else:
                    raise TypeError("propRel = " + str(propRel) + " with type" + str(type(propRel)))                                       
    
                cutoff = propRel.cutoff
                testRefQueryProp = testQueryProp.filter(Propagation.sourceChannelName==refChannel)   
                validQuery       = testRefQueryProp.filter(Propagation.similarity >= cutoff)                
                            
                # Rejection because of a too low similarity                                                                                   
                N      = testRefQueryProp.count()
                Nvalid = validQuery.count()                 
                if Nvalid == 0: 
                    continue
                
                # Expected false detection rate
                propRel.FDR = (alpha/100.0)/(float(Nvalid)/float(N))            
                
                outlierMin, outlierMax = getOutlierThresholds([prop.delay for prop in validQuery.all()])
                
                nbOut   = 0
                nbValid = 0
                for prop in validQuery.all(): 
                    prop.isFP       = False
                    prop.isOutlier  = prop.delay < outlierMin or  prop.delay > outlierMax
                    
                    if prop.isOutlier:
                        nbOut   += 1
                    else:
                        nbValid += 1                        

                propRel.nbValid = nbValid    
                propRel.nbOut   = nbOut
                self.session.flush()

        self.session.commit()
        

    
            
            
            
    """        
       Identify the spindle propagation fields (SPF). That is, it makes links between
       detected spindles on separated electrodes to postulate, from different observations
       on different channel, a same physiological entity (i.e., the SPF).
    """     
    def identifySPF(self) :
            
            
        
        # Retourne la liste des indices, parmis candidateInds, qui sont des électrodes
        # liées avec l'électrode référence de refRow.
        def getLinkedInd(refRow, candidateInds):
        

            refs  = sources[candidateInds]
            tests = sinks[candidateInds]
            
            
            lastLinked = array([props[refRow].sourceChannelName], dtype=str)
            lastLinked  = unique(concatenate((lastLinked, 
                                               tests[in1d(refs, lastLinked)],
                                               refs[in1d(tests, lastLinked)]  
                                             )))  
                                                  
                                          
            newLinked = unique(concatenate((lastLinked, 
                                               tests[in1d(refs, lastLinked)],
                                               refs[in1d(tests, lastLinked)]  
                                             )))  
                                          
            # Les lignes ci-dessus pourraient être intégrées dans la boucle ci-dessous. Cependant, cela 
            # causerait de évaluation inutile de la condition du while, ce que l'on veut éviter puisque le
            # temps d'exécution de cette fonction est critique.
            while not all(in1d(newLinked, lastLinked)) :
                lastLinked = newLinked
                newLinked = unique(concatenate((newLinked, 
                                               tests[in1d(refs, lastLinked)],
                                               refs[in1d(tests, lastLinked)]  
                                             )))   
                                          
                                         
                                          
            # Non nécessaire de vérifier les références ET les électrodes tests étant donné la définition même
            # de ce que sont des électrodes liées.
            #return(candidateInds[inlist(filteredData$elect, lastLinked) | inlist(filteredData$ref, lastLinked)])
            return array(candidateInds[in1d(tests, lastLinked)])
        
        
        
        
        
        # Vérifie s'il y a une relation bidirectionnelle entre l'électrode ref et source
        # de noRow et si oui retourne l'indice de la ligne bidirectionnelle avec refRow
        def findBidirectionnality(refRow, candidateInds):
        
            # We make sure refRow is not in CandidateInds, else the function would trivially
            # return refRow
            if refRow in candidateInds:
                candidateInds = candidateInds[candidateInds != refRow]
                
            if len(candidateInds) == 0:
                return None
        

            # S'il y a plusieurs candidats, on ne garde que celui dont l'occurence est 
            # la plus proche de celle de noRow. Cela correspond
            # au premier indice de la liste.  
            IND = where(array((sources[refRow] == sinks[candidateInds])*
                            (sinks[refRow] == sources[candidateInds])))[0]
                                            
            if len(IND) == 0:
                return None
            else:
                return candidateInds[IND[0]]
        





        #######################################################################
        # START OF THE FUNCTION BODY

        # Getting all the records for the night
        nightQueryProp_unsorted   = self.session.query(Propagation)\
                                        .join(TransientEvent, TransientEvent.ID == Propagation.transientEventID)\
                                        .filter(TransientEvent.psgNight == self.night)\
                                        .filter(TransientEvent.eventName == self.eventName)\
                                        .filter(Propagation.isFP == False)\
                                        .filter(Propagation.isOutlier == False)     

        nightQueryProp   = nightQueryProp_unsorted.order_by(TransientEvent.startTime)           

        props = nightQueryProp.all()     

        sources  = array([prop.sourceChannelName for prop in props], dtype=str)
        sinks = array([prop.sinkChannelName for prop in props], dtype=str)

        delta                = 0.5
        
        # Compute the widht of the window in which we will be looking for coocuring 
        # sleep spindles. 
        startTimes = array([prop.transientEvent.startTime for prop in props])
        N = nightQueryProp.count()
        if N == 0:
            return
        
        Window = min(10, N)
        while min(startTimes[Window:N] - startTimes[0:(N-Window)]) < delta :
          Window += 1
          if Window >= N:
              break

        for noRow in range(N):
            
            rowBiDir = props[noRow].bidirect        
            if rowBiDir > -1 :
    
                # On prend le temps moyen de propagation calculé de X à Y ce qui correspond à la moyenne
                # du délais de X à Y et de la valeur négative du délais de Y à X.
                assert(props[noRow].delay > 0 and props[rowBiDir].delay > 0)
                props[rowBiDir].delay += props[noRow].delay
                props[rowBiDir].delay /= 2.0
         
            else:
    
                # Une formulation comme 
                # inds = which(DataInd$timeDiff[noRow:length(DataInd$timeDiff)] - DataInd$timeDiff[noRow] <= delta) + noRow-1
                # utilise beaucoup de temps de processeur, d'où l'importance d'utiliser une fenêre. Celle-ci nécessite
                # bien sûr que les valeurs de DataInd$timeDiff soit triées en ordre croissant.
                inds = where(startTimes[noRow:min(noRow+Window, N)] - 
                                                        startTimes[noRow] <= delta )[0] + noRow
                
                # On doit garder seulements les entrées liées à l'électrode Ref de noRow.
                inds = getLinkedInd(noRow, inds)            
    
    
                # Si une ou plusieurs sources sont présumés pour noRow, on conserve celle dont 
                # le temps d'occurence est le plus près de celui de noRow.  Celle-ci sera la dernière 
                # découverte, ce qui ne nous oblige qu'à garder la dernière en mémoire.  
                # Notons qu'une source doit avoir une référence différente de celle de noRow
                for ind in inds:
                    if props[noRow].sourceChannelName != props[ind].sourceChannelName :
                        props[ind].source = props[noRow].no


      
                # Enlever bidirectionnalité :
                # Si une ou plusieurs relations bidirectionnelles ont étées détectées avec
                # noRow, on considère celle dont le temps d'occurence est le plus près de celui
                # de noRow. Celle-ci sera la dernière découverte, ce qui ne nous oblige 
                # qu'à garder la dernière en mémoire. On réécrit donc directement sur toute relations
                # bidirectionnelles ayant été déjà trouvées.
                indBidir = findBidirectionnality(noRow, inds)
                if not indBidir is None :
                    # We mark one (the first one) row as being bidirectionnal
                    # by putting its bidirect attribute to -2. This row will be kept
                    # and its delay will be the average of the two delays of the
                    # bidirectionnal relationship. Labelling it -2 distinguished by
                    # the non bidirectionnal channel pairs which are all labels -1
                    props[noRow].bidirect    = -2   
                    
                    # To mark the other (the second one) ow of the bidirectionnal
                    # relationship as being linked to the first on by indicating
                    # the row of the first row as its bidirect attribute. Once it
                    # has been used to compute the average delay, this second row
                    # will be removed.
                    props[indBidir].bidirect = noRow  
    
    
        self.session.flush()    
       
       
        # SQLAlchemy do not allow to perform Query.update() when order_by() has 
        # been called. Thus, for the following part we ue unsorted data.
        """
        DataSourceM1 = nightQueryProp_unsorted.filter(Propagation.source == -1) 
        dataTmp = set(((prop.transientEvent.startTime, prop.sourceChannelName) for prop in DataSourceM1.all() ))
        SPFs = [SPF(psgNight = self.night)]*len(dataTmp)
        self.session.add_all(SPFs)
        self.session.flush()    
        
        for aSPF, (startTime, ref) in zip(SPFs, dataTmp) :
            props = DataSourceM1\
                        .filter(and_(TransientEvent.startTime == startTime,\
                                     Propagation.sourceChannelName == ref))\
                        .all()
            for prop in props:
                prop.noSPF = aSPF.no


            #DataSourceM1\
            #    .filter(and_(TransientEvent.startTime == startTime,\
            #                 Propagation.sourceChannelName == ref))\
            #    .update({Propagation.noSPF: aSPF.no})

        self.session.flush()    
        print datetime.now()
                
    
        whileQuery = nightQueryProp_unsorted.filter(Propagation.noSPF == -1)
        while whileQuery.count() :
            for prop in whileQuery.all():   
                prop.noSPF = nightQueryProp_unsorted.filter(Propagation.no == prop.source).one().noSPF            
            whileQuery = nightQueryProp_unsorted.filter(Propagation.noSPF == -1)
        """
        self.session.commit()
               
           
    
    
    """
     When a spindle has been detected on two channels and that for both these
     detection, a propagation has been detected between these channels, we are
     left with two records related to a same physiological propagation. We therefore
     remove this duplicate. The duplicate record to remove have been labeled by 
     a value -1 in the field "bidirect" during the execution of identifySPF().
    """
    def removeBidirectionnalDuplicates(self):
        
        self.session.query(Propagation).filter(Propagation.bidirect >= 0).delete()
        self.session.query(Propagation).update({Propagation.bidirect: -(Propagation.bidirect +1)},
                                               synchronize_session=False)
        self.session.commit()        
        

    
    
    
    """
     CORRECTING DELAYS TO HAVE ONLY PROPAGATION WITH POSITIVE TIME DELAYS 
    """
    def negativeDelayCorrection(self): 
        
        
        negQuery = self.session.query(Propagation).filter(Propagation.delay < 0.0 )
        negQuery.update({Propagation.inverted           : True, 
                         Propagation.delay              : -Propagation.delay,
                         Propagation.sourceChannelName  : Propagation.sinkChannelName,
                         Propagation.sinkChannelName    : Propagation.sourceChannelName},
                         synchronize_session=False)
        self.session.commit()        





    """
      Compute robust central tendencies of observationVariables, aggregated at 
      aggeragationlevels.
    """         
    """
    def computeAveragePropagation(self, aggeragationlevels, observationVariables, minNbValid=40):
    
        def computeCellAverages(data, observationVariables, 
                                aggeragationlevels, levelInstanciation):
    
            dataOut = {}
            # For the N variables ....
            NbTotal = len(data)
            for observationVariable in observationVariables:
                X          = array(data[observationVariable])
                notOutInd  = getInliers(X, coef=1.5)
                nbOut      = NbTotal-len(notOutInd)
                
                validX     = X[notOutInd]    
                validX     = validX[logical_not(isnan(validX))]                    
                
                #valids = !(X %in% outTau$out) & !is.na(X)
                if len(validX) >= minNbValid :   # criterion c4
                    try:
                        # MinCovDet().fit() crashes if we don't make this reshape...
                        validX = reshape(validX, (len(validX), 1))
                        covX    = MinCovDet().fit(validX)
                        meanX   = covX.location_[0]
                        sdX     = sqrt(covX.covariance_[0, 0])
                    except ValueError:
                        # When the spindle detection resulted in discretized values
                        # of spindle characteristics, the distribution of the variable
                        # might not be acceptable for MCD algorithm which craches
                        # with a "ValueError: expected square matrix" error. In these
                        # time we fall back to simple statistics. 
                        meanX   = median(validX)
                        sdX     = std(validX)                        
                else:
                    meanX   = nan
                    sdX     = nan     
                    
                dataOut[observationVariable + "_mean"]    = meanX
                dataOut[observationVariable + "_sd"]      = sdX
                dataOut[observationVariable + "_nbOut"]   = nbOut
                dataOut[observationVariable + "_nbValid"] = len(validX)
            
            for key in levelInstanciation:
                dataOut[key] = levelInstanciation[key]
            
            for level in aggeragationlevels:
                dataOut[level] = str(data[level].iloc[0])
                
            return DataFrame(dataOut, [1])   
            
            
                   
        ##########################################################################
        # Code execution   

        levelInstanciation = {"night":self.night}


        query = self.session.query(Propagation)\
                        .join(TransientEvent, TransientEvent.ID == Propagation.transientEventID)\
                        .filter(and_(TransientEvent.psgNight  == self.night,
                                     TransientEvent.eventName == self.eventName,  
                                     Propagation.isFP == False,
                                     Propagation.isOutlier == False))        

        df  = rows2df(query.all())
        dat = df.groupby(aggeragationlevels).apply(lambda x: computeCellAverages(x, observationVariables, aggeragationlevels, levelInstanciation))
        
        for norow in range(dat.shape[0]):
            row = dict(dat.iloc[norow])       
            
            propRel = self.session.query(PropagationRelationship)\
                            .filter_by(sinkChannelName   = row["sinkChannelName"],
                                       sourceChannelName = row["sourceChannelName"],  
                                       psgNight          = row["night"],
                                       eventName         = self.eventName).one()             

            del row["sinkChannelName"]
            del row["sourceChannelName"]
            del row["night"]

            propRel.update(self.session, row, behavior="addMissingFields")
    """  




    """
      Compute robust central tendencies of observationVariables, aggregated at 
      aggeragationlevels.
    """            
    def computeAveragePropagation(self, deltaWindow = 0.5, alphaSD = 0.2, minNbValid=40):
    
        propRels = self.session.query(PropagationRelationship)\
                                        .filter_by(psgNight  = self.night,
                                                   eventName = self.eventName).all() 

        for propRel in propRels:
            
            if propRel.nbValid >= 2 :
                try:
                    props   = propRel.getValidPropagations(self.session)
                    delays  = [prop.delay for prop in props]
                    covX    = MinCovDet().fit(array(delays))
                    meanX   = covX.location_[0]
                    sdX     = sqrt(covX.covariance_[0, 0])
                except ValueError:
                    # When the spindle detection resulted in discretized values
                    # of spindle characteristics, the distribution of the variable
                    # might not be acceptable for MCD algorithm which craches
                    # with a "ValueError: expected square matrix" error. In these
                    # time we fall back to simple statistics. 
                    meanX   = median(delays)
                    sdX     = std(delays)                        
            else:
                meanX   = nan
                sdX     = nan     
            
            propRel.delay_mean = meanX
            propRel.delay_sd   = sdX
            
            propRel.testRejectionC3(deltaWindow, alphaSD)
            propRel.testRejectionC4(minNbValid)
                  
        self.session.commit()




    def getAveragePropagation(self, applyC3=True, applyC4=True):
        query = self.session.query(PropagationRelationship)\
                        .filter_by(psgNight  = self.night,
                                   eventName = self.eventName)                  
        
        if applyC3:
            query = query.filter_by(isValidC3 = True)
        if applyC4:
            query = query.filter_by(isValidC4 = True)
        
        return rows2df(query.all())
        
        """
        import collections
        result = collections.defaultdict(list)
        for d in [propRel.getDict(self.session) for propRel in propRels]:
            for k, v in d.items():
                result[k].append(v)
                
        return pd.DataFrame(result)
        """
    
    
    
    
    
    
    
    
    
    
    
    
    """
     Helper function used to provide an easy interface for computing SPF of a 
     set of files.
    """
    """
    def computeSPF_allFiles(path, offsetFilePattern="spindleDelays_offset_*.bdf_*.txt", 
                   includeElectrodePatterns=['Fp1', 'Fp2', 'F7', 'F8', 'F3', 'F4', 
                                             'T3', 'T4', 'C3', 'C4', 'T5', 'T6', 'P3', 
                                             'P4', 'O1', 'O2', 'Fz', 'Cz', 'Pz', 'Oz'],
                   excludeElectrodePatterns=[], alpha = 10.0, verbose = True):
    
        
        files = glob(path + offsetFilePattern)
        rePattern =  offsetFilePattern.replace("*", "|").replace(".", '\.')
            
        tmp = map(lambda f: re.split(rePattern, f)[1:3], files)
        nights = [record[0] for record in tmp]
        channels = [record[1] for record in tmp]
        
        # Keeping only the part of the channel name identifying the active electrode.
        electrodes = map(lambda channel: getActiveElectrode(channel, includeElectrodePatterns, excludeElectrodePatterns), channels)
    
        tokens =  offsetFilePattern.split("*")
        
        # for each subject in the file list
        for night in unique(nights):
            
            channelList     = [channel for strNight, channel in zip(nights, channels) if strNight == night ]
            electrodeList   = [electrode for strNight, electrode in zip(nights, electrodes) if strNight == night ]
            
            fileNameSyncList  = []
            fileNameAsyncList = [] 
            # For eletrodes of this recording night (these are the reference electrodes)
            for refElectrode, refChannel in zip(electrodeList, channelList):
                
                if verbose :
                    print "ref:", refElectrode        
                
                fileNameSyncList.append(tokens[0][:-7] + night + tokens[1] + refChannel + tokens[2])
                fileNameAsyncList.append(tokens[0] + night + tokens[1] + refChannel + tokens[2])
        
            computeSPF(path, night, fileNameSyncList, fileNameAsyncList, channelList, electrodeList, alpha, verbose)        
            
           
    """



            

    