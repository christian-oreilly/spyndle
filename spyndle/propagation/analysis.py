# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 14:06:58 2013

@author: REVESTECH
"""


import re

from pandas import DataFrame, read_csv
from glob import glob
from scipy import unique, where, percentile, concatenate, array, mod, arange
from scipy import isnan, sqrt, nan, logical_not, reshape
from sklearn.covariance import MinCovDet


###############################################################################
# Internal functions
###############################################################################


# Retourne la liste des indices, parmis candidateInds, qui sont des électrodes
# liées avec l'électrode référence de refRow.
def getLinkedInd(Data, refRow, candidateInds):
    filteredData =  Data.iloc[candidateInds, ]
    lastLinked = array([Data["ref"].iloc[refRow]], dtype=str)
    lastLinked  = unique(concatenate((lastLinked, 
                               array(filteredData["test"][[ref  in lastLinked for ref  in filteredData["ref" ]]], dtype=str),
                               array(filteredData["ref" ][[test in lastLinked for test in filteredData["test"]]], dtype=str)
                              )))  
                                  
    newLinked = unique(concatenate((lastLinked, 
                           array(filteredData["test"].iloc[[ref  in lastLinked for ref  in filteredData["ref" ]]], dtype=str),
                           array(filteredData["ref" ].iloc[[test in lastLinked for test in filteredData["test"]]], dtype=str)
                          )))  
                                  
    # Les lignes ci-dessus pourraient être intégrées dans la boucle ci-dessous. Cependant, cela 
    # causerait de évaluation inutile de la condition du while, ce que l'on veut éviter puisque le
    # temps d'exécution de cette fonction est critique.
    while not all([newL in lastLinked for newL in newLinked]) :
        lastLinked = newLinked
        newLinked = unique(concatenate((newLinked, 
                           array(filteredData["test"].iloc[[ref  in lastLinked for ref  in filteredData["ref" ]]], dtype=str),
                           array(filteredData["ref" ].iloc[[test in lastLinked for test in filteredData["test"]]], dtype=str)
                          )))  
                                  
    # Non nécessaire de vérifier les références ET les électrodes tests étant donné la définition même
    # de ce que sont des électrodes liées.
    #return(candidateInds[inlist(filteredData$elect, lastLinked) | inlist(filteredData$ref, lastLinked)])
    return array([candidateInd for candidateInd, test in zip(candidateInds, filteredData["test"]) if test in lastLinked])



# Vérifie s'il y a une relation bidirectionnelle entre l'électrode ref et source
# de noRow et si oui retourne l'indice de la ligne bidirectionnelle avec refRow
def findBidirectionnality(Data, refRow, candidateInds):
  # S'il y a plusieurs candidats, on ne garde que celui dont l'occurence est 
  # la plus proche de celle de noRow. Cela correspond
  # au premier indice de la liste.  
  IND = where(array((Data["ref" ].iloc[refRow] == Data["test"].iloc[candidateInds])*(Data["test"].iloc[refRow] == Data["ref" ].iloc[candidateInds])))[0]
  if len(IND) == 0:
      return []
  else:
      return candidateInds[IND[0]]





def getActiveElectrode(electrodeNames, excludePatterns, channelLabel):
    
    retName = ""
    for electrodeName in electrodeNames:
        if electrodeName in channelLabel :
            if retName == "":
                retName = electrodeName
            else:
                raise "Conficting electrode names."

    for pattern in excludePatterns:
        if pattern in channelLabel :
            return ""
            
    return retName



"""
  Compute outliers according to the nonparametric test proposed by Tukey (1977).
  
  Return the indexes of the outliers. 
  
  Tukey JW. Exploratory data analysis. Addison-Wesley Pub. Co.: Reading, Mass. ; 
  Don Mills, Ont., 1977.
"""
def getOutliers(data, coef=1.5):
    Q1   = percentile(data, 25.0)
    Q3   = percentile(data, 75.0)
    return where((data > Q3 + coef*(Q3-Q1))*(data < Q1 - coef*(Q3-Q1)))[0]

def getNotOutliers(data, coef=1.5):
    Q1   = percentile(data, 25.0)
    Q3   = percentile(data, 75.0)
    return where((data <= Q3 + coef*(Q3-Q1))*(data >= Q1 - coef*(Q3-Q1)))[0]


###########################################################################    
# FALSE DETECTION AND OUTLIERS REJECTION
###########################################################################
def propagationRejection(path, night, nights, channels, electrodes, alpha, verbose = True):        

    if verbose :
        print "night:", night
    
    DataProp        = DataFrame()
    channelList     = [channel for strNight, channel in zip(nights, channels) if strNight == night ]
    electrodeList   = [electrode for strNight, electrode in zip(nights, electrodes) if strNight == night ]
    # For eletrodes of this recording night (these are the reference electrodes)
    for refElectrode, refChannel in zip(electrodeList, channelList):
        
        if verbose :
            print "ref:", refElectrode        
        
        sync  = read_csv(path + "spindleDelays_" + night + ".bdf_" + refChannel + ".txt", sep=";")
        async = read_csv(path + "spindleDelays_offset_" + night + ".bdf_" + refChannel + ".txt", sep=";") 


        for testElectrode, testChannel in zip(unique(electrodeList), unique(channelList)) :
                     
            if testElectrode != refElectrode: # and not any(is.na(async[nameSim])) :
                
                nameSim   = "similarity_" + testChannel   
                nameDelay = "delay_" + testChannel   
                   
                # Rejection because of a too low similarity
                cutoff    = percentile(async[nameSim], 100-alpha)
                IND1      = where(sync[nameSim] >= cutoff)[0] 
                sim       = sync[nameSim][IND1]
                delays    = sync[nameDelay][IND1]
                
                # Expected false detection rate
                FDR       = alpha/(float(len(sim))/len(sync[nameSim]))
                
                # Rejection of delays outliers
                IND2      = getNotOutliers(delays)
                sim       = sim.iloc[IND2]
                delays    = delays.iloc[IND2]                          
                          
                if len(sim):  
                    # Keep all the fields exctep the delay and similarity related ones
                    colNames =  [colName for colName in sync.columns.tolist() if not ("delay_" in colName or "similarity_" in colName)] 
                    newData = sync[colNames].iloc[IND1[IND2]]
                    
                    # Add the other fields
                    newData["night"]        = night
                    newData["ref"]          = refElectrode
                    newData["test"]         = testElectrode
                    newData["delay"]        = delays
                    newData["similarity"]   = sim
                    newData["FDR"]          = FDR

  
                    DataProp = DataProp.append(newData, ignore_index=True)
        
    return DataProp

        
        
        
"""        
   Identify the spindle propagation fields (SPF). That is, it makes links between
   detected spindles on separated electrodes to postulate, from different observations
   on different channel, a same physiological entity (i.e., the SPF).
"""     
def identifySPF(DataProp) :
        
    DataProp             = DataProp.sort("startTime")
    DataProp             = DataProp.rename(dict(zip(DataProp.index.tolist(), range(len(DataProp)))))
        
    delta                = 0.5
    DataProp["source"]   = -1
    DataProp["bidirect"] = -1
    
    # Compute the widht of the window in which we will be looking for coocuring 
    # sleep spindles. 
    N = len(DataProp)
    Window = min(10, N)
    while min(DataProp["startTime"].iloc[Window:N] - DataProp["startTime"].iloc[0:(N-Window)]) < delta :
      Window += 1
      if Window >= N:
          break
    
    
    for noRow in range(N):
        if mod(noRow, 1000) == 0:
            print noRow, N
        
        rowBiDir = DataProp["bidirect"].iloc[noRow]        
        if rowBiDir > -1 :

            # On prend le temps moyen de propagation calculé de X à Y ce qui correspond à la moyenne
            # du délais de X à Y et de la valeur négative du délais de Y à X.
            DataProp["delay"].iloc[rowBiDir] = (DataProp["delay"].iloc[rowBiDir] - DataProp["delay"].iloc[noRow])/2.0
     
        else:

            # Une formulation comme 
            #inds = which(DataInd$timeDiff[noRow:length(DataInd$timeDiff)] - DataInd$timeDiff[noRow] <= delta) + noRow-1
            # utilise beaucoup de temps de processeur, d'où l'importance d'utiliser une fenêre. Celle-ci nécessite
            # bien sûr que les valeurs de DataInd$timeDiff soit triées en ordre croissant.
            inds = where(DataProp["startTime"].iloc[noRow:min(noRow+Window, N)] - 
                         DataProp["startTime"].iloc[noRow] <= delta)[0] + noRow
            
            # On doit garder seulements les entrées liées à l'électrode Ref de noRow.
            inds = getLinkedInd(DataProp, noRow, inds)            



            # Si une ou plusieurs sources sont présumés pour noRow, on conserve celle dont 
            # le temps d'occurence est le plus près de celui de noRow.  Celle-ci sera la dernière 
            # découverte, ce qui ne nous oblige qu'à garder la dernière en mémoire.  
            # Notons qu'une source doit avoir une référence différente de celle de noRow
            condition = array(DataProp["ref"].iloc[noRow] != DataProp["ref"].iloc[inds])
            if isinstance(condition, bool) or isinstance(condition, list) :
                DataProp["source"].iloc[inds[condition]] = noRow
  
  
            # Enlever bidirectionnalité pour éviter que les délais se cancellent entre eux :
            # Si une ou plusieurs relations bidirectionnelles ont étées détectées avec
            # noRow, on considère celle dont le temps d'occurence est le plus près de celui
            # de noRow. Celle-ci sera la dernière découverte, ce qui ne nous oblige 
            # qu'à garder la dernière en mémoire. On réécrit donc directement sur toute relations
            # bidirectionnelles ayant été déjà trouvées.
            DataProp["bidirect"].iloc[findBidirectionnality(DataProp, noRow, inds)] = noRow

    
    DataProp["noSpindle"] = -1
    DataSourceM1 = DataProp[DataProp["source"] == -1]
    dataTmp = DataSourceM1[["startTime", "ref"]].drop_duplicates()
    for iRowTmp in range(len(dataTmp)) :
        DataSourceM1["noSpindle"][ (DataSourceM1["startTime"] == dataTmp["startTime"].iloc[iRowTmp])* 
                      (DataSourceM1["ref"]       == dataTmp["ref"].iloc[iRowTmp])] = iRowTmp   
    
    DataProp[DataProp["source"] == -1] = DataSourceM1

    while any(DataProp["noSpindle"] == -1) :
        iRowTmp = arange(len(DataProp))[DataProp["noSpindle"] == -1]
        DataProp["noSpindle"].iloc[iRowTmp] = array(DataProp["noSpindle"].iloc[array(DataProp["source"].iloc[iRowTmp], dtype=int)])

    return DataProp


"""
 When a spindle has been detected on two channels and that for both these
 detection, a propagation has been detected between these channels, we are
 left with two records related to a same physiological propagation. We therefore
 remove this duplicate. The duplicate record to remove have been labeled by 
 a value -1 in the field "bidirect" during the execution of identifySPF().
"""
def removeBidirectionnalDuplicates(dataProp):
    
    return dataProp.drop(dataProp["bidirect"] == -1)




"""
 CORRECTING DELAYS TO HAVE ONLY PROPAGATION WITH POSITIVE TIME DELAYS 
"""
def negativeDelayCorrection(dataProp): 
            
    indNeg                     = dataProp["delay"] < 0.0   
    test                       = dataProp["test"][indNeg]
    ref                        = dataProp["ref"][indNeg]
    dataProp["delay"][indNeg]  = -dataProp["delay"][indNeg]
    dataProp["ref"][indNeg]    = test
    dataProp["test"][indNeg]   = ref

    return dataProp










###############################################################################
# External functions.
###############################################################################

"""
 Compute spindle propagation fields.
 TODO:  Maybe object which contain the rules for obtaining the files and 
        its information (suject, nights, etc.) should be passed instead.
"""
def computeSPF(path, offsetFilePattern="spindleDelays_offset_*.bdf_*.txt", 
               includeElectrodePatterns=['Fp1', 'Fp2', 'F7', 'F8', 'F3', 'F4', 
                                         'T3', 'T4', 'C3', 'C4', 'T5', 'T6', 'P3', 
                                         'P4', 'O1', 'O2', 'Fz', 'Cz', 'Pz', 'Oz'],
               excludeElectrodePatterns=[], alpha = 10.0, verbose = True):

    
    files = glob(path + offsetFilePattern)
    rePattern =  offsetFilePattern.replace("*", "|").replace(".", '\.')
        
    tmp = map(lambda file: re.split(rePattern, file)[1:3], files)
    nights = [record[0] for record in tmp]
    channels = [record[1] for record in tmp]
    
    # Keeping only the part of the channel name identifying the active electrode.
    electrodes = map(lambda electrode: getActiveElectrode(includeElectrodePatterns, 
                                                 excludeElectrodePatterns, electrode), channels)

    # for each subject in the file list
    for night in unique(nights):
        
        # FALSE DETECTION AND OUTLIERS REJECTION
        dataProp = propagationRejection(path, night, nights, channels, electrodes, alpha, verbose)        
    
        # COMPUTING SPINDLE PROPAGATION FIELD      
        dataProp = identifySPF(dataProp)            
    
    
        # REMOVING DUPLICATE PROPAGATION DUE TO BIDIRECTIONNALITY
        if verbose: 
            print "removing duplicates"        
        dataProp = removeBidirectionnalDuplicates(dataProp)  
    
                   
        # CORRECTING DELAYS TO HAVE ONLY PROPAGATION WITH POSITIVE TIME DELAYS 
        if verbose: 
            print "correctin delays"   
        dataProp = negativeDelayCorrection(dataProp)

        # Saving the result.
        dataProp.to_csv(path + "correctedData_" + night + ".csv", sep=";")
        
               
               
               
               
               
               
               
"""
  Compute robust central tendencies of observationVariables, aggregated at 
  aggeragationlevels.
"""            
def computeAveragePropagation(path, aggeragationlevels,
                              observationVariables, pattern="correctedData_*.csv", verbose=True):

    def computeCellAverages(data, observationVariables, aggeragationlevels, levelInstanciation):

        dataOut = {}
        # For the N variables ....
        NbTotal = len(data)
        for observationVariable in observationVariables:
            X          = array(data[observationVariable])
            notOutInd  = getNotOutliers(X, coef=1.5)
            nbOut      = NbTotal-len(notOutInd)
            
            validX     = X[notOutInd]    
            validX     = validX[logical_not(isnan(validX))]                    
            
            #valids = !(X %in% outTau$out) & !is.na(X)
            if len(validX) >= 10 :
                # MinCovDet().fit() crashes if we don't make this reshape...
                validX = reshape(validX, (len(validX), 1))
                
                covX    = MinCovDet().fit(validX)
                meanX   = covX.location_[0]
                sdX     = sqrt(covX.covariance_[0, 0])
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
        
        
    """                  
    def computeCellAverages(data, observationVariables, levelInstanciation):

        dataOut = {}
        # For the N variables ....
        NbTotal = len(data)
        for observationVariable in observationVariables:
            X          = array(data[observationVariable])
            notOutInd  = getNotOutliers(X, coef=1.5)
            nbOut      = NbTotal-len(notOutInd)
            
            validX     = X[notOutInd]    
            validX     = validX[logical_not(isnan(validX))]                    
            
            #valids = !(X %in% outTau$out) & !is.na(X)
            if len(validX) >= 10 :
                # MinCovDet().fit() crashes if we don't make this reshape...
                validX = reshape(validX, (len(validX), 1))
                
                covX    = MinCovDet().fit(validX)
                meanX   = covX.location_[0]
                sdX     = sqrt(covX.covariance_[0, 0])
            else:
                meanX   = nan
                sdX     = nan     
                
            dataOut[observationVariable + "_mean"]    = meanX
            dataOut[observationVariable + "_sd"]      = sdX
            dataOut[observationVariable + "_nbOut"]   = nbOut
            dataOut[observationVariable + "_nbValid"] = len(validX)
        
        for key in levelInstanciation:
            dataOut[key] = levelInstanciation[key]
        
        return dataOut   
        


        

        
        
    def recursiveExecution(data, aggeragationlevels, observationVariables, levelInstanciation={}):                              
    
        if len(aggeragationlevels) :
            meanData = []
            for item in unique(data[aggeragationlevels[0]]) :    
                levelInstanciation[aggeragationlevels[0]] = item
                meanData.extend(recursiveExecution(data[data[aggeragationlevels[0]] == item], 
                                               aggeragationlevels[1:], observationVariables, 
                                                levelInstanciation))
            return meanData
        else:
            return [computeCellAverages(data, observationVariables, levelInstanciation)]
    """
                 
                        
    ##########################################################################
    # Code execution   
    files = glob(path+pattern)
    
    meanData = DataFrame() 
    for f in files:            
        rePattern =  pattern.replace("*", "|").replace(".", '\.')
        night = re.split(rePattern, f)[1]         
        
        if verbose:
            print night        
        levelInstanciation = {"night":night}
        #from datetime import datetime        
        #a= datetime.now()
        #dataTmp = recursiveExecution(read_csv(f, sep=";"), aggeragationlevels, 
        #                             observationVariables, levelInstanciation) 
        #meanData =  meanData.append(DataFrame(dataTmp), ignore_index=True)               
        #b= datetime.now()
        dat          = read_csv(f, sep=";").groupby(aggeragationlevels).apply(lambda x: computeCellAverages(x, observationVariables, aggeragationlevels, levelInstanciation))
        meanData     = meanData.append(dat, ignore_index=True)           
        #print len(dat), dat
        
        #c = datetime.now()
        #print b-a, c-b
         
    meanData.to_csv(path + "meanData.csv", sep=";")
    
path= 'C:\\DATA\\Julie\\Results\\'
computeAveragePropagation(path, aggeragationlevels = ["ref", "test"],
                              observationVariables = ["RMSamp", "delay", "duration", "meanFreq", 
                                                      "similarity", "slope", "slopeOrigin"], 
                                                      pattern="correctedData_*.csv",  verbose=True)    