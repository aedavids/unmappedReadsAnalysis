'''
Created on Jun 7, 2021

@author: andrewdavidson
'''
import copy
from itertools import combinations
import logging
import numpy as np
import pandas as pd

################################################################################
class UpSetPlotData(object):
    '''
    classdocs
    
    ref: 
        Notebook pape 219 " monday 6/7 upset plot data"
        https://pandas.pydata.org/pandas-docs/stable/user_guide/advanced.html
        https://upsetplot.readthedocs.io/en/stable/index.html
        https://github.com/jnothman/UpSetPlot
    '''
    logger = logging.getLogger(__name__)


    ################################################################################
    def __init__(self, setDict):
        '''
        Constructor
        '''
        # use deep copy to prevent side effects
        # we change the initial set so that then only contain
        # elements that do not interesect with other sets
        self.intersectionDict = copy.deepcopy(setDict) 
        self.setNames = sorted( list( self.intersectionDict.keys() ) )
        self.plotData = None
        
        self._findPossibleSubSets()
        self._findElementsThatOnlyOccurInASingleSet()
        self._createPlotData()
        
    ################################################################################
    def _calculateIntersection(self, s1Key, s2Key):
        self.logger.debug("BEGIN")
        self.logger.debug("s1Key:{} s2Key:{}".format(s1Key, s2Key))
        s1 = self.intersectionDict[s1Key]
        s2 = self.intersectionDict[s2Key]
        s1IntersectS2 = s1.intersection( s2 )
        self.logger.debug("s1:{}".format(s1))
        self.logger.debug("s2:{}".format(s2))
        self.logger.debug("s1IntersectS2:{}".format(s1IntersectS2))
        
        if len(s1IntersectS2) > 0:
            interesectionKey = s1Key + "," + s2Key
            self.logger.debug("adding interesectionKey:{}".format(interesectionKey))
            self.intersectionDict[ interesectionKey ] = s1IntersectS2
        
        self.logger.debug("END\n")

    ################################################################################
    def _createPlotData(self):
        intersectionSetKeys = sorted(list( self.intersectionDict.keys()))
        numInterections = len(intersectionSetKeys )
        numSets = len(self.setNames)
        indexNP = np.zeros((numInterections, numSets),  dtype=bool)
        
        cardinalityNP = np.zeros(numInterections)
        for i in range(numInterections):
            key = intersectionSetKeys[i]
            tokens = key.split(",")
            for token in tokens:
                j = self.setNames.index(token)
                indexNP[i,j] = True
                
            s = self.intersectionDict[key]
            cardinalityNP[i] = len(s)
            
        df = pd.DataFrame(indexNP, columns=self.setNames)
        idx = pd.MultiIndex.from_frame(df)
        self.plotData = pd.Series(cardinalityNP, index=idx)
        
    ################################################################################
    def _findElementsThatOnlyOccurInASingleSet(self):
        self.logger.debug("BEGIN")
         
        # check elements in original sets against all the pairwise intersections
        choose = 2
        # you can not iterate over generator more than once
        combGenerator = combinations(self.setNames, choose)
        comb = list(combGenerator)
        for key in self.setNames:
            self.logger.debug("\n!!!!!!!!!!!!!!!")            
            self.logger.debug("key:{}".format(key))
            elementsNotInInterSection = self.intersectionDict[key]
            for pairedSet in comb:
                pairedKey = ",".join(pairedSet)
                if pairedKey in self.intersectionDict:
                    pairedIntersection = self.intersectionDict[pairedKey]
                    self.logger.debug("\n****")
                    self.logger.debug("pairedKey:{} {}".format(pairedKey, pairedIntersection))
                    self.logger.debug("       e:{}".format(elementsNotInInterSection))
                    elementsNotInInterSection = elementsNotInInterSection - pairedIntersection
                    self.logger.debug("       e:{}".format(elementsNotInInterSection))                    
                if len(elementsNotInInterSection) == 0:
                    # no need to check remaining paired interesection
                    self.logger.debug("{} no single elements".format(key))
                    break
                
            if len(elementsNotInInterSection) > 0:
                self.intersectionDict[key] = elementsNotInInterSection
             
        self.logger.debug("END\n")
        
    ################################################################################
    def _findPossibleSubSets(self):
        '''
        AEDWIP
        '''
        self.logger.debug("BEGIN")
        
        # avoid side effects while debugging
        for stepi in range(1, len(self.setNames)):
            self.logger.debug("stepi:{}".format(stepi))
            choose = stepi + 1
            combGenerator = combinations(self.setNames, choose)
            # you can not iterate over generator more than once
            # self.logger.debug(list(comb))
            
            for c in combGenerator:
                self.logger.debug( "c:{}".format( c ) )
                # s1Key is set from previous set
                s1Key = ",".join( c[:-1])
                s2Key = c[-1]
                
                self.logger.debug("s1Key:{} s2Key:{}".format(s1Key, s2Key))

                # make sure sets exist
                if s1Key in self.intersectionDict :
                    assert s2Key in self.intersectionDict, "ERROR s2Key:{} not found".format(s2Key)
                    self._calculateIntersection(s1Key, s2Key)
                else:
                    self.logger.debug("skipping intersection key:{} set was null ".format(s1Key))

        self.logger.debug("END\n")

