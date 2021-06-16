'''
Created on Jun 7, 2021

@author: andrewdavidson
'''
import logging
from matplotlib import pyplot as plt

import pandas as pd
import numpy as np
from   setupLogging import setupLogging
import upsetplot as up
from utils.upsetPlotData import UpSetPlotData
import unittest

################################################################################
class TestUpsetPlotData(unittest.TestCase):
    configFilePath = setupLogging( default_path='logging.test.ini.json')
    logger = logging.getLogger(__name__)
    logger.info("using logging configuration file:{}".format(configFilePath))   

    ################################################################################
    def createTestData(self):
        self.logger.info("BEGIN")
        
        s1 = set( ['a', 'b', 'c', 'd'] )
        s2 = set( [     'b', 'c',      'e', 'f', 'g', 'h'] )
        s3 = set( [                         'f', 'g',           'j', 'k'] )
        s4 = set( [          'c', 'd',           'g', 'h', 'i', 'j'] )
        
        # make sure we do not have order depend bugs
        setList = [s1,     s4,   s2,   s3]
        setNames = ['s1', 's4', 's2', 's3']
        dataDict = { setNames[i] : setList[i] for i in range(len(setList)) }
        self.logger.info("END\n")  
        
        return( dataDict )
             
    ################################################################################
    def testTemplate(self):
        self.logger.info("BEGIN")
            
        self.logger.info("END\n")  

    ################################################################################
    def testFindSubSets(self):
        self.logger.info("BEGIN")
        
        dataDict = self.createTestData()
        upData = UpSetPlotData( dataDict )
        print(upData.intersectionDict)
        
        expectedSets = {'s1': {'a'}, 
             's4': {'i'}, 
             's2': {'e'}, 
             's3': {'k'}, 
             's1,s2': {'c', 'b'}, 
             's1,s4': {'c', 'd'}, 
             's2,s3': {'g', 'f'}, 
             's2,s4': {'c', 'h', 'g'}, 
             's3,s4': {'j', 'g'}, 
             's1,s2,s4': {'c'}, 
             's2,s3,s4': {'g'}
         }

        for key,intSet in upData.intersectionDict.items():
            self.logger.info( "key:{} i:{}".format(key, sorted(intSet) ) )

        # check interesections are correct
        self.assertDictEqual(expectedSets, upData.intersectionDict, "ASSERT FAILED")
        
        # check the index of the pandas series is correct
        self.logger.info("idx:\n{}".format(upData.plotData.index))
        expectedIndex = pd.MultiIndex.from_tuples([( True, False, False, False),
                                        ( True,  True, False, False),
                                        ( True,  True, False,  True),
                                        ( True, False, False,  True),
                                        (False,  True, False, False),
                                        (False,  True,  True, False),
                                        (False,  True,  True,  True),
                                        (False,  True, False,  True),
                                        (False, False,  True, False),
                                        (False, False,  True,  True),
                                        (False, False, False,  True)],
                                       names=['s1', 's2', 's3', 's4'])
        self.assertTrue( (expectedIndex == upData.plotData.index).all(), "index miss match" )
        
        # check the cardinality is correct. ie. the size of the interesections
        self.logger.info("plotData:\n{}".format(upData.plotData))
        retNumpy = upData.plotData.to_numpy()
        expectedCardinality = np.array([1, 2, 1, 
                                        2, 1, 2, 
                                        1, 3, 1, 
                                        2, 1])
        
        self.assertTrue((expectedCardinality == retNumpy).all(), 'cardinality is wrong')
        
        self.logger.info("END\n")  
        
    ################################################################################
    def testPlot(self):
        self.logger.info("BEGIN")  
        
        dataDict = self.createTestData()
        upData = UpSetPlotData( dataDict )
        
        figureWidthInInches = 8
        figureHeightInInches = 3
        fig = plt.figure(figsize=(figureWidthInInches,figureHeightInInches))
        subPlotDict = up.plot(upData.plotData, fig)
        self.logger.info("subPlotDict keys:{}".format(subPlotDict.keys()))
        # [INFO] subPlotDict keys:dict_keys(['matrix', 'shading', 'totals', 'intersections'])

        
        
        outfile = "testPlot.png"
        fig.savefig(outfile, dpi=300, bbox_inches='tight')

        self.logger.info("END\n")  


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()