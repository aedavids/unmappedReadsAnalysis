
# set -x
# #
# # make script work for both anaconda and minconda
# #
# # $ conda info | grep -i 'base environment'
# # base environment : /Users/andrewdavidson/anaconda3  (writable)
# condaBase=`conda info | grep -i 'base environment' | cut -d : -f 2 | cut '-d ' -f 2`
# # source ~/anaconda3/etc/profile.d/conda.sh
# source ${condaBase}/etc/profile.d/conda.sh
# conda activate extraCellularRNA


#!/usr/bin/env conda run -n extraCellularRNA
# aedavids@ucsc.edu
# 5/14/21

import argparse
import pandas as pd
import sys

########################################################################
class CommandLine() :
    '''
    Handle the command line, usage and help requests.
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        
        arguments:
            inOpst: a list of cli arguments. pass None if you want to use the the
                    true CLI arguments. pass a list if you want to use from a juypter notebook
        '''
        
        proLog = 'BME 237 Applied RNA Bioinformtics'
        epilog = 'sorts a tsv file created by mineSalmonLogs.sh'        
        self.parser = argparse.ArgumentParser(description = proLog,
                                             epilog = epilog, 
                                             add_help = True, #default is True 
                                             prefix_chars = '-' 
                                             # autogenerate usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('inFile', action = 'store', help='input file name')
        #self.parser.add_argument('-o', '--outputFile', action = 'store', help='output file name') 
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)
            
########################################################################
def sortSalmonTSV(inFile):
    df = pd.read_csv(inFile, sep="\t")
    sortCols=['sampleName', 'mappingRate']
    sortedDF = df.sort_values(by=sortCols, ascending=[True, False])
    #print( sortedDF.loc[:,sortCols].head() )
    #sortedDF.loc[:,sortCols].tail()
    stdout = sys.stdout
    sortedDF.to_csv(stdout, sep="\t", index=False)
    
########################################################################
def main(inComandLineArgsList=None):
    '''
    process command line arguments can call createPlot()  
    '''
    if inComandLineArgsList is None:
        cli = CommandLine()
    else :
        cli = CommandLine(inComandLineArgsList)
        
    inFile = cli.args.inFile
    print("inFile:{}".format(inFile))
    
    if inFile == None:
        cli.parser.print_help(sys.stderr)
        sys.exit(-1)
        
    sortSalmonTSV(inFile)
   
########################################################################     
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    main()
