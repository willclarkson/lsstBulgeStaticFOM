#
# calcFOM.py
#

# WIC started 2020-02-21 - methods to compute the bulge static figure
# of merit by comparing two evaluations stored on disk. (2020-02-21
# Might need to update the namespace for the various modules... should
# this module import the other two, for example?)

import os
import numpy as np
from astropy.table import Table

class fomCalc(object):

    def __init__(self, pathJoined='DUMMY', \
                 lFilts=['g','r','i','z','y'], \
                 magSurplus = -5., \
                 pmMax = 0.8, \
                 spatial_bMax = -1.5, \
                 spatial_bMin = -10., \
                 Verbose=True):

        """Performs the comparison for the 'MAF' evaluations. Expects a map
with MAF evaluations interpolated to the positions of the MSTO
evaluations, and all on the same system (i.e. each row has both parts
of the comparison)."""

        # Note that for the moment the default arguments are set VERY
        # generously to test the selections...
        
        # user input
        self.pathJoined  = pathJoined[:]
        self.Verbose = Verbose

        # Internal variables
        self.tJoined = Table()

        # String stems for the various pieces
        self.stemCrowd = 'Crowding_to_Precision_0.05_'
        self.stemMSTO = 'msto_'
        self.colPM = 'properMotion_i'

        # magnitiude minus msto
        self.stemCompare = 'mMSTO_'

        # booleans for assessment
        self.stemGood = 'bMag_'
        self.colGoodPM = 'bPM'
        
        # minimum apparent magnitude excess to accept
        self.magSurplus = magSurplus
        self.pmMax = pmMax

        # for spatial selection
        self.spatial_bMin = -20.
        self.spatial_bMax = spatial_bMax
        
        # list of filters to assess
        self.lFilts = lFilts[:]

        # table for results
        self.tRes = Table()
        self.colNames = {}  # column name mapping
        self.colsJudge = {} # columns for judgement
        
        # badval for msto, extinction
        self.badval = -0.1

        # results table column for master boolean: which fields passed?
        self.colMaster = 'passed'

        # OUR RESULTS:
        # standard deviation (in l, b) of the fields that passed the
        # criteria, and the number of "good" fields. 
        self.s_l = 0.
        self.s_b = 0.
        self.nGood = 0
        self.fomRaw = 0
        self.fomIdeal = 1.
        self.fom = 0.

        # temp dump fits name
        self.outFits='tmp_fom.fits'
        
        # by default, run on initialization
        self.setupColnames()
        self.loadJoined()
        self.trimSpatially()

        if len(self.tJoined) < 1:
            if self.Verbose:
                print("fomCalc.init WARN - no table rows.")
            return

        self.initTableResults()
        self.checkForBadvals()

        # perform the comparison
        self.compareMags()
        self.assessFields()
        self.calcFieldDistns()

        # calculate the figure of merit
        self.calcFom()

        # write the results to file
        self.writeResults()
        
    def loadJoined(self):

        """Loads the MAF evaluations"""

        # One method per source because that way the error message
        # will be a bit more informative.
        
        if not os.access(self.pathJoined, os.R_OK):
            if self.Verbose:
                print("fomCalc.loadJoined WARN - cannot read path %s" \
                      % (self.pathJoined))
            return

        self.tJoined = Table.read(self.pathJoined)

    def trimSpatially(self):

        """Trims the region under consideration by spatial area"""

        # For the moment, this is a simple cutoff with limits. We
        # could go more sophisticated later.
        if not 'b' in self.tJoined.colnames:
            return
        goodB = (self.tJoined['b'] > self.spatial_bMin) & \
                (self.tJoined['b'] < self.spatial_bMax)

        # just trim the table for the moment
        self.tJoined = self.tJoined[goodB]

        self.tJoined.meta['spatial_bMin'] = self.spatial_bMin
        self.tJoined.meta['spatial_bMax'] = self.spatial_bMax
        
    def initTableResults(self):

        """Initializes the results table, copying across only the coordinates"""

        self.tRes = Table()
        for sCol in ['l', 'b', 'ra', 'dec']:
            if sCol in self.tJoined.colnames:
                self.tRes[sCol] = self.tJoined[sCol]

        self.tRes.meta = self.tJoined.meta.copy() 
                
    def setupColnames(self):

        """Utility - set up the column names for the comparison quantities for
each filter"""

        # The intention here is to wrap all the column names for
        # magnitudes into a single dictionary, whose properties we
        # then inherit later on.
        
        self.colNames = {}

        self.colsAssess = {}

        for sFilt in self.lFilts:
            colCrowd = '%s%s' % (self.stemCrowd, sFilt)
            colMSTO  = '%s%s' % (self.stemMSTO, sFilt) 

            # column name for the msto magnitude comparison
            colComp = '%s%s' % (self.stemCompare, sFilt)
            
            self.colNames[colComp] = {'crowd':colCrowd, 'msto':colMSTO}

            # column name for the judgement
            self.colsJudge[colComp] = '%s%s' % (self.stemGood, sFilt)

    def checkForBadvals(self):

        """Checks the table for bad values"""

        self.tRes['valid'] = np.repeat(True, len(self.tRes))
        for sFilt in self.colNames.keys():
            for sTyp in self.colNames[sFilt].keys():
                sCol = self.colNames[sFilt][sTyp]
                bValid = np.isfinite(self.tJoined[sCol]) & \
                         (self.tJoined[sCol] > self.badval)

            self.tRes['valid'] = (self.tRes['valid']) & bValid
                

    def compareMags(self):

        """Does the comparison by apparent magnitude"""

        for sComp in self.colNames.keys():
            colCrowd = self.colNames[sComp]['crowd']
            colMSTO = self.colNames[sComp]['msto']

            self.tRes[sComp] = self.tJoined[colCrowd] - self.tJoined[colMSTO]

    def assessFields(self):

        """Having computed the various columns, assess them by whether they
meet the selection criteria"""

        # boolean for matching ALL the magnitude criteria
        bAll = np.copy(self.tRes['valid'])

        # loop through the apparent magnitudes...
        for sMag in self.colsJudge.keys():
            bThis = (self.tRes[sMag] > self.magSurplus)
            colThis = self.colsJudge[sMag]

            self.tRes[colThis] = bThis

            # perform the compound conditional here
            bAll = (bAll) & (bThis)

        # ... and then do the proper motion assessment
        bPropm = self.tJoined[self.colPM] < self.pmMax
        self.tRes[self.colGoodPM] = bPropm # record in table

        # perform the compound conditional
        bAll = (bAll) & (bPropm)
        
        self.tRes[self.colMaster] = bAll

    def calcFieldDistns(self):

        """Simple standard deviation of the fields that passed our criteria"""

        if np.sum(self.tRes[self.colMaster]) < 3:
            if self.Verbose:
                print("fomCalc.calcFieldDistns WARN - <3 fields passed")
            return
        
        bPassed = self.tRes[self.colMaster]
        self.s_l = np.std(self.tRes[bPassed]['l'])
        self.s_b = np.std(self.tRes[bPassed]['b'])
        self.nGood = np.sum(bPassed)

    def calcFom(self):

        """Calculates the raw, ideal figures of merit and finds the ratio"""

        self.calcFomRaw()
        self.calcFomIdeal()
        self.fom = self.fomRaw / self.fomIdeal
        
    def calcFomRaw(self):

        """Calculates the raw figure of merit"""
        
        self.fomRaw = np.float(self.nGood) * self.s_l * self.s_b
        
    def calcFomIdeal(self):

        """Calculates the 'ideal' figure of merit"""

        # for the moment, we assume that every field passes.
        self.fomIdeal = np.float(len(self.tJoined)) * \
                        np.std(self.tRes['l']) * \
                        np.std(self.tRes['b'])

    def writeResults(self):

        """Dumps the results to disk"""

        if len(self.tRes) < 1:
            return

        if len(self.outFits) < 4:
            return
        
        self.tRes.write(self.outFits, overwrite=True)
        
# =====

def testFindFom(magSurplus=0, pmMax=0.5, \
                    pathJoined='TEST_interp_joined.fits', \
                    summFil='DUMMY.lis'):

    """Tests the various stages of the comparison"""

    if not os.access(pathJoined, os.R_OK):
        print("calcFOM.testFindFom WARN - path not readable: %s" \
                  % (pathJoined))

        return

    # prepare the summary file name out of the input file name
    summFil = 'SUMMARY_%s.txt' % (pathJoined.split('.fits')[0])

    fC = fomCalc(pathJoined, magSurplus=magSurplus, pmMax=pmMax)

    # for the moment, let's just output
    #fC.tRes.write('TEST_compoMags.fits', overwrite=True)

    if len(fC.tJoined) < 1:
        return

    # dump the summary statistics to disk for now
    with open(summFil, 'w') as wObj:
        wObj.write('# Summary statistics: %s\n' % (pathJoined))
        wObj.write('# magSurplus=%.2f, pmMax=%.2f\n' \
                       % (magSurplus, pmMax))
        wObj.write('#\n')
        wObj.write("INFO: RAW: %.2f, %.2f, %i, %.2e \n" \
          % (fC.s_l, fC.s_b, fC.nGood, fC.fomRaw))
        wObj.write("INFO: IDEAL: %.2f, %.2f, %i, %.2e\n" \
                       % (np.std(fC.tRes['l']), np.std(fC.tRes['b']), \
                              len(fC.tJoined), fC.fomIdeal))

        wObj.write("INFO: RATIO: %.2e" % (fC.fom))
    
        

    print("INFO: RAW: %.2f, %.2f, %i, %.2e" \
          % (fC.s_l, fC.s_b, fC.nGood, fC.fomRaw))

    print("INFO: IDEAL: %.2f, %.2f, %i, %.2e" \
          % (np.std(fC.tRes['l']), np.std(fC.tRes['b']), \
             len(fC.tJoined), fC.fomIdeal))

    print("INFO: RATIO: %.2e" % (fC.fom))

  
