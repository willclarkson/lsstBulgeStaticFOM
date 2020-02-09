#
# fomStatic.py 
#

# WIC 2020-02-09

# Set of metrics to implement the Gonzalez et al. Bulge static science metric.

# Pseudo:
#
# For each filter, the 5-sigma depth is computed. The i-band PM metric
# is computed. Then the output is sent to a 2D table for comparison to
# an external table (e.g. a new reddening map).

# NOTE - the computation of each data selection separately is handled
# here in python rather than within sims_maf to be more transparent.

# (Note: I think the bundleList idiom actually allows different
# selectors to be used for each bundle, which might be another way to
# implement a diverse set of metrics in a single bundle.)

# 2020-02-09 - test this only with the M5Metric since I don't have the
# crowding metric on my laptop atm...

# Saved for later: https://github.com/LSST-nonproject/sims_maf_contrib/blob/master/tutorials/Spatial_Coordinates.ipynb

import os

# bring in the LSST pieces
import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.maps as maps

class singleMetric(object):

    """Single metric (might be a metric bundle)"""

    def __init__(self, dbFil='sourceDatabase.db', filters=['u'], \
                     dayFirst=-1, dayLast=10000, \
                     NSIDE=32, \
                     metrics=[metrics.CountMetric(col='observationStartMJD')], \
                     dirOut='testMetric', \
                     Verbose=True):

        # control variables
        self.Verbose = Verbose

        # filter for selection
        self.filters = filters[:]
        
        # time selection criteria
        self.dayFirst = dayFirst
        self.dayLast = dayLast

        # nside for the healpix
        self.nside = NSIDE

        # source database
        self.dbFil = dbFil[:]

        # output directory
        self.dirOut = dirOut[:]

        # operatives
        self.canRun = True

        # set of metrics
        self.metrics=metrics[:]

        # MAF quantities common to all the metrics
        self.bundleList = []
        self.slicer=None
        self.sql='DUMMY'
        

        # operations to run on initialization
        self.checkDbReadable()
        self.buildSelString()

    def checkDbReadable(self):

        """Checks that the database is readable"""
        
        if not os.access(self.dbFil, os.R_OK):
            self.canRun = False
            if self.Verbose:
                print("singleMetric.checkDbReadable WARN - database not readable: %s" % (self.dbFil))
           

    def buildSelString(self):

        """Builds the selection string for the metric"""

        self.sql='night > %i and night < %i' % \
            (self.dayFirst, self.dayLast)

        # add on the filter string
        self.sql = '%s and %s' % \
            (self.sql, self.getStringFilters())

    def getStringFilters(self):

        """Utility - returns selection string for filters"""

        strFilter = ''
        if not isinstance(self.filters,list):
            return 'filter="%s"' % (self.filters)

        strFilter = 'filter="%s"' % (self.filters[0])
        if len(self.filters) < 2:
            return strFilter
                            
        strFilter = '(%s' % (strFilter)
        for iFilt in range(1, len(self.filters)):
            strFilter='%s and filter="%s"' % \
                (strFilter, self.filters[iFilt])
        
        strFilter = '%s)' % (strFilter)

        return strFilter


    def setupBundleDict(self):

        """Sets up the metric bundle"""

        # self.slicer = slicers.HealpixSlicer(nside=self.nside, useCache=False)
        self.slicer = slicers.HealpixSlicer(latCol='galb', lonCol='gall', latLonDeg=False, nside=self.nside)


        self.bundleList=[]
        for iMetric in range(len(metrics)):
            thisMetric = None
            thisMetric = self.metrics[iMetric]
            thisBundle = metricBundles.metricBundle(\
                metric, self.slicer, self.sql)
        
            self.bundleList.append(thisBundle)

        self.bundleDict = metricBundles.makeBundlesDict(self.bundleList)

    def setupGroupAndRun(self):

        """Sets up and runs the bundle group"""

        self.checkDbReadable()
        if not self.canRun:
            print("setupGroupAndRun WARN - failed pre-run check. Not running.")
            return
        
        # set up the database connection, ensure the output
        # destination exists
        opsdb = self.db.OpsimDatabase(self.db)
        self.ensureOutdirExists()        
        self.resultsDb = db.ResultsDb(outDir=self.dirOut)

        bgroup = metricBundles.MetricBundleGroup(\
            self.bundleDict, opsdb, \
                outDir=self.dirOut, resultsDb=self.resultsDb)
       
        bgroup.runAll()
 
    def ensureOutdirExists(self):

        """Utility: ensures output directory exists"""
        
        if not os.access(self.outDir, os.R_OK):
            os.makedirs(self.outDir)

# =====

def TestSel(filtr='u'):

    sM = singleMetric(filters=filtr)


    print sM.sql
