#
# fomStatic.py 
#

# WIC 2020-02-09

# Set of metrics to implement the Gonzalez et al. Bulge static science metric.

# Pseudo:
#
# For each filter, the 5-sigma depth is computed (the standard
# "CrowdingM5Metric" is used). The i-band PM metric is computed. The
# result is one table per metric. These tables are then joined into a
# single 2D table for comparison with external tables (like a map of
# MSTO apparent magnitude predictions.)

# NOTE - the computation of each data selection separately is handled
# here in python rather than within sims_maf to be more transparent.

# (Note: I think the bundleList idiom actually allows different
# selectors to be used for each bundle, which might be another way to
# implement a diverse set of metrics in a single bundle.)


import os, glob
import numpy as np
from astropy.table import Table, join
from astropy import units as u

# bring in the LSST pieces
import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.maps as maps
import lsst.sims.maf.plots as plots

# UPDATE 2020-08-13: when running on sciserver, the updated crowding
# metrics are available locally. Handle that here.
blah = None
try:
    blah = metrics.NstarsMetric
except:
    try:
        import crowdingMetric
        blah = crowdingMetric.NstarsMetric
    except:
        NO_NSTARS = True
    

class singleMetric(object):

    """Single metric (might be a metric bundle)"""

    def __init__(self, dbFil='baseline_v1.4_10yrs.db', filters=['r'], \
                     dayFirst=-1, nightMax=10000, \
                     NSIDE=64, \
                 #metrics=[metrics.CountMetric(col='observationStartMJD')], \
                 metrics=[metrics.CrowdingM5Metric(crowding_error=0.05, filtername='r', maps=['TrilegalDensityMap']) ], \
                 dirOut='testMetric', \
                     Verbose=True, \
                 getFilterFromMetric=True, kwargs={}):

        # WATCHOUT - this trusts the user to input sensible arguments
        # for the metric. If changing the filter for the crowding map,
        # need to change BOTH the filters and the filter name argument
        # in the metrics value. Return to this later!
        
        # control variables
        self.Verbose = Verbose

        # input filter (could be superseded below for consistency)
        self.filters = filters[:]
        
        # try to get the filter from the zeroth metric?
        if getFilterFromMetric:
            if hasattr(metrics[0],'filtername'):
                foundFilter = metrics[0].filtername
                if self.Verbose:
                    print("INFO: using filter passed to metric[0]: %s" % (foundFilter))
                self.filters = foundFilter[:]

        # time selection criteria
        self.dayFirst = dayFirst
        self.nightMax = nightMax

        # nside for the healpix
        self.nside = NSIDE

        # source database
        self.dbFil = dbFil[:]

        # output directory for temp files
        self.dirOut = dirOut[:]

        # list of output paths
        self.pathsOut = []

        # operatives
        self.canRun = True

        # set of metrics
        self.metrics=metrics[:]

        # MAF quantities common to all the metrics
        self.bundleList = []
        self.slicer=None
        self.sql='DUMMY'
        
        # additional keyword arguments to pass to the metric
        # 2020-08-13 CURRENTLY UNUSED
        self.kwargs = kwargs
        
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
            (self.dayFirst, self.nightMax)

        # add on the filter string
        self.sql = '%s and %s' % \
            (self.sql, self.getStringFilters())

    def pithyFilterString(self):

        """Generates single string with all the filters"""

        if isinstance(self.filters, str):
            return self.filters[:]

        # If we've been given a list of filter strings, return joined.
        return ''.join(self.filters)
        
    def getStringFilters(self):

        """Utility - returns selection string for filters."""

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

        # 2020-02-09 - slice by spatial coordinates
        self.slicer = slicers.HealpixSlicer(nside=self.nside, useCache=False)
        ##self.slicer = slicers.HealpixSlicer(latCol='galb', lonCol='gall', \
        ##                                    latLonDeg=False, nside=self.nside, useCache=False)
        #self.slicer = slicers.HealpixSlicer(latCol='ditheredDec', lonCol='ditheredRA', latLonDeg=False, nside=self.nside, useCache=False)
        
        # (For reference, here is the spatial slicer MAF tutorial:)
        # https://github.com/LSST-nonproject/sims_maf_contrib/blob/master/tutorials/Spatial_Coordinates.ipynb

        # sets the map for a given filter
        filtername = self.pithyFilterString()
        mapsList = [maps.TrilegalDensityMap(filtername=filtername, nside=self.nside)]

        # 2020-02-11 cargo-cult copy to try to get this to plot galL, galB
        plotFuncs = [plots.HealpixSkyMap()]
        
        self.bundleList=[]
        self.outNPZlist=[]
        self.outBundNames=[]
        for iMetric in range(len(self.metrics)):
            thisMetric = None
            thisMetric = self.metrics[iMetric]
            
            # find or generate names we'll need later                                                                                    
            thisName=thisMetric.name.replace(" ","_")
            thisName = '%s_%s' % (thisName, filtername)

            # If we're not doing the crowding metric, we don't need the maps. This is maybe a bit clunky...

            # 2020-08-13 updated for nstarsmetric
            mapsListUse = mapsList[:]
            if thisName.find('rowd') < 0 and thisName.find('NstarsMetric') < 0:
                mapsListUse = []

            thisBundle=metricBundles.MetricBundle(thisMetric,self.slicer,self.sql, \
                                                  mapsList=mapsListUse, plotFuncs=plotFuncs)
                
                
            #print("INFO: this Metric Name:",thisName)
            #print("INFO:", thisBundle.fileRoot)
            
            self.outNPZlist.append('%s/%s.npz' % \
                                       (self.dirOut, thisBundle.fileRoot))
            self.outBundNames.append(thisName)
            # generate output file path for the .npz file out of this
            
            
            self.bundleList.append(thisBundle)

        self.bundleDict = metricBundles.makeBundlesDictFromList(self.bundleList)

    def setupGroupAndRun(self):

        """Sets up and runs the bundle group"""

        self.checkDbReadable()
        if not self.canRun:
            print("setupGroupAndRun WARN - failed pre-run check. Not running.")
            return
        
        # set up the database connection, ensure the output
        # destination exists
        opsdb = db.OpsimDatabase(self.dbFil)
        self.ensureOutdirExists()        
        self.resultsDb = db.ResultsDb(outDir=self.dirOut)

        self.bgroup = metricBundles.MetricBundleGroup(\
            self.bundleDict, opsdb, \
                outDir=self.dirOut, resultsDb=self.resultsDb)
       
        self.bgroup.runAll()
        self.bgroup.plotAll()
        
    def ensureOutdirExists(self):

        """Utility: ensures output directory exists"""

        if len(self.dirOut) < 3:
            self.dirOut='testOut'
            if self.Verbose:
                print("ensureOutdirExists INFO - refusing to use <3 char output directory.")
                print("ensureOutdirExists INFO - defaulted to %s" % (self.dirOut))
        
        if not os.access(self.dirOut, os.R_OK):
            os.makedirs(self.dirOut)

    def translateResultsToArrays(self):

        """Utility - translates the .npz output to a flatter format"""

        print("translateResults INFO - output %s" % (self.resultsDb))

        print(self.outNPZlist)

        self.loadMetricValues(self.outNPZlist[0], "test_%s.fits" % (self.outBundNames[0]), \
                              self.outBundNames[0])
        
    def loadMetricValues(self, pathIn='BLAH', tableName='TEST.fits', metricName='metricValues'):

        """Loads the metric values into memory"""

        if not os.access(pathIn, os.R_OK):
            if self.Verbose:
                print("loadMetricValues WARN - cannot read input path %s" % (pathIn))
            return

        metricValues = None
        slicePoints = None
        sliceMask = None
        
        with np.load(pathIn) as resNPZ:
            metricValues = resNPZ.f.metricValues
            slicePoints = resNPZ.f.slicePoints.item()
            sliceMask = resNPZ.f.mask
            
        # I'm open to suggestions for how to do this more
        # nicely... For the moment, let's use an astropy table.
        tRes = Table()

        # 2020-02-09 WIC - these 'ra' and 'dec' might actually already be [l,b]...
        
        tRes['sid'] = slicePoints['sid']
        tRes['ra'] = np.degrees(slicePoints['ra'])
        tRes['dec'] = np.degrees(slicePoints['dec'])
        # tRes['metricValues'] = metricValues
        tRes[metricName] = metricValues
        
        # I don't trust the 'mask' value. make our own boolean
        # variable. Make our own
        tRes['%s_gtr0' % (metricName)] = metricValues > 1e-5
        tRes['%s_finite' % (metricName)] = np.isfinite(metricValues)
        
        # set units as radians
        tRes['ra'].units = u.deg
        tRes['dec'].units = u.deg
        
        # try writing this to disk
        
        pathRes = '%s/%s' % (self.dirOut, tableName)

        tRes.write(pathRes, format='fits', overwrite=True)
        
        # append this to the list of output paths
        self.pathsOut.append(pathRes)

# =====

def TestSel(filtr='r', nside=64):

    """Test the runthru metric"""

    # 2020-02-10 - updated with Peter Yoachim's help: sets the crowding metric for the given filter.
    # mapsList = [maps.TrilegalDensityMap(filtername=filtr, nside=nside)]
    metric = metrics.CrowdingM5Metric(crowding_error=0.05, filtername=filtr, maps=[])
    
    #metric = metrics.CrowdingM5Metric(crowding_error=0.05,filtername='%s' % (filtr))

    sM = singleMetric(metrics=[metric])
    # print(sM.sql)

    sM.setupBundleDict()
    sM.setupGroupAndRun()
    sM.translateResultsToArrays()


def TestFewMetrics(dbFil='baseline_v1.4_10yrs.db', nside=128, \
                       nightMaxCrowd=365, nightMaxPropm=1e4, \
                       filtersCrowd = ['g', 'r', 'i', 'z', 'y'], \
                       cleanTmpDir=True, tmpDir='./tmpMetrics', \
                       buildPathJoined = True, \
                       dirOut='tmpProds', \
                       Verbose=True, \
                   crowdingUncty=0.05):

    """Test routine to test a few metrics with different
    selections. Returns the path to the joined metric file. 

    dbfil -- OPSIM database file to assess metrics

    nside -- HEALPIX resolution (2020-02-23: nside > 128 nonstandard in MAF)

    nightMaxCrowd -- maximum obsdate for the crowding metrics

    nighMaxPropm -- maximum obsdate for the proper motion metric

    filtersCrowd -- list of filters for which to perform the crowding checks

    cleanTmpDir -- remove intermediate results

    tmpDir -- directory to hold intermediate products

    buildPathJoined -- construct path for joined-metrics from the
    database filename

    dirOut -- directory for output files

    crowdingUncty -- crowding error

    Verbose -- provide 'informative' terminal output"""

    # Path to return (to the joined MAF outcomes. One for success, one
    # for failure.) Do the path building up-front so that it can be
    # easily tested without waiting for evaluation of everything
    # else.
    if not os.access(dirOut, os.R_OK) and len(dirOut) > 2:
        os.makedirs(dirOut)

    pathFail = 'NONE'
    pathSuccess = '%s/tmp_joined.fits' % (dirOut)
    if buildPathJoined:
        pathSuccess = '%s/METRICS_%s.fits' % (dirOut, dbFil.split('.db')[0])

    # Is the opsim database file accessible?
    if not os.access(dbFil, os.R_OK):
        print("fomStatic.TestFewMetrics WARN - dbfile not accessible: %s" \
                  % (dbFil))
        return pathFail

    # use the same dbfile throughout
    # dbFil='baseline_v1.4_10yrs.db'

    # list of single-metrics
    listMetrics = []

    if nside < 64:
        print("fomStatic.TestFewMetrics WARN - maps won't work for nside < 64")
        return pathFail
        
    for filt in filtersCrowd:
#        metricThis = metrics.CrowdingM5Metric(crowding_error=0.05, \
#                                                  filtername=filt)

        metricThis = metrics.CrowdingM5Metric(crowding_error=crowdingUncty, \
                                              filtername=filt, maps=[])

        
        sM = singleMetric(metrics=[metricThis], \
                              nightMax=nightMaxCrowd, \
                              NSIDE=nside, dbFil=dbFil, \
                              getFilterFromMetric=True, \
                              dirOut=tmpDir[:])
    
        sM.setupBundleDict()
        sM.setupGroupAndRun()
        sM.translateResultsToArrays()

        listMetrics.append(sM)

    # Now do the proper motion metric
    metricPropmI = metrics.ProperMotionMetric()
    filterPropmI = 'i'

    sP = singleMetric(metrics=[metricPropmI], nightMax=nightMaxPropm, \
                          NSIDE=nside, dbFil=dbFil, \
                          getFilterFromMetric=False, \
                          filters=[filterPropmI], \
                          dirOut=tmpDir[:])

    # Now do the NstarsMetric if we have access to it... This is just
    # a little awkward:
    hasNstars = False
    try:
        nstarsMetric = metrics.NstarsMetric
        hasNstars = True
    except:
        try:
            import crowdingMetric
            nstarsMetric = crowdingMetric.NstarsMetric
            hasNstars = True
        except:
            hasNstars = False


    if hasNstars:

        print("fomStatic.TestFewMetrics INFO - trying density metrics")
        
        # set up the nstars metric and send to our horrendous wrapper
        densMetricCrowd = nstarsMetric(maps=['TrilegalDensityMap'], \
                                       crowding_error=crowdingUncty, \
                                       ignore_crowding=False)

        densMetricNoCrowd = nstarsMetric(maps=['TrilegalDensityMap'], \
                                         crowding_error=crowdingUncty, \
                                         ignore_crowding=True)

        # 2020-08-13 we set up a separate directory for the
        # non-crowded output just so we can see if we're otherwise
        # overwriting our output directory
        dirOutDensCrowd = '%s_crowd' % (sM.dirOut[:])
        dirOutDensNoCrowd = '%s_noCrowd' % (sM.dirOut[:])


        # We'll write this out for the moment...
        
        SNcrowd = singleMetric(metrics=[densMetricCrowd], nightMax=nightMaxPropm, \
                               NSIDE=nside, dbFil=dbfil, \
                               getFilterFromMetric=False, \
                               filters=['r'], \
                               dirOut=dirOutDensCrowd)

        SNnocrowd = singleMetric(metrics=[densMetricNoCrowd], nightMax=nightMaxPropm, \
                               NSIDE=nside, dbFil=dbfil, \
                               getFilterFromMetric=False, \
                               filters=['r'], \
                                 dirOut=dirOutDensNoCrowd)

        for dens in [SNcrowd, SNnocrowd]:
            dens.setupBundleDict()
            dens.setupGroupAndRun()
            dens.translateResultsToArrays()
            listMetrics.append(dens)
        
        
    # ensure the same output tmp directory is used as for the sM
    # object
    sP.dirOut = sM.dirOut[:]

    sP.setupBundleDict()
    sP.setupGroupAndRun()
    sP.translateResultsToArrays()

    listMetrics.append(sP)

    # list of output paths for results
    listOutPaths = []
    for metricObj in listMetrics:
        listOutPaths.append(metricObj.pathsOut[0])

    # remove any old instances of pathSuccess - we don't want to
    # accidentally pick up a prior analysis here.
    if os.access(pathSuccess, os.R_OK):
        os.remove(pathSuccess)
        
    # now join the tables together
    mergeTables(listOutPaths, pathJoined=pathSuccess[:])

    # if asked, remove the contents of the temporary directory. To be
    # safe, refuse to do this if the tmpdir has fewer than 3
    # characters (e.g. if the user is trying to erase the contents of
    # "./").
    if cleanTmpDir:

        # since we've specified in this method the tempdir to use, we
        # could just pass the same variable in here. I prefer to use
        # the tmpdir actually used, though...
        dirIntermed = sM.dirOut[:]
        
        # all the housekeeping, safety checks, etc., are ported to a
        # separate method below.
        removeExtraFiles(dirIntermed, ['.fits', 'db', 'npz'], \
                             safetyCheck=True)

    # if we've got here, return the path to the output file
    pathRet = pathFail[:]
    if os.access(pathSuccess, os.R_OK):
        pathRet = pathSuccess[:]

    return pathRet

def mergeTables(paths=[], pathJoined='./TEST_joined.fits'):

    """Merges tables specified one per file. Written for use on tables
    built from exactly the same spatial slice."""
        
    tJoined = Table()
    if len(paths) < 1:
        return

    # read the first one in to start the joined table
    tJoined = Table.read(paths[0])
    for iTable in range(1, len(paths)):
        tThis = Table.read(paths[iTable])
        tJoined = join(tJoined, tThis)
        

    # write the joined table to disk
    tJoined.write(pathJoined, overwrite=True)
        
def removeExtraFiles(dirTop='DUMMY', lTails = ['.fits', '.db'], \
                         safetyCheck=True):

    """From dirTop, removes files that end in lTails entries. 

    safetyCheck = perform some conservative safety checks (e.g. don't
    do anything if the directory to check has fewer than 4
    characters (like ./, ../) )"""

    print("INFO - cleaning intermediate products")
    print(dirTop, os.access(dirTop, os.W_OK))

    # directory must actually be accessible
    if not os.access(dirTop, os.W_OK):
        return

    # safety check?
    if safetyCheck:
        if len(dirTop) < 4:
            return

    # lTails must actually be a list. We don't want to iterate through
    # all the characters in a string.
    if not isinstance(lTails, list):
        return

    for sTail in lTails:
        lJunk = glob.glob('%s/*%s' % (dirTop, sTail))
        if len(lJunk) < 1:
            continue

        for thisPath in lJunk:
            os.remove(thisPath)
        
