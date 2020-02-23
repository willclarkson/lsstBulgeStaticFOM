#
# mapRead.py
# 
# Started 2020-02-16 WIC - use HEALPIX to (rapidly??) interpolate one
# map onto another. Application in mind: matching two coverage maps
# for the bulge-static LSST figure of merit.

# 2020-02-16 - using the "standard" healpy rather than the new astropy
# healpy because I think the former is the one that's in the LSST stack.
import healpy as hp

from astropy.table import Table
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import os

class pointSet(object):

    """Sets of points"""

    # Try to be as agnostic as possible about the input format. It
    # must have RA, DEC or it must be straightforward to produce RA,
    # DEC from the coordinates we do have.

    def __init__(self, tMap=Table(), \
                     colRA='ra', colDEC='dec', \
                     Verbose=True, \
                     assignHealpix=True, \
                     colGL = 'l', colGB='b', \
                     nside=64, guessNSIDE=True, \
                     isDeg=True, guessDeg=True, \
                     nested=True, \
                     colHEALPIX='healpix', addHEALPIXtoTable=True):

        # Control variables
        self.Verbose = Verbose
        self.assignHealpix = assignHealpix
        
        # if the tMap argument is a string, treat it as a table file
        # and try to load it
        if isinstance(tMap, str):
            self.loadTable(tMap[:])
        else:
            self.tMap = tMap

        # Columns for equatorial and Galactic coordinates
        self.colRA = colRA[:]
        self.colDEC = colDEC[:]

        self.colGL = colGL[:]
        self.colGB = colGB[:]

        # user-passed keyword for degrees or not (needs to come before
        # call to getRAfromGalactic)
        self.isDeg = isDeg

        # find RA, DEC if present, attempt to convert from Galactic
        # longitude, latitude if not.
        self.hasRADEC=False
        self.checkRADEC()
        if not self.hasRADEC:
            self.getRAfromGalactic()

        # By this point we should have our RA, DEC. Determine if
        # degrees.
        self.guessIsDegrees()

        # parse the coords into colatitude and longitude
        # colatitude, longitude.
        self.theta = np.array([])
        self.phi = np.array([])
        self.parseRADECforHP()

        # some housekeeping parameters
        self.healpix = np.array([])

        ## Parameters and action for optional HEALPIX assignment
        ## follow below.
        self.nested=nested
        self.nside = nside # later on, we can guess NSIDE from the
                           # input

        # useful for some applications), what do we call it?
        self.colHEALPIX = colHEALPIX[:]

        # assign healpix?
        if assignHealpix:
            if guessNSIDE:
                self.guessNSIDE()
            self.assignHEALPIXid()
            self.sortByHEALPIX()

    def guessNSIDE(self):

        """Guesses the NSIDE from the size of the input array."""
        
        # (Note that the array needs to cover the whole sky for this
        # to be meaningful). We could add syntax to perform this check
        # later.

        if len(self.tMap) < 1:
            if Verbose:
                print("getHEALPIX WARN - table not populated.")
            return

        nsideGuess = 0
        try:
            nSideGuess = hp.npix2nside(len(self.tMap))
            self.nside = nSideGuess
        except:
            if self.Verbose:
                print("getHEALPIX WARN - table length does not correspond to an NSIDE.")
            return

    def guessIsDegrees(self):

        """Use the bounds on the RA, DEC to guess whether the input
        coordinates are in degrees"""

        self.checkRADEC()
        if not self.hasRADEC:
            return

        maxRA = np.max(self.tMap[self.colRA])
        maxDE = np.max(self.tMap[self.colDEC])

        self.isDeg = False
        if maxRA > 2.01*np.pi or np.maxDE > 0.501*np.pi:
            self.isDeg = True

    def loadTable(self, fileIn='DUMMY'):

        """Convenience-method: load the table from disk"""

        if not os.access(fileIn, os.R_OK):
            if self.Verbose:
                print("loadTable WARN - cannot read path %s" % (fileIn))
            return

        try:
            self.tMap = Table.read(fileIn)
        except:
            if self.Verbose:
                print("loadTable WARN - problem loading table %s" \
                          % (fileIn))

    def checkRADEC(self):

        """Checks that the RA, DEC columns are present"""

        self.hasRADEC = self.colRA in self.tMap.colnames and \
            self.colDEC in self.tMap.colnames

    def assignHEALPIXid(self):

        """Assigns the healpix ID for the coordinates"""

        # (re-)perform the check here.
        self.checkRADEC()
        if not self.hasRADEC:
            return

        self.healpix = hp.ang2pix(self.nside, self.theta, self.phi, \
                                      nest=self.nested)

        # add to the table (must do this, since we may later want to
        # sort by HEALPIX ID)
        self.tMap[self.colHEALPIX] = self.healpix
        self.tMap.meta['HPnested'] = self.nested

    def sortByHEALPIX(self):

        """Sorts the table by HEALPIX ID"""

        # uses astropy's built-in sort in the expectation (hope??)
        # that it will be fast
        self.tMap.sort(self.colHEALPIX)

    def parseRADECforHP(self):
        
        """Older versions of healpy do not have the lonlat
        keyword. This utility routine converts the RA, DEC into
        radians, and converts DEC into the healpy-expected format (0
        at the 'north pole,' pi at the 'south pole')"""

        self.checkRADEC()
        if not self.hasRADEC:
            if self.Verbose:
                print("parseRADECforHP WARN - table does not have RA, DEC columns %s, %s" % (self.colRA, self.colDEC))
            return

        raRad = np.copy(self.tMap[self.colRA])
        deRad = np.copy(self.tMap[self.colDEC])
        if np.max(raRad) > 2.01*np.pi:   # if input is degrees
            raRad = np.radians(raRad)
            deRad = np.radians(deRad)
            
        # convention conversion for the dec: zero at the north pole,
        # +pi at the south pole
        deRad = 0.5*np.pi - deRad

        self.theta = deRad
        self.phi = raRad

    def getRAfromGalactic(self, clobber=False):

        """Populate RA, DEC columns from Galactic longitude and
        latitude columns"""

        # We do this so that we can compute the colatitude and
        # longitude that HEALPY will need.

        # do nothing if the table already has RA, DEC and
        # clobber=False
        self.checkRADEC()
        if self.hasRADEC and not clobber:
            return

        if not self.colGL in self.tMap.colnames or \
                not self.colGB in self.tMap.colnames:
            if self.Verbose:
                print("getRAfromGalactic WARN - columns not present: %s, %s" \
                          % (self.colGL, self.colGB))
            return

        # now we do the conversion. If astropy's functionality turns
        # out to be a little slow, we can recode later.
        unitGalac = u.degree
        if not self.isDeg:
            unitGalac = u.rad
            
        gal = SkyCoord(self.tMap[self.colGL], self.tMap[self.colGB], \
                           unit=(unitGalac, unitGalac), frame='galactic')
        
        # 2020-02-16: since I suspect this might be a little slow to
        # parse, provide input on status
        if self.Verbose:
            print("getRAfromGalactic INFO - parsing coords...")
        
        self.tMap[self.colRA]  = gal.icrs.ra.degree
        self.tMap[self.colDEC] = gal.icrs.dec.degree

        if not self.isDeg:
            self.tMap[self.colRA]  = np.radians(self.tMap[colRA])
            self.tMap[self.colDEC] = np.radians(self.tMap[colDEC]) 

class MapPair(object):

    """Pair of maps with values to interpolate. One is assumed to be a
    regular grid, the other is a set of points which will need
    interpolation."""

    def __init__(self, pathGrid='TEST_joined_radec.fits', \
                     pathPoints='lb_MSTO_ugriz_cutbm2.fits', \
                     pathJoined='TEST_interp_joined.fits', \
                     Verbose=True, \
                     nearestNeighbor=True):

        # control variables
        self.Verbose = Verbose

        # set up the pair of catalog objects. By default, objGrid will
        # attempt to assign a HEALPIX ID to each row.
        self.objGrid = pointSet(pathGrid, Verbose=self.Verbose, \
                                    nested=True)
        self.objPoints = pointSet(pathPoints, assignHealpix=False, \
                                      Verbose=self.Verbose)

        print("mapPair INFO - reading file %s" % (pathPoints))

        # path to write the output file
        self.pathJoined = pathJoined[:]

        # list of table keywords to interpolate from the original grid
        self.cols2interp = ['Crowding_to_Precision_0.05_g'] # default

        self.filters = ['u', 'g', 'r', 'i', 'z', 'y']
        self.strStem = 'Crowding_to_Precision_0.05_'
        self.buildCols2interp()

        # colname string extension used to determine mask for the grid
        self.strFinite = '_finite'
        self.strNonzero = '_gtr0'

        # for debugging, it may be useful to use a nearest-neighbor
        # approach. This flag sets the option.
        self.nearestNeib = nearestNeighbor

    def buildCols2interp(self):

        """Builds a list of columns to interpolate"""

        self.cols2interp = []
        for filtr in self.filters:
            self.cols2interp.append('%s%s' % (self.strStem, filtr))

        self.cols2interp.append('properMotion_i')

    def doInterpolation(self):

        """Interpolates through all the quantities we've asked for"""

        for colInterp in self.cols2interp:
            self.interpOneColumn(colInterp)

    def interpOneColumn(self, colInterp='DUMMY'):

        """Interpolates the positions of the 'grid' object onto the
        'points' object, for a particular keyword"""

        # Notice that the original array is assumed to be sorted in
        # order by healpix. 

        if not colInterp in self.cols2interp:
            return

        if not colInterp in self.objGrid.tMap.colnames:
            return

        
        # set up the data and mask
        hpData = self.objGrid.tMap[colInterp]
        bMask = np.isnan(hpData)

        # construct the masked array in healpix order from the
        # original grid
        colFinite = '%s%s' % (colInterp, self.strFinite)
        colNonzer = '%s%s' % (colInterp, self.strNonzero) 

        for goodCol in [colFinite, colNonzer]:
            if goodCol in self.objGrid.tMap.colnames:
                bMask = (bMask) | (self.objGrid.tMap[goodCol] < 1)

        hpVec = np.ma.masked_array(hpData, bMask)

        # Now try the interpolation
        if not self.nearestNeib:
            interpVal = hp.pixelfunc.get_interp_val(\
                hpVec, self.objPoints.theta, self.objPoints.phi, \
                    nest=self.objGrid.nested)
        else:
            lPix = hp.pixelfunc.ang2pix(\
                self.objGrid.nside, \
                    self.objPoints.theta, self.objPoints.phi,\
                    nest=self.objGrid.nested)
            interpVal = hpVec[lPix]

        # add the interpolated values to the POINTS table (take the
        # unit across as well)
        self.objPoints.tMap[colInterp] = interpVal
        self.objPoints.tMap[colInterp].unit = \
            self.objGrid.tMap[colInterp].unit

        # add metadata keyword to the table
        sMethod = 'bilinear'
        if self.nearestNeib:
            sMethod = 'nearestNeibour'

        self.objPoints.tMap.meta['IntrpTyp'] = sMethod 

    def writeJoined(self):

        """Writes the interpolated table to a new file w/ joined columns"""

        self.objPoints.tMap.write(self.pathJoined, overwrite=True)

def testLoadTable(tablTest='TEST_joined_radec.fits', \
                      tablMSTO = 'lb_MSTO_ugriz_cutbm2.fits'):

    """Tests our framework on an example map"""

    OPSIM = pointSet(tablTest)

    MSTO = pointSet(tablMSTO, assignHealpix=False)

    #print(OPSIM.nside)
    #print(OPSIM.hasRADEC)
    #print(OPSIM.isDeg, np.max(OPSIM.tMap[OPSIM.colRA]))

    print(OPSIM.tMap)

def testPair(nneib=False, fullVVV=False):

    """Tests loading the pair of maps"""

    pathOut = 'TEST_interp_joined.fits'
    if nneib:
        pathOut = 'TEST_interp_joined_nneib.fits'

    pathVVV = 'lb_MSTO_ugriz_cutbm2.fits'
    if fullVVV:
        pathVVV = 'lb_MSTO_ugriz.fits'

    mp = MapPair(pathPoints = pathVVV, \
                     nearestNeighbor=nneib, pathJoined=pathOut)
    mp.doInterpolation()
    mp.writeJoined()

def TestInterpMAF(pathMSTO='lb_MSTO_ugriz.fits', \
                      pathMAF='DUMMY', \
                      nneib=False):

    """Interpolates the MAF evaluations onto the positions in the MSTO
    map Returns the path to the interpolated file.

    nneib -- use nearest-neighbor interpolation? (Default is bilinear)."""

    # Again, build the success/fail paths to return
    pathFail = 'DUMMY'
    pathSuccess = 'MERGED_%s' % (pathMAF)
    if nneib:
        pathSuccess = 'MERGED_NNEIB_%s' % (pathMAF)

    if not os.access(pathMSTO, os.R_OK):
        print("mapRead.TestInterpMAF WARN - cannot read MSTO path %s" \
                  % (pathMSTO))
        return pathFail

    if not os.access(pathMAF, os.R_OK):
        print("mapRead.TestInterpMAF WARN - cannot read MAF path %s" \
                  % (pathMAF))

        return pathFail

    # ensure we don't accidentally pick up an older analysis
    if os.access(pathSuccess, os.W_OK):
        os.remove(pathSuccess)

    # ok now actually do the operations...
    mp = MapPair(pathMSTO=pathMSTO[:], pathMAF=pathMAF[:], \
                     nneib=nneib, \
                     pathJoined = pathSuccess)
    
    mp.doInterpolation()
    mp.writeJoined()

    # return the "success" path if that actually produced the output
    # we need.
    if os.access(pathSuccess, os.R_OK):
        return pathSuccess

    return pathFail
