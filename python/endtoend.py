#
# endtoend.py
#

# Started 2020-02-23 WIC - wrapper to perform the end-to-end
# calculation of the figure of merit

import os, glob

# our ingredients
import fomStatic
import mapRead
import calcFOM

def runSeveral(nside=128, nMax=3):

    """Convenience-wrapper to run on several opsims"""

    # by default, runs on all the .db files in the current working
    # directory
    lDb = glob.glob("*1.4*.db")

    # how many do we run?
    if nMax < 0:
        iMax = len(lDb)
    else:
        iMax = min([nMax, len(lDb)])

    for iDb in range(iMax):
        thisDb = lDb[iDb]
        go(thisDb, nside=nside)

def go(dbFil='baseline_v1.4_10yrs.db', nside=128, \
           nightMaxCrowd=365, nightMaxPropm=1e4, \
           filtersCrowd=['g', 'r', 'i', 'z', 'y'], \
           magSurplus=1., pmMax=0.8, \
           pathMSTO='lb_MSTO_ugriz.fits'):

    """End-to-end run of bulge figure of merit"""

    # At each step in the process, we 'gracefully' fail out if the
    # output path is not readable.

    # 2020-03-12 produce subdirectory for output files
    dirSub = 'nside%i' % (nside)
    if not os.access(dirSub, os.R_OK):
        os.makedirs(dirSub)

    # 1. Evaluate the LSST MAF crowding and proper motion metrics
    pathMAF = fomStatic.TestFewMetrics(\
        dbFil, nside, nightMaxCrowd, nightMaxPropm, \
            filtersCrowd=filtersCrowd, \
            dirOut=dirSub[:])

    if not os.access(pathMAF, os.R_OK):
        print("endtoend.go FATAL - joined-MAF path not readable: %s" \
                  % (pathMAF))
        return

    if not os.access(pathMSTO, os.R_OK):
        print("endtoend.go WARN - MSTO path not readable. Not proceeding further.")
        return

    # 2. Interpolate the LSST MAF metrics to the MSTO spatial locations
    pathInterpol = mapRead.TestInterpMAF(pathMSTO, pathMAF, nneib=False)

    if not os.access(pathInterpol, os.R_OK):
        print("endtoend.go FATAL - interpolated MAF path not readable: %s" \
                  % (pathInterpol))

        return

    # 

    # 3. Now that the two tables are merged, compute the figure of merit
    thisSumm = calcFOM.testFindFom(magSurplus, pmMax, pathInterpol)
    print("DONE: %s" % (thisSumm))
    print("###########################")
