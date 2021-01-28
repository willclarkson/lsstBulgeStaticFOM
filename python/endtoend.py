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

def runSeveral(nside=128, nMax=3, sciserver=True, \
                   lSims=['baseline_2snaps_v1.5_10yrs.db', \
                              'bulges_cadence_bulge_wfd_v1.5_10yrs.db'], \
                   filSims='NONE', fileHasFullPaths=True):

    """Convenience-wrapper to run on several opsims. The
    sciserver=True sets the opsim path appropriately. 

    To specify input lists: either specify a list of filenames as the
    input argument lSims, or populate filSims with one line per
    database. If the file has been given, its contents supersede the
    contents of lSims."""

    ## 2021-01-28 update: set the root directory rather than copying
    ## the opsims in here.
    if sciserver:
        dbroot = '/home/idies/workspace/lsst_cadence/FBS_1.5/'
    else:
        dbroot = '/sims_maf/fbs_1.5/short_exp/'

    # 2021-01-28 WIC: There are so many opsims that it's not
    # immediately clear to me how pass the inputs in correctly. I
    # suspect a file is going to be the way forward. For the moment,
    # though, we will also set a default entry.

    # by default, runs on all the .db files in the current working
    # directory
    # lDb = glob.glob("*1.4*.db")

    # read the file list from file if given. This list may well have
    # leading paths.
    if os.access(filSims, os.R_OK):
        with open(filSims, 'r') as rObj:
            lSims = rObj.read().splitlines()

        # If the file already lists the full paths to the opsims, then
        # we don't need to add it to the individual paths.
        if fileHasFullPaths:
            dbroot = ''

    # how many do we run?
    if nMax < 0:
        iMax = len(lSims)
    else:
        iMax = min([nMax, len(lSims)])

    for iDb in range(iMax):
        thisDb = "%s%s" % (dbroot, lSims[iDb])
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
