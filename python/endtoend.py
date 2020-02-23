#
# endtoend.py
#

# Started 2020-02-23 WIC - wrapper to perform the end-to-end
# calculation of the figure of merit

import os

# our ingredients
import fomStatic
import mapRead
import calcFOM

def go(dbFil='baseline_v1.4_10yrs.db', nside=128, \
           nightMaxCrowd=365, nightMaxPropm=1e4, \
           filtersCrowd=['g', 'r', 'i', 'z', 'y'], \
           pathMSTO='lb_MSTO_ugriz.fits'):

    """End-to-end run of bulge figure of merit"""

    # At each step in the process, we 'gracefully' fail out if the
    # output path is not readable.

    # 1. Evaluate the LSST MAF crowding and proper motion metrics
    pathMAF = fomStatic.TestFewMetrics(\
        dbFil, nside, nightMaxCrowd, nightMaxPropm, \
            filtersCrowd=filtersCrowd)

    if not os.access(pathMAF, os.R_OK):
        print("endtoend.go FATAL - joined-MAF path not readable: %s" \
                  % (pathMAF))
        return

    # 2. Interpolate the LSST MAF metrics to the MSTO spatial locations
    pathInterpol = mapRead.TestInterpMAF(pathMSTO, pathMAF, nneib=False)

    if not os.access(pathInterpol, os.R_OK):
        print("endtoend.go FATAL - interpolated MAF path not readable: %s" \
                  % (pathInterpol))

        return

    # 3. Now that the two tables are merged, compute the figure of merit
    calcFOM.testFindFom(-4, 0.8, pathInterpol)
    
