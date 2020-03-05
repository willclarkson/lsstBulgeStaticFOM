# python #

Contents:

* **endtoend.py** -- end-to-end wrapper. It calls, in turn:

* **fomStatic.py** -- uses MAF to compute the crowding and proper motion metrics.

* **mapRead.py** -- uses HEALPIX to interpolate the metrics onto the (Galactic) coordinates of the (fits-format) MSTO map. Produces 
a FITS file with the metrics and MSTO magnitudes side-by-side for each location in the MSTO map.

* **calcFOM.py** -- Performs the calculation of the figure of merit using the fused output from **mapRead.py**
