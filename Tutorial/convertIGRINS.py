import Moog960
import astropy.io.fits as pyfits


#This program reads in an IGRINS spectrum from the .fits files produced by the
#reduction pipeline, and saves it in a Moog960-readable .fits file.  It merges
#the spectral orders, weighting by signal-to-noise ratio.
