import Moog960
import SpectralTools
import astropy.io.fits as pyfits
import sys
import glob
import AstroUtils
import numpy

#This program unpacks a Moog960-compatible .fits file into wavelengths, fluxes, and 
#signal-to-noise (for observed spectra), and saves the result in an ascii file.

config = AstroUtils.parse_config(sys.argv[1])

filenames = numpy.array([x.strip() for x in config["desiredFiles"].split(',')])

if 'outputDir' in config.keys():
    outdir = config["outputDir"]
else:
    outdir = './PlainFITS/'

#Loads in desired file
Score = Moog960.Score(files = filenames)
labels = Score.getLabels(keySignature='CONVOLVED')

num = 0
for label in labels:
    print asdf
    label.Spectrum.savePlainFits(outfileName=outdir+label.Spectrum.filename.split('/')[-1])
    num += 1
