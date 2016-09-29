import Moog960
import SpectralTools
import glob
import AstroUtils
import sys

"""
This program should be run second in the tutorial sequence.
Once the generateGrid.py finishes creating your grid, you 
can then call this program to read in the raw models
and convolve them with your favorite vsini and instrumental 
resolving power to produce a grid of convolved spectra.  
While you only have to run the generateGrid.py program once
to produce the grid of raw data, you can run convolveGrid.py 
as many times as you like to produce synthetic spectra with
different vsini and instrumental resolving powers.

Instructions:

1) Edit the configuration file to reflect the location/naming convention
of the grid of raw data produced by generateGrid, as well as the vsini,
instrumental resolving power, and output name of the your desired convovled
spectra.  For this tutorial, the name of this configuration file is 

convolveGrid.cfg

2) run convolveGrid.py convolveGrid.cfg

The program will then read in, one by one, all of the raw models from 
grid_data_dir, convolve them with a rotational broadening kernel and 
instrumental resolving power, and save the result in the convolved_out_dir
directory, with the convolved_out_basename as the basename for the file.

"""

config = AstroUtils.parse_config(sys.argv[1])
vsini = config["vsini"]
R = config["resolving_power"]
datadir = config["grid_data_dir"]

outputdir = config["convolved_out_dir"]
outputbase = config["convolved_out_basename"]

gridFiles = glob.glob(datadir+'*_raw.fits')

for filename in gridFiles:
    print filename
    
    Score = Moog960.Score()
    syntheticMelody = Moog960.SyntheticMelody(filename=filename, Score=Score)
    syntheticMelody.selectPhrases(selectAll=True)
    syntheticMelody.rehearse(vsini=vsini, R=R)
    print("%d" % len(syntheticMelody.phrases[0].convolvedData[0].wl))
    convolved_Labels = Score.getLabels(keySignature="CONVOLVED")
    for label in convolved_Labels:
        label.Melody.record(labels=[label], outdir=outputdir, basename=outputbase)

    del(Score)
    del(syntheticMelody)
    del(convolved_Labels)




