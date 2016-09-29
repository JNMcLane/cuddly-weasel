import Moog960
import AstroUtils
import sys
import numpy

"""
1) This program should be run third in the tutorial sequence.  Once the grid has been
convovled with the previous program, this program reads in the grid of convolved
models and interpolates to give a synthetic spectrum corresponding to the desired
parameters (Teff, log g, and Bfield).  The blended spectrum is then saved in a
.fits file for further processing

1) Edit the configuration parameters in blendModels.cfg to reflect your desired
blends

2) run blendModels.py blendModels.cfg

"""

config = AstroUtils.parse_config(sys.argv[1])

wlStart = config["wlStart"]
wlStop = config["wlStop"]

try:
    desiredTeff = numpy.array(config['desiredTeff'].split(','), dtype=numpy.int)
except:
    desiredTeff = numpy.array([config['desiredTeff']], dtype=numpy.int)
try:
    desiredlogg = numpy.array(config['desiredlogg'].split(','), dtype=numpy.float)
except:
    desiredlogg = numpy.array([config['desiredlogg']], dtype=numpy.float)
try:
    desiredBfield = numpy.array(config['desiredBfield'].split(','), dtype=numpy.float)
except:
    desiredBfield = numpy.array([config['desiredBfield']], dtype=numpy.float)

blendedOutput = config["output_dir"]
blendedBase = config["output_base"]

Score = Moog960.Score(directory=config["convolved_dir"], suffix='', observed='../Theremin/TWHydra.fits')

blah = Score.listen()

"""
If your grid has multiple wavelength regions, you will first need to master the 
phrases into a single phrase which contains all the data.
"""

mastered = Score.master()
Score.selectMelodies(wlRange=[wlStart, wlStop])
#convolved = orchestra.getLabels(keySignature='CONVOLVED', selected=True)
Score.selectEnsemble(selectedLabels=mastered)

for T in desiredTeff:
    for G in desiredlogg:
        for B in desiredBfield:
            desiredParams = {}
            desiredParams["TEFF"] = T
            desiredParams["LOGG"] = G
            desiredParams["BFIELD"] = B
            blendedLabel = Score.blend(desiredParameters=desiredParams)

            filename = blendedOutput+blendedBase+'_T%d_G%.2f_B%.2f.fits' % (T, G, B)
            # Save fits files for further processing  - Very messy!  I know!!!
            #blendedLabels[0].Spectrum.bin(Score.compositeObservedLabel.Spectrum.wl)
            blendedLabel.Phrase.saveConvolved(label=blendedLabel, 
                      filename=filename, header=blendedLabel.Melody.header)

