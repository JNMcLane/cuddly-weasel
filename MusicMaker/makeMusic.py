import scipy
import numpy
import matplotlib.pyplot as pyplot
import Moog960
import SpectralTools
import glob

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

vsini = 10.0
R = 1000.0

datadir = '/Users/JNMcLane/data/CorrectedRawData/'


files = glob.glob(datadir+'*_raw.fits')

for filename in files:
    print filename

    Score = Moog960.Score()
    syntheticMelody = Moog960.SyntheticMelody(filename=filename, Score=Score)
    syntheticMelody.selectPhrases(selectAll=True)
    syntheticMelody.rehearse(vsini=vsini, R=R)
    conv_Phrases = Score.convolved_labels[syntheticMelody.ID].keys()
    for phrase in conv_Phrases:
        conv_Labels = Score.convolved_labels[syntheticMelody.ID][phrase]
        syntheticMelody.record(labels=conv_Labels, basename='TWHydra')
    del(Score)
    del(syntheticMelody)
    del(conv_Phrases)
    del(conv_Labels)

