import Moog960
import matplotlib.pyplot as pyplot
import AstroUtils
import sys
import numpy

"""

This program should be run fourth in the tutorial sequence.  This program allows the
user to read in a previously stored model to play around with it.

"""

#Loads in configuration
config = AstroUtils.parse_config(sys.argv[1])
filenames = numpy.array([x.strip() for x in config["desiredFiles"].split(',')])

#Loads in desired file
Score = Moog960.Score(files = filenames)
labels = Score.getLabels(keySignature='CONVOLVED')

#Sets up the plot

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

for label in labels:
    label.Spectrum.plot(ax=ax)

ax.legend(loc=3)
fig.show()

