import Moog960
import numpy
import SpectralTools
import matplotlib.pyplot as pyplot

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
ax.clear()


syntheticScore = Moog960.Score(suffix='', observed = '../../Theremin/TWHydra.fits')


models = syntheticScore.getLabels(keySignature='CONVOLVED')

observed = syntheticScore.getLabels(keySignature='OBSERVED')

for model in models:
    parameters = model.parameters
    label = "T=%d log g=%.1f B=%.1f" % (parameters["TEFF"], parameters["LOGG"], parameters["BFIELD"])
    model.Spectrum.plot(ax=ax, label=label)

for obs in observed:
    # Account for wavelength shift
    obs.Spectrum.wl -= 5.5
    obs.Spectrum.plot(ax=ax)

ax.legend()
fig.show()
