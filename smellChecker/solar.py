import MoogTools
import astropy.io.fits as pyfits
import Moog960
import matplotlib.pyplot as pyplot

fig = pyplot.figure(0)
fig.clear()
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

#solar = Moog960.ObservedMelody.fromFile(filename='SolarSpectrum.fits')
#solar.loadData()
solar = Moog960.Score(observed='SolarSpectrum.fits')
arcturus = Moog960.Score(observed='ArcturusSpectrum.fits')

wlStart = 17000
wlStop = 17100

observed = solar.getLabels(keySignature='OBSERVED')

for obs in observed:
    obs.Spectrum.plot(ax=ax, color = 'g')

observed = arcturus.getLabels(keySignature='OBSERVED')

for obs in observed:
    obs.Spectrum.plot(ax=ax, color = 'b')



ArcSynth = MoogTools.MoogStokes('Arcturus.cfg', wlStart=wlStart, wlStop=wlStop)
SolarSynth = MoogTools.MoogStokes('Solar.cfg', wlStart=wlStart, wlStop=wlStop)

SolarSynth.run()
ArcSynth.run()

SolarSynth.Spectra[0].plot(ax=ax, color = 'r')
SolarSynth.Spectra[1].plot(ax=ax, color = 'm')

ax.set_xbound(wlStart, wlStop)
ax.set_ybound(0.0, 1.0)
fig.show()
