import Moog960
import SpectralTools
import matplotlib.pyplot as pyplot
import scipy
import scipy.optimize
import numpy

def calc_VeilingSED(Thot, Tint, Tcool, fi_fh, fc_fh):
    wave = numpy.arange(1.0, 2.5, 0.01)
    bm = numpy.argsort(abs(wave-2.2))[0]
    Bhot = SpectralTools.blackBody(wl = wave/10000.0, T=Thot,
               outUnits='Energy')
    Bint = SpectralTools.blackBody(wl = wave/10000.0, T=Tint,
               outUnits='Energy')*fi_fh
    Bcool = SpectralTools.blackBody(wl = wave/10000.0, T=Tcool,
               outUnits='Energy')*fc_fh
    composite = Bhot+Bint+Bcool
    composite /= composite[bm]
    return scipy.interpolate.interp1d(wave*10000.0, composite)


def fit(synthetic, observed):
    veiling = 0.1
    slope = 0.0001
    continuum = 1.0
    veilingSED = calc_VeilingSED(10000, 5000, 1500, 0.1, 0.7)

    params = [veiling, slope, continuum]

    def fitfunc(pars, synth):
        xpts = numpy.arange(len(synth))

        retval = (synth + pars[0])/(1.0+pars[0])*pars[2] + pars[1]*xpts

        return retval
        

    def errfunc(pars, synth, obs):
        return numpy.abs(fitfunc(pars, synth) - obs)

    bestFit, success = scipy.optimize.leastsq(errfunc, params, args=(synthetic, observed))

    #return (synthetic+bestFit[0])/(1.0+bestFit[0])*bestFit[2] + bestFit[1]*xpts
    return bestFit

fig = pyplot.figure(0)
fig.clear()
ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.4])
ax2 = fig.add_axes([0.1, 0.5, 0.8, 0.4])

#wlStart = 22880.0
#wlStop = 24000.0
#wlStart = 22700.0
#wlStop = 22890.0
#wlStart = 22100.0
#wlStop = 22550.0
#wlStart = 21700.0
#wlStop = 21980.0
#wlStart = 22010.0
#wlStop = 22550.0
#wlStart = 20000.0
#wlStop = 21000.0
#wlStart = 21000.0
#wlStop = 22000.0
wlStart = 22800.0
wlStop = 23800.0

filename = '../Transpose/TWHydra_20150128_composite.fits'

orchestra = Moog960.Score(directory='./TWHydra', suffix='', observed=filename)
observedSpectra, observedLabel = orchestra.listen()
observedSpectra.wl *= 10.0
observedSpectra.wl = observedSpectra.wl - 5.5
median = numpy.median(observedSpectra.flux_I)
observedSpectra.flux_I /= median
observedSpectra.dflux_I /= median
region = (observedSpectra.wl > wlStart) & (observedSpectra.wl < wlStop)
observedSpectra.wl = observedSpectra.wl[region]
observedSpectra.flux_I = observedSpectra.flux_I[region]
observedSpectra.dflux_I = observedSpectra.dflux_I[region]


mastered=orchestra.master()
orchestra.selectMelodies(wlRange=[wlStart, wlStop])
convolved = orchestra.getLabels(keySignature='CONVOLVED', selected=True)
orchestra.selectEnsemble(selectedLabels=mastered)

DSMParams = {"TEFF":4180, "LOGG":4.0, "BFIELD":2.3}
CMJParams = {"TEFF":4180, "LOGG":4.8, "BFIELD":2.3}
WVParams = {"TEFF":3600, "LOGG":3.5, "BFIELD":0.0}
SpectrumCMJ, LabelCMJ = orchestra.blend(desiredParameters=CMJParams)
SpectrumWV, LabelWV = orchestra.blend(desiredParameters=WVParams)
SpectrumDSM, LabelDSM = orchestra.blend(desiredParameters=DSMParams)


# Save fits files for further processing  - Very messy!  I know!!!
LabelCMJ[0].Phrase.saveConvolved(label=LabelCMJ[0], filename='./Output/CMJ.fits', header=LabelCMJ[0].Melody.header)
LabelDSM[0].Phrase.saveConvolved(label=LabelDSM[0], filename='./Output/DSM.fits', header=LabelDSM[0].Melody.header)
LabelWV[0].Phrase.saveConvolved(label=LabelWV[0], filename='./Output/WV.fits', header=LabelWV[0].Melody.header)

#spectra, params = orchestra.perform(selectedLabels = [LabelCMJ, LabelWV])

#colors = ['g', 'm', 'c', 'r']
colors = ['g', 'b', 'r']

bestFit= []
veiled = []
observedSpectra.plot(ax=ax1, color='k')
for label, color in zip([LabelCMJ[0], LabelWV[0], LabelDSM[0]], colors):
    spectrum = label.Spectrum
    spectrum.bin(observedSpectra.wl)
    bestFit.append(fit(spectrum.flux_I, observedSpectra.flux_I))
    print bestFit
    xpts = numpy.arange(len(spectrum.flux_I))
    spectrum.flux_I = (spectrum.flux_I+bestFit[-1][0])/(1.0+bestFit[-1][0])*bestFit[-1][2] + xpts*bestFit[-1][1]
    #ax1.clear()
    spectrum.plot(ax=ax1, color = color, lw=2.0)
    #observedSpectra.plot(ax=ax1, color = 'k')
    diffSpectra = observedSpectra - spectrum
    diffSpectra.plot(ax=ax2, color = color, lw=0.5)
    v = {}
    v["spectrum"] = spectrum
    v["Teff"] = label.parameters["TEFF"]
    v["log g"] = label.parameters["LOGG"]
    v["B"] = label.parameters["BFIELD"]
    v["color"] = color
    v["veiling"] = bestFit[-1][0]
    veiled.append(v)
    #fig.show()
    #raw_input()

fig.show()
