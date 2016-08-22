import Moog960
import SpectralTools
import matplotlib.pyplot as pyplot
import scipy
import scipy.optimize
import numpy

def fit(synthetic, observed):
    veiling = 0.1
    slope = 0.0001
    continuum = 1.0

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
wlStart = 21700.0
#wlStop = 21980.0
#wlStart = 22010.0
#wlStop = 22550.0
#wlStart = 20000.0
#wlStop = 21000.0
#wlStart = 21000.0
#wlStop = 22000.0
#wlStart = 22000.0
wlStop = 22800.0

filename = '../Theremin/TWHydra.fits'
#TWHya = Moog960.Score(observed=filename)
#TWHya = Moog960.ObservedMelody.fromFile(filename=filename, label='IGRINS TWHydra')

#TWHya.selectPhrases(wlRange = [wlStart, wlStop])
#TWHya.loadData()

orchestra = Moog960.Score(directory='./TWHydra', suffix='', observed=filename)
#raw, interpolated, integrated, convolved, observed = orchestra.getMelodyParams()
observedSpectra, observedLabel = orchestra.listen()
observedSpectra.wl = observedSpectra.wl - 5.5
region = (observedSpectra.wl > wlStart) & (observedSpectra.wl < wlStop)
observedSpectra.wl = observedSpectra.wl[region]
observedSpectra.flux_I = observedSpectra.flux_I[region]


mastered=orchestra.master()
orchestra.selectMelodies(wlRange=[wlStart, wlStop])
convolved = orchestra.getLabels(keySignature='CONVOLVED', selected=True)
orchestra.selectEnsemble(selectedLabels=mastered)

CMJParams = {"TEFF":4180, "LOGG":4.8, "BFIELD":2.3}
WVParams = {"TEFF":3600, "LOGG":3.5, "BFIELD":0.0}
SpectrumCMJ, LabelCMJ = orchestra.blend(desiredParameters=CMJParams)
SpectrumWV, LabelWV = orchestra.blend(desiredParameters=WVParams)


#spectra, params = orchestra.perform(selectedLabels = [LabelCMJ, LabelWV])

#colors = ['g', 'm', 'c', 'r']
colors = ['g', 'r']

bestFit= []
veiled = []
observedSpectra.plot(ax=ax1, color='k')
#for T, G, B, color in zip(Ts, Gs, Bs, colors):
for label, color in zip([LabelCMJ[0], LabelWV[0]], colors):
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
raw_input()


bestFit = numpy.average(numpy.array(bestFit), axis=0)
#bestContinuum = 
#bestFit[2] + bestFit[1]*xpts
observed_spectra[0][0].flux_I = observed_spectra[0][0].flux_I/bestFit[2] - bestFit[1]*numpy.arange(len(observed_spectra[0][0].wl))
ax1.plot(observed_spectra[0][0].wl, observed_spectra[0][0].flux_I, color = 'k', lw=2.0)

ax1.set_xbound(wlStart, wlStop)
ax1.set_ybound(0.5, 1.1)

for v in veiled:
    sp = v["spectrum"]
    difference = sp - observed_spectra[0][0]
    ax2.plot(difference.wl, difference.flux_I, color = v["color"], lw = 2.0)
    print v["Teff"], v["B"], v["veiling"]
fig.show()

#del(orchestra)
#del(spectra)
#del(params)
#del(labels)
