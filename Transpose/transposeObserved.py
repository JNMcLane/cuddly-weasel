import Moog960
import SpectralTools
import numpy
import astropy.io.fits as pyfits
import matplotlib.pyplot as pyplot

fig1 = pyplot.figure(1)
fig1.clear()
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])


wlsolfiles = ['/home/cdeen/Data/20150128/SDCH_20150128_0131.wave.fits', '/home/cdeen/Data/20150128/SDCK_20150128_0131.wave.fits']
datafiles = ['/home/cdeen/Data/20150128/SDCH_20150128_0125.spec_a0v.fits', '/home/cdeen/Data/20150128/SDCK_20150128_0125.spec_a0v.fits']
snfiles = ['/home/cdeen/Data/20150128/SDCH_20150128_0125.sn.fits', '/home/cdeen/Data/20150128/SDCK_20150128_0125.sn.fits']
varfiles = ['/home/cdeen/Data/20150128/SDCH_20150128_0125.variance.fits', '/home/cdeen/Data/20150128/SDCK_20150128_0125.variance.fits']

#wlsolfiles = ['/home/cdeen/Data/EL24/SKY_SDCH_20140713_0021.wvlsol_v1.fits', '/home/cdeen/Data/EL24/SKY_SDCK_20140713_0021.wvlsol_v1.fits']
#datafiles = ['/home/cdeen/Data/EL24/SDCH_20140713_0010.spec_a0v.fits', '/home/cdeen/Data/EL24/SDCK_20140713_0010.spec_a0v.fits']
#snfiles = ['/home/cdeen/Data/EL24/SDCH_20140713_0010.sn.fits', '/home/cdeen/Data/EL24/SDCK_20140713_0010.sn.fits']
#aovfiles = ['/home/cdeen/Data/EL24/SDCH_20140713_0010.spec.fits', '/home/cdeen/Data/EL24/SDCK_20140713_0010.spec.fits']

phrases = []
for datafile, wlsolfile, snfile, varfile in zip(datafiles, wlsolfiles, snfiles, varfiles):
    Score = Moog960.Score()
    Melody = Moog960.ObservedMelody(Score=Score)
    wlsol = pyfits.getdata(wlsolfile)
    data = pyfits.getdata(datafile)
    snr = pyfits.getdata(snfile)
    var = pyfits.getdata(varfile)
    observedData = []
    observedLabels = []
    for wl, I, sn, v in zip(wlsol, data, snr, var):

        header = pyfits.Header()
        header.set("WLSTART", numpy.min(wl))
        header.set("WLSTOP", numpy.max(wl))

        #noise = SpectralTools.Spectrum(wl=wl, I=sn, header=header, spectrum_type="SNR")
        

        Melody.addPhrase(Moog960.ObservedPhrase(observedData=[SpectralTools.Spectrum(wl=wl, 
                         I=I, dI=I/sn, header=header, spectrum_type="OBSERVED")], Melody=Melody))

        ID = Melody.phrases[-1].ID
        parameters = {}
        parameters["PHRASE"] = ID
        parameters["MELODY"] = Melody.ID
        parameters["SCORE"] = Melody.Score.ID
        parameters["WLSTART"] = header.get('WLSTART')
        parameters["WLSTOP"] = header.get('WLSTOP')
        parameters["SELECTED"] = False
        if not(Melody.ID in Melody.Score.observed_labels.keys()):
            Melody.Score.observed_labels[Melody.ID] = {}
        Melody.Score.observed_labels[Melody.ID][ID] = []

        Melody.Score.observed_labels[Melody.ID][ID].append(Moog960.Label(parameters, 
                  reference=Melody.phrases[-1].observedData[-1], Phrase=Melody.phrases[-1],
                  Spectrum=Melody.phrases[-1].observedData[-1], Score=Score, Melody=Melody))
        
    Score.record(keySignature='OBSERVED', filename = 'TWHydra_20150128.fits', dI=True)
    composite, label = Score.listen()
    label.Phrase.record(filename = 'TWHydra_20150128_composite.fits', dI=True)
    composite.plot(ax=ax1)
    fig1.show()
    raw_input()



