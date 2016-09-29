import Moog960
import astropy.io.fits as pyfits
import SpectralTools
import numpy
import sys
import AstroUtils
import glob


#This program reads in an IGRINS spectrum from the .fits files produced by the
#reduction pipeline, and saves it in a Moog960-readable .fits file.  It merges
#the spectral orders, weighting by signal-to-noise ratio.


config = AstroUtils.parse_config(sys.argv[1])

datadir = config["datadir"]
outdir = config['outdir']

wlsolfiles = [glob.glob(datadir+'SKY_SDCH*wvlsol_v1.fits')[0], glob.glob(datadir+'SKY_SDCK*wvlsol_v1.fits')[0]]
datafiles = [glob.glob(datadir+'SDCH*spec_a0v.fits')[0], glob.glob(datadir+'SDCK*spec_a0v.fits')[0]]
snfiles = [glob.glob(datadir+'SDCH*sn.fits')[0], glob.glob(datadir+'SDCK*sn.fits')[0]]
varfiles = [glob.glob(datadir+'SDCH*variance.fits')[0], glob.glob(datadir+'SDCK*variance.fits')[0]]

phrases = []
for wlsolfile, datafile, snfile, varfile in zip(wlsolfiles, datafiles, snfiles, varfiles):
    wlsol = pyfits.getdata(wlsolfile)
    data = pyfits.getdata(datafile)
    dataheader = pyfits.getheader(datafile)
    snr = pyfits.getdata(datafile)
    var = pyfits.getdata(varfile)
    observedData = []
    observedLabels = []
    headerKWs = {}
    Name = dataheader.get('OBJECT')
    headerKWs['OBJECT'] = Name
    headerKWs["DATE-OBS"] = dataheader.get('DATE-OBS')
    headerKWs["INSTRUMENT"] = dataheader.get('INSTRUME')
    headerKWs["OBSERVER"] = dataheader.get('OBSERVER')
    headerKWs["EXPTIME"] = dataheader.get('EXPTIME')
    Score = Moog960.Score()
    Melody = Moog960.ObservedMelody(Score=Score, name=Name)
    for wl, I, sn, v in zip(wlsol, data, snr, var):

        header = pyfits.Header()
        header.set('NAME', Name)
        header.set("WLSTART", numpy.min(wl))
        header.set("WLSTOP", numpy.max(wl))

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

    Score.record(keySignature='OBSERVED', filename = outdir+Name+'.fits', dI=True,
                 primaryHeaderKWs=headerKWs)
    label = Score.listen()
    label.Phrase.record(filename = outdir+Name+'_composite.fits', dI=True,
                        primaryHeaderKWs=headerKWs)

                                                                                   
del(phrases)
del(Score)
del(Melody)
del(label)                                         
                                                                                   
                                                                                     