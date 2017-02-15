import numpy
import scipy
import matplotlib.pyplot as pyplot
import AstroUtils
import MoogTools
import Moog960
import SpectralTools
import astropy.io.fits as pyfits
import random
import string
import sys

def parseParams(labels):
    keys = labels[0].parameters.keys()
    params = {}
    for key in keys:
        params[key] = []
        
    for label in labels:
        for key in keys:
            params[key].append(label.parameters[key])
            
    for key in keys:
        params[key] = numpy.unique(params[key])
        
    return params


class ControlMatrix( object ):
    def __init__(self, IMs = [], factors=[], coordinates = {}, Gain=0.5, observed = None, ax=None):
        self.IM_Spectra = IMs
        self.observed = observed
        self.IMs = []
        self.factors = factors
        self.ax = ax
        self.nFiltModes = 4
        self.calcCM()
        self.Gain = Gain
        self.coordinates = coordinates

    def calcCM(self):
        for spec, factor in zip(self.IM_Spectra, self.factors):
            #self.IMs.append(spec.flux_I/factor)
            self.IMs.append(spec.flux_I)
        self.IMs = numpy.array(self.IMs)
        #"""
        #lineIMs = self.IM[:-4:,]
        dims = self.IMs.shape
        U,S,V = scipy.linalg.svd(self.IMs)
        if self.ax != None:
            self.ax.clear()
            self.ax.plot(numpy.log10(S))
            self.ax.figure.show()
            blah = raw_input("Enter number of modes to keep: ")
            try:
                nFiltModes = dims[0] - int(blah)
                self.nFiltModes = nFiltModes
            except:
                pass
            self.ax.clear()
        if self.nFiltModes == 0:
            D = 1.0/S
        else:
            D = 1.0/(S[0:-self.nFiltModes])
            S[-self.nFiltModes:] = 0.0
        newS = numpy.zeros((dims[0], dims[1]))
        I = [i for i in range(dims[1])]
        for i in range(len(D)):
            newS[i][i] = D[i]

        S = newS.copy()
        self.CM = numpy.array(scipy.matrix(V.T.dot(S.T.dot(U.T)),dtype=numpy.float32)).T
        #"""

        #self.CM = scipy.linalg.pinv(self.IMs).T
        #blah = scipy.linalg.pinv(self.IMs).T
        #if self.ax != None:
        #    self.ax.clear()
        #    self.ax.matshow(self.CM)
        #    self.ax.figure.show()
        #    self.ax.set_aspect('auto')

        #print asdf
        """
        globalIMs = self.IM[-4:,]
        dims = globalIMs.shape
        U,S,V = scipy.linalg.svd(globalIMs)
        D = 1.0/(S[0:-1])
        S[-1:] = 0.0
        newS = numpy.zeros((dims[0], dims[1]))
        I = [i for i in range(dims[1])]
        for i in range(len(D)):
            newS[i][i] = D[i]

        S = newS.copy()
        globalCM = numpy.array(scipy.matrix(V.T.dot(S.T.dot(U.T)),dtype=numpy.float32)).T

        self.CM = numpy.zeros((lineCM.shape[0]+globalCM.shape[0], globalCM.shape[1]))
        for i in range(lineCM.shape[0]):
            self.CM[i] = lineCM[i]
        for j in range(globalCM.shape[0]):
            self.CM[i+j+1] = globalCM[j]

        """
        

    def dot(self, observed, difference):
        overlap_start = numpy.max([numpy.min(observed.wl), numpy.min(difference.wl)])
        overlap_stop = numpy.min([numpy.max(observed.wl), numpy.max(difference.wl)])
        overlap = scipy.where((observed.wl >= overlap_start) & (observed.wl <= overlap_stop))

        diff = numpy.zeros(len(observed.wl))

        diff[overlap] = difference.flux_I

        return self.CM.dot(diff)

    def getCommand(self, observed, synthesized, ax=None):
        difference = observed - synthesized

        if ax!=None:
            ax.clear()
            difference.plot(ax=ax)
            ax.plot(self.observed.wl, self.CM[-2])
            ax.plot(self.observed.wl, self.IM[-2])
            ax.figure.show()
            raw_input()
        command =  self.Gain*(self.dot(observed, difference))*self.factors
        #command[self.factors < 1e-3] = 0.0
        return command

def computeCMs(orchestra, ax):
    convolved = orchestra.convolved_labels

    mastered = orchestra.master(selectedLabels=convolved, keySignature='CONVOLVED')
    test = orchestra.blend(desiredParameters={"TEFF":3000, "LOGG":3.00, "BFIELD":0.7},
                           appendTheBlend=False)
    test.Spectrum.nyquistSample(R=45000.0)
    test.Spectrum.rv(5.32)
    zero = orchestra.blend(desiredParameters={"TEFF":3000, "LOGG":3.00, "BFIELD":0.0},
                           appendTheBlend=False)
    zero.Spectrum.rv(5.32)
    one = orchestra.blend(desiredParameters={"TEFF":3000, "LOGG":3.00, "BFIELD":0.1},
                           appendTheBlend=False)
    one.Spectrum.rv(5.32)
    two = orchestra.blend(desiredParameters={"TEFF":3000, "LOGG":3.00, "BFIELD":0.2},
                           appendTheBlend=False)
    two.Spectrum.rv(5.32)
    three = orchestra.blend(desiredParameters={"TEFF":3000, "LOGG":3.00, "BFIELD":0.3},
                           appendTheBlend=False)
    three.Spectrum.rv(5.32)
    orchestra.compositeObservedLabel=test
    observed = orchestra.compositeObservedLabel

    params = parseParams(mastered)

    wlStart = numpy.min(params["WLSTART"])
    wlStop = numpy.max(params["WLSTOP"])

    for blah in mastered:
        blah.parameters["SELECTED"] = True

    TeffStroke = 50.0
    loggStroke = 0.5
    BfieldStroke = 1.0
    Tmin = numpy.min(params["TEFF"])
    Tmax = numpy.max(params["TEFF"])
    Gmin = numpy.min(params["LOGG"])
    Gmax = numpy.max(params["LOGG"])
    Bmin = numpy.min(params["BFIELD"])
    Bmax = numpy.max(params["BFIELD"])
    
    CMs = []

    for T in params["TEFF"]:
        for G in params["LOGG"]:
            for B in params["BFIELD"]:
                print("T=%d, G=%.2f B=%.2f" % (T, G, B))
                IMs = []
                factors = []
                orig = orchestra.blend(desiredParameters = {"TEFF":T, "LOGG":G,
                         "BFIELD":B}, appendTheBlend=False)
                orig.Spectrum.bin(newWl=observed.Spectrum.wl)

                #Change in Teff
                plusT = numpy.min([T+TeffStroke, Tmax])
                plus = orchestra.blend(desiredParameters = {"TEFF":plusT, 
                    "LOGG":G, "BFIELD":B}, appendTheBlend=False)

                minusT = numpy.max([T-TeffStroke, Tmin])
                minus = orchestra.blend(desiredParameters = {"TEFF":minusT, 
                    "LOGG":G, "BFIELD":B}, appendTheBlend=False)

                IMs.append(plus.Spectrum - minus.Spectrum)
                IMs[-1].bin(newWl=observed.Spectrum.wl)
                factors.append(plus.parameters["TEFF"] - minus.parameters["TEFF"])

           
                #Change in Surface Gravity
                plusG = numpy.min([G+loggStroke, Gmax])
                plus = orchestra.blend(desiredParameters = {"TEFF":T, 
                    "LOGG":plusG, "BFIELD":B}, appendTheBlend=False)

                minusG = numpy.max([G-loggStroke, Gmin])
                minus = orchestra.blend(desiredParameters = {"TEFF":T, 
                    "LOGG":minusG, "BFIELD":B})

                IMs.append(plus.Spectrum-minus.Spectrum)
                IMs[-1].bin(newWl=observed.Spectrum.wl)
                factors.append(plus.parameters["LOGG"] - minus.parameters["LOGG"])

                #Change in B-field
                plusB = numpy.min([B+BfieldStroke, Bmax])
                plus = orchestra.blend(desiredParameters = {"TEFF":T, 
                    "LOGG":G, "BFIELD":plusB}, appendTheBlend=False)

                minusB = numpy.max([B-BfieldStroke, Bmin])
                minus = orchestra.blend(desiredParameters = {"TEFF":T, 
                    "LOGG":G, "BFIELD":minusB}, appendTheBlend=False)

                IMs.append(plus.Spectrum-minus.Spectrum)
                IMs[-1].bin(newWl=observed.Spectrum.wl)
                factors.append(plus.parameters["BFIELD"] - minus.parameters["BFIELD"])

                #Radial Velocity shift
                red = orig.Spectrum * 1.0
                red.rv(10.0)
                blue = orig.Spectrum * 1.0
                blue.rv(-10.0)
                rv = red - blue
                IMs.append(rv)
                IMs[-1].bin(newWl=observed.Spectrum.wl)
                factors.append(20.0)
                del(red)
                del(blue)

                #Continuum Shift
                continuum_shift = (orig.Spectrum * 1.1) - (orig.Spectrum * 0.9)
                IMs.append(continuum_shift)
                factors.append(0.2)

                #Veiling

                CMs.append(ControlMatrix(IMs=IMs, factors=factors,
                           coordinates = {"TEFF":T, "LOGG":G, "BFIELD":B}))

                params = {"TEFF":3000, "LOGG": 3.0, "BFIELD":0.0}
                rv = 5.32
                while True:
                    ax1.clear()
                    orchestra.compositeObservedLabel.Spectrum.plot(ax=ax1, color = 'k')
                    zero.Spectrum.plot(ax=ax1, color = 'b')
                    guess = orchestra.blend(desiredParameters=params, appendTheBlend=False)
                    guess.Spectrum.rv(rv)
                    guess.Spectrum.bin(newWl = orchestra.compositeObservedLabel.Spectrum.wl)
                    guess.Spectrum.plot(ax=ax1, color = 'm')

                    command = CMs[-1].getCommand(orchestra.compositeObservedLabel.Spectrum, guess.Spectrum)
                    print params
                    print command
                    params["TEFF"] += command[0]
                    params["TEFF"] = numpy.max([3000.0, params["TEFF"]])
                    params["LOGG"] += command[1]
                    params["LOGG"] = numpy.max([3.0, params["LOGG"]])
                    params["BFIELD"] += command[2]
                    params["BFIELD"] = numpy.max([0.0, params["BFIELD"]])
                    rv += command[3]
                    f1.show()
                    #print params
                    raw_input()
                del(plus)
                del(minus)
                del(orig)

    return CMs


def findClosestCM(CMs, params):
    index = 0
    for i in range(len(CMs)):
        for key in params.keys():
            if (numpy.abs(CMs[i].coordinates[key] - params[key]) < numpy.abs(CMs[index].coordinates[key] - params[key])):
                index = i
    return index

def fitParams(orchestra, CMs, ax):
    params = {"TEFF":3600, "LOGG":4.2, "BFIELD":2.1}
    rv = 3.5
    continuum = 1.0

    while True:
        guess = orchestra.blend(desiredParameters=params, appendTheBlend=False)
        guess.Spectrum.rv(rv)
        guess.Spectrum.bin(newWl = orchestra.compositeObservedLabel.Spectrum.wl)
        
        i = findClosestCM(CMs, params)
        print CMs[i].coordinates
        print params
        print rv
        command = CMs[i].getCommand(orchestra.compositeObservedLabel.Spectrum, guess.Spectrum)
        print command
        params["TEFF"] += command[0]
        params["LOGG"] += command[1]
        params["BFIELD"] += command[2]
        rv += command[3]
        ax.clear()
        guess.Spectrum.plot(ax=ax)
        orchestra.compositeObservedLabel.Spectrum.plot(ax=ax)
        ax.figure.show()
        raw_input()
        del(guess)

configFile = sys.argv[1]
calcIM = bool(sys.argv[2]=='True')
moogInstance = sys.argv[3]
try:
    contFloat = bool(sys.argv[4]== 'True')
except:
    contFloat = False
#try:
#    rotFloat = bool(sys.argv[5] == 'True')
#except:
#    rotFloat = False

f1 = pyplot.figure(0)
f1.clear()
ax1 = f1.add_axes([0.1, 0.1, 0.8, 0.8])

config = AstroUtils.parse_config(configFile)

orchestra = Moog960.Score(directory='../MusicMaker/Convolved/Kband_Convolved_T3*', observed=config["inFile"], suffix='')

#mastered = orchestra.master(selectedLabels=orchestra.convolved_labels, keySignature='CONVOLVED')

#one = mastered[0].Spectrum.copy()

#one.plot(ax=ax1)
#one.nyquistSample(R=45000.0)
#one.plot(ax=ax1)
#mastered[0].Spectrum.plot(ax=ax1)
#f1.show()

#print asdf
if calcIM:
    CMs = computeCMs(orchestra, ax1)

    #HDU = pyfits.PrimaryHDU(CMs)
    #HDU.writeto("ControlMatrices.fits")
    #pickle.dump(CMs, "ControlMatrices.dat")
else:
    CMs = pickle.load("ControlMatrices.dat")

fit = fitParams(orchestra, CMs, ax1)

print asdf
twhya = observed[0][0]
twhya.wl += wlOffset
twhya.flux_I += continuum

window = ((twhya.wl > wlRange[0]) & (twhya.wl < wlRange[1]))

twhya.wl = twhya.wl[window]
twhya.flux_I = twhya.flux_I[window]

Synth.run()
nominal = Synth.Spectra[0]
smoothed = nominal.resample(R=resolution)
resampled = nominal.resample(R=resolution, observedWl=twhya.wl, pad=1.0)
Spectra = [resampled]

if calcIM:
    calculateIM(Synth, resampled, resolution)

Gain = 0.15
IM = pyfits.getdata(Synth.config["outputDir"]+Synth.config["baseName"]+"_IM.fits")
factor = pyfits.getdata(Synth.config["outputDir"]+Synth.config["baseName"]+"_factors.fits")
gfIndices = pyfits.getdata(Synth.config["outputDir"]+Synth.config["baseName"]+'_gf.fits')
vdWIndices = pyfits.getdata(Synth.config["outputDir"]+Synth.config["baseName"]+'_vdW.fits')
CM = ControlMatrix(IM, factor, solar, Gain, ax=ax1)
#CM = ControlMatrix(IM, factor, solar, nFiltModes, Gain)

hdu = pyfits.PrimaryHDU(CM.CM)
hdu.writeto(Synth.config["outputDir"]+Synth.config["baseName"]+"_CM.fits", clobber=True)

cb = ax1.matshow(CM.CM, aspect='auto')
blah = pyplot.colorbar(cb)
f1.show()

f2 = pyplot.figure(1)
f2.clear()
ax2 = f2.add_axes([0.1, 0.1, 0.8, 0.8])

f3 = pyplot.figure(2)
f3.clear()
ax3 = f3.add_axes([0.1, 0.1, 0.8, 0.8])

nLines = Synth.lineList.numLines

continuum = 0.00
wlOffset = 0.0
rotation = 0.0

for j in range(50):
    #command = CM.getCommand(Spectra[-1], ax=ax3)
    command = CM.getCommand(Spectra[-1])
    
    for i in range(len(gfIndices)):
        if gfIndices[i] != -1:
            Synth.lineList.perturbGf(i, command[gfIndices[i]], push=True)
    for i in range(len(vdWIndices)):
        if vdWIndices[i] != -1:
            Synth.lineList.perturbVdW(i, command[vdWIndices[i]], push=True)

    if contFloat:
        continuum = continuum+command[-4]
    else:
        continuum = 0.0
    wlOffset = wlOffset+command[-3]
    #resolution = resolution*(1.0+command[-2])
    if rotFloat:
        rotation = rotation+command[-1]
    else:
        rotation = 0.0
    Synth.run()
    output = Synth.Spectra[-1].resample(R=resolution).rotate(angle=rotation)
    Synth.Spectra = []
    output.wl += wlOffset
    output.flux_I += continuum
    output.bin(solar.wl)
    Spectra.append(output)

    """
    ax2.clear()
    if diffplot:
        difference = solar.diff_spectra(Spectra[-1], pad=True)
        difference.flux_I[Spectra[-1].flux_I == 0] = 0.0
        ax2.plot(difference.wl, difference.flux_I)
    else:
        ax2.plot(solar.wl, solar.flux_I, lw=2.0)
        for spec in Spectra:
            ax2.plot(spec.wl, spec.flux_I)
    f2.show()
    #"""
    """
    ax3.clear()
    for line in range(Synth.lineList.nStrong):
        ax3.plot(Synth.lineList.strongLines[line].loggfHistory)
        ax3.plot(Synth.lineList.strongLines[line].VdWHistory)
    for line in range(Synth.lineList.numLines - Synth.lineList.nStrong):
        ax3.plot(Synth.lineList.weakLines[line].loggfHistory)
        ax3.plot(Synth.lineList.weakLines[line].VdWHistory)
    f3.show()
    #"""

    print j, continuum, wlOffset, resolution, rotation
    power = numpy.mean((command**2.0)**(0.5))
    print Synth.config["baseName"], numpy.std(command), numpy.mean((command**2.0)**(0.5))
    #raw_input()

ax2.clear()
ax2.plot(solar.wl, solar.flux_I, lw=2.0)
ax2.plot(Spectra[0].wl, Spectra[0].flux_I)
ax2.plot(Spectra[-1].wl, Spectra[-1].flux_I)
f2.show()
f2.savefig(Synth.config["outputDir"] + Synth.config["baseName"]+'_results.png')
Synth.lineList.saveLineList(filename=outFile, changed=True)
