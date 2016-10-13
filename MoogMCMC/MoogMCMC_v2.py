'''
Free Parameters:
B = B_field in kG, allow to vary between 0.0 and 4.0
T = T_eff in K, allow to vary between 3000.0 and 5000.0
logg = logg in dex, allow to vary between 3.0 and 5.0
v = velocity shift (RV + BVC) in km/s, allow to varry between -60.0 and 60.0

maybe veiling and continuum?

'''
import Moog960
import numpy as np
import emcee
import corner
import time
import matplotlib.pyplot as pyplot

fig = pyplot.figure(0)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

def lnprior(theta): #Priors for Bayesian Model, assuming flat priors
    B, T, logg, v = theta
    if 0.0 < B < 4.0 and 3000.0 < T < 5000.0 and 3.0 < logg < 5.0 and -60.0 < v <60.0:
        return 0.0
    return -np.inf

def lnlike(theta): #Log-likelyhood for the posterior, using chi-squared
    B, T, logg, v = theta
    retval = Score.calc_lnlike(Bfield=B, Teff=T, logg=logg, rv=v, ax=None)
    #synth = #casey's code to call the synth model
    #return -0.5(np.sum( ((obs-synth)/obs) ** 2.))
    return retval

def lnprob(theta): #The probability of accepting the new point. 
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta)


#Score = Moog960.Score(directory='/home/cdeen/Code/Python/cuddly-weasel/MusicMaker/Convolved/', suffix='')
Score = Moog960.Score(directory='/data/cdeen/Data/MoogPy/Convolved/', suffix='')
mastered = Score.master()

T_low = 3000.0
T_high = 5000.0
G_low = 3.0
G_high = 5.0
B_low = 0.0
B_high = 4.0
V_low = -60.0
V_high = 60.0

npoints = 50
temps = np.random.rand(npoints)*(T_high-T_low) + T_low
gravs = np.random.rand(npoints)*(G_high-G_low) + G_low
bs = np.random.rand(npoints)*(B_high-B_low) + B_low
vs = np.random.rand(npoints)*(V_high-V_low) + V_low

n = 0
for t, g, b, v in zip(temps, gravs, bs, vs):
    n = n + 1
    desiredParams = {}
    desiredParams["TEFF"] = t
    desiredParams["LOGG"] = g
    desiredParams["BFIELD"] = b
    
    compositeObservedLabel = Score.blend(desiredParameters=desiredParams, appendTheBlend=False)
    compositeObservedLabel.Spectrum.rv(v)
    compositeObservedLabel.Spectrum.nyquistSample(R=45000.0)
    Score.compositeObservedLabel = compositeObservedLabel

    #combine observed phrases into a master phrase
    #compositeObservedLabel = Score.listen()

    #MCMC parameters
    nwalkers = 50       #number of walkers used in fitting
    ndim = 4             #number of parameters being fit


    #Initialize the chain
    pos_min = np.array([0.0, 3000.0, 3.0, -60.0])
    pos_max = np.array([4.0, 5000.0, 5.0, 60.0])
    psize = pos_max - pos_min
    pos = [pos_min + psize*np.random.rand(ndim) for i in range(nwalkers)]

    #Setup ensamble sampler from emcee
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=())

    #Preform and time burn-in phase
    time0 = time.time()
    pos, prob, state  = sampler.run_mcmc(pos, 100)
    sampler.reset()
    time1=time.time()
    print 'Burn-in time was '+str(time1-time0)+' seconds.'

    #Perform MCMC fit
    time0 = time.time()
    pos, prob, state  = sampler.run_mcmc(pos, 500)
    time1=time.time()
    print 'Fitting time was '+str(time1-time0)+' seconds.'

    samples = sampler.flatchain
    outfile = open("Test_%d.out" % n, 'w')
    outfile.write("Input Parameters #%d: \n" % n)
    outfile.write("Teff = %dK\n" % t)
    outfile.write("log g = %.3fK\n" % g)
    outfile.write("Bfield = %.3fK\n" % b)
    outfile.write("vsin i = %.3fK\n" % v)
    for i in range(4):
        devs = corner.quantile(samples[:,i], [0.16, 0.5, 0.84])
        outfile.write("Pameter %d : %.3f %.3f %.3f\n" % (i, devs[0], devs[1], devs[2]))
    outfile.close()
    figure = corner.corner(samples)
    figure.savefig("test_%d.png" % n)
    del(sampler)
    del(desiredParams)
    del(samples)
