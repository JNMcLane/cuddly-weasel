'''
Free Parameters:
B = B_field in kG, allow to vary between 0.0 and 4.0
T = T_eff in K, allow to vary between 3000.0 and 5000.0
logg = logg in dex, allow to vary between 3.0 and 5.0
v = velocity shift (RV + BVC) in km/s, allow to varry between -60.0 and 60.0

maybe vieling and continuum?

'''
import Moog960
import numpy as np
import emcee
import corner
import time
import matplotlib.pyplot as pyplot

wlStart = 22000.0
wlStop = 22100.0

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


Score = Moog960.Score(directory='../MusicMaker/Convolved/', suffix='', observed='../MusicMaker/Blended/blended_T3550_G3.75_B0.50.fits')
mastered = Score.master()

#combine observed phrases into a master phrase
compositeObservedLabel = Score.listen(wlRange=[wlStart, wlStop])

#MCMC parameters
nwalkers = 50 	#number of walkers used in fitting
ndim = 4 		#number of parameters being fit


#Initialize the chain
pos_min = np.array([0.0, 3000.0, 3.0, -60.0])
pos_max = np.array([4.0, 5000.0, 5.0, 60.0])
psize = pos_max - pos_min
pos = [pos_min + psize*np.random.rand(ndim) for i in range(nwalkers)]

#Setup ensamble sampler from emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(), threads=6)

#Preform and time burn-in phase
time0 = time.time()
pos, prob, state  = sampler.run_mcmc(pos, 300)
time1=time.time()
print 'Burn-in time was '+str(time1-time0)+' seconds.'

chi = sampler.lnprobability
for i in chi:
    ax.plot(i)

fig.show()
raw_input()

sampler.reset()

#Perform MCMC fit
time0 = time.time()
pos, prob, state  = sampler.run_mcmc(pos, 100)
time1=time.time()
print 'Fitting time was '+str(time1-time0)+' seconds.'


samples = sampler.flatchain
figure = corner.corner(samples)
figure.show()
figure.savefig("test.png")
