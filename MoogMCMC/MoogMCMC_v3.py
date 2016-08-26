'''
Free Parameters:
B = B_field in kG, allow to vary between 0.0 and 4.0
T = T_eff in K, allow to vary between 3000.0 and 5000.0
logg = logg in dex, allow to vary between 3.0 and 5.0
v = velocity shift (RV + BVC) in km/s, allow to varry between -60.0 and 60.0

maybe vieling and continuum?

'''
import Moog961
import numpy as np
import emcee
import corner
import time
import matplotlib.pyplot as plt

b_true=1.62
t_true=4570
g_true=3.09
v_true=0.0

fig1 = plt.figure(0)
fig1.clear()
ax1 = fig1.add_axes([0.1, 0.1, 0.8, 0.8])

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


Score = Moog961.Score(directory='../MusicMaker/vsini_10/Model_T4*', suffix='', observed='../Tutorial/Blended/Model_T4570_G3.09_B1.62.fits')
mastered = Score.master()

#combine observed phrases into a master phrase
compositeObservedSpectrum, compositeObservedLabel = Score.listen()

blended, blendedLabel = Score.blend(desiredParameters={"TEFF":4570, "LOGG":3.09, "BFIELD":1.62})

compositeObservedLabel.Spectrum.plot(ax=ax1)
blendedLabel[0].Spectrum.plot(ax=ax1)

fig1.show()
raw_input()

"""

#MCMC parameters
nwalkers = 8 	#number of walkers used in fitting, must be even and >= 2x ndim
ndim = 4 		#number of parameters being fit


#Initialize the chain
pos_min = np.array([1.62, 4570.0, 3.09, 0.0])
pos_max = np.array([1.62, 4570.0, 3.09, 0.0])
psize = pos_max - pos_min
pos = [pos_min + psize*np.random.rand(ndim) for i in range(nwalkers)]

#Setup ensamble sampler from emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=())

#Preform and time burn-in phase
print 'Running first burn-in.'
time0 = time.time()
pos, prob, state  = sampler.run_mcmc(pos, 1)
sampler.reset()
time1=time.time()
print 'First burn-in finished, '+str(time1-time0)+' seconds elapsed.'
print 'Running second burn-in.'
p=pos[np.argmax(prob)]
pos=[p+(psize/100)*np.random.randn(ndim) for i in xrange(nwalkers)]
pos, prob, state  = sampler.run_mcmc(pos, 1)
sampler.reset()
time2=time.time()
print 'Second burn-in finished, '+str(time2-time1)+' seconds elapsed.'
print 'Burn-in time was '+str(time2-time0)+' seconds.'

#Perform MCMC fit
time0 = time.time()
pos, prob, state  = sampler.run_mcmc(pos, 1)
time1=time.time()
print 'Fitting time was '+str(time1-time0)+' seconds.'

samples = sampler.flatchain

fig =corner.corner(samples, labels=["$B$","$T_{eff}$","$logg$","$RV$"], truths=[b_true,t_true,g_true,v_true])
fig.savefig('test_2.png')

"""
