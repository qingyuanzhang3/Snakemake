import matplotlib.pyplot as plt
import numpy as np
import ticktack
from ticktack import fitting
import os
from matplotlib.lines import Line2D
import matplotlib as mpl
mpl.style.use('seaborn-colorblind')

lat_prolonged = [50 + 18/60, 54 + 28/60, 40 + 6/60,
       67 + 29/60, 31 + 7/60,
       72 + 13/60, 69,
       48 + 10/60, 67 + 31/60,
       35 + 10/60, 68 + 15/60, 68.5]
events_prolonged = ["Büntgen18_ALT01",
          "Büntgen18_CAN06", "Büntgen18_GRE02",
          "Büntgen18_RUS04", "Büntgen18_TIB01",
          "Büntgen18_RUS15", "Büntgen18_RUS17",
          "Büntgen18_MON09", "Büntgen18_RUS20",
          "Büntgen18_PAK04", "Büntgen18_SWE02",
          "Uusitalo18_Pine"]
lat_sharp = [48 + 43/60, 63 + 9/60,
             48 + 48/60, 37 + 27/60,
             30 + 20/60, 35 + 57/60,
            46 + 40/60, 58 +38/60,
            37 + 77/60, ]
events_sharp = ["Büntgen18_GER01", "Büntgen18_SWE05",
                "Büntgen18_GER07", "Büntgen18_CHN01",
                "Büntgen18_JAP01", "Büntgen18_USA02",
               "Büntgen18_MON03", "Büntgen18_USA11",
               "Büntgen18_USA18", ]
events = events_prolonged + events_sharp
lat = lat_prolonged + lat_sharp
cbm_models = ["Guttler15", "Buntgen18", "Brehm21",]
steady_state = []
for cbm_model in cbm_models:
    sf = fitting.SingleFitter(cbm_model, cbm_model=cbm_model)
    steady_state.append(sf.steady_state_production)
Buntgen18 = np.zeros((6000, len(events)))
for i, event in enumerate(events_prolonged):
    Buntgen18[:, i] = np.load("individual_chain/Prolonged_chain/{}_Buntgen18.npy".format(event))[:, 2]
    Buntgen18[:, i] = 10**Buntgen18[:, i]
for i, event in enumerate(events_sharp):
    Buntgen18[:, i + len(events_prolonged)] = np.load("individual_chain/Sharp_chain/{}_Buntgen18.npy".format(event))[:, 2]
    Buntgen18[:, i + len(events_prolonged)] = 10**Buntgen18[:, i + len(events_prolonged)]
Buntgen18 = Buntgen18 / steady_state[1]
def log_likelihood(theta, x, y, yerr):
    m, b, log_f = theta
    model = m * x + b
    sigma2 = yerr ** 2 + model ** 2 * np.exp(2 * log_f)
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))
def log_prior(theta):
    m, b, log_f = theta
    if -5.0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < log_f < 1.0:
        return 0.0
    return -np.inf
def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)

from scipy.optimize import minimize
np.random.seed(42)
f_true =0.1
x, y, yerr = np.array(lat), np.mean(Buntgen18, axis=0), np.std(Buntgen18, axis=0)
x0 = np.linspace(x.min(),x.max(),1000)

nll = lambda *args: -log_likelihood(*args)
initial = np.array([0.1, 3, np.log(f_true)]) + 0.1 * np.random.randn(3)
soln = minimize(nll, initial, args=(x, y, yerr))

import emcee
pos = soln.x + 1e-4 * np.random.randn(32, 3)
nwalkers, ndim = pos.shape

sampler = emcee.EnsembleSampler(
    nwalkers, ndim, log_probability, args=(x, y, yerr)
)
sampler.run_mcmc(pos, 5000, progress=True);
flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
nsamps = 1001
all_models = np.zeros((nsamps,len(x0)))

inds = np.random.randint(len(flat_samples), size=nsamps)
for j, ind in enumerate(inds):
    sample = flat_samples[ind]
    all_models[j,:] = np.dot(np.vander(x0, 2), sample[:2])
lower, mid, upper = np.percentile(all_models,(16,50,84),axis=0)

lines = np.percentile(flat_samples[:,0],(16,50,84))
print('%.1f +- %.1f' % (1e3*np.mean(flat_samples[:,0]), 1e3*np.std(flat_samples[:,0])))
plt.figure(figsize=(9, 6), dpi=80)

x0 = np.linspace(x.min()-5,x.max()+5,1000)

plt.errorbar(np.array(lat), np.mean(Buntgen18, axis=0), yerr=np.std(Buntgen18, axis=0), fmt="o", capsize=3)
plt.plot(x0,all_models.mean(axis=0),color='k')
plt.fill_between(x0,lower,upper,alpha=0.2,color='k')
plt.xlabel("Latitude (degrees N)", fontfamily="serif", fontsize=14)
plt.ylabel("Spike production in SS yr", fontfamily="serif", fontsize=14)
plt.text(31, 3.2, "JAP")
plt.text(66, 4.4, "LAP")
plt.text(39, 2.7, "CAL")
plt.text(51, 3.8, "ALT")
plt.text(46, 3.93, "GER")
plt.text(65, 3.93, "YAM")
plt.xlim(x0.min(),x0.max())
plt.savefig(snakemake.output[0], bbox_inches='tight')
