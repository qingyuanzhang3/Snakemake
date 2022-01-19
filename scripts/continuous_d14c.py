import matplotlib.pyplot as plt
import numpy as np
import ticktack
from ticktack import fitting
import os

models = snakemake.params.cbm_model
event = snakemake.params.event
colors = ["g", "b", "r", "c"]
size = 100
size2 = 30
fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 8), sharex=True, gridspec_kw={'height_ratios': [2, 1]})
min_prod = np.inf; max_prod = -np.inf
i = 0
for i in range(len(snakemake.input[1:])):
    model = models[i]
    cbm = ticktack.load_presaved_model(model, production_rate_units = 'atoms/cm^2/s')
    sf = fitting.SingleFitter(cbm, cbm_model=model)
    sf.load_data(snakemake.input[0])
    chain = np.load(snakemake.input[i+1])
    nwalkers = chain.shape[1] * 2
    sf.compile_production_model(model=snakemake.params.production_model)

    idx = np.random.randint(len(chain), size=size)
    params = np.zeros((size, chain[0].size))
    d14cs = np.zeros((size, sf.time_data_fine.size))
    params[:, :] = chain[idx]
    for j in range(size):
        dc14 = sf.dc14_fine(params=params[j, :])
        d14cs[j, :] = dc14

    for d14c in d14cs[:-1]:
        ax1.plot(sf.time_data_fine, d14c, alpha=0.05, color=colors[i])

    ax1.plot(sf.time_data_fine, d14cs[-1], color=colors[i], label=model)
    ax1.set_ylabel("$\Delta^{14}$C (‰)")
    fig.subplots_adjust(hspace=0.05)

    for param in params[:size2]:
        production_rate = sf.production(sf.time_data_fine, *param)
        ax2.plot(sf.time_data_fine, production_rate, alpha=0.2, color=colors[i])
        if np.max(production_rate) > max_prod:
            max_prod = np.max(production_rate)
        elif np.min(production_rate) < min_prod:
            min_prod = np.min(production_rate)
    i += 1

path = "data/" + event
file_names = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
for file in file_names:
    sf.load_data(path + '/' + file)
    ax1.errorbar(sf.time_data, sf.d14c_data, fmt="o", color="gray", yerr=sf.d14c_data_error, capsize=3, alpha=0.2)
    ax1.plot(sf.time_data, sf.d14c_data, "o",  color="gray",  alpha=0.2)
sf.load_data(snakemake.input[0])
ax1.errorbar(sf.time_data, sf.d14c_data, yerr=sf.d14c_data_error, fmt="ok", capsize=3,
             markersize=6.5, elinewidth=3, label="average $\Delta^{14}$C data")
ax1.legend();
ax2.set_ylim(-0.1, 20);
ax2.set_xlabel("Calendar Year (CE)");
ax2.set_xlim(sf.start-0.2, sf.end+0.2);
ax2.set_ylabel("Production rate ($cm^2s^{-1}$)");
ax2.legend(loc="upper left")
plt.suptitle(event);
plt.tight_layout();
plt.savefig(snakemake.output[0])
