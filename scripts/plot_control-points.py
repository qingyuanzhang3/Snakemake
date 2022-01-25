import matplotlib.pyplot as plt
import numpy as np
import ticktack
from ticktack import fitting
import os
from matplotlib.lines import Line2D
import matplotlib as mpl
import pickle
mpl.style.use('seaborn-colorblind')

models = snakemake.params.cbm_model
event = snakemake.params.event
colors = mpl.rcParams['axes.prop_cycle'].by_key()['color']
fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 8), sharex=True, gridspec_kw={'height_ratios': [2, 1]})

for i, model in enumerate(snakemake.input[1:]):
    model = models[i]
    cbm = ticktack.load_presaved_model(model, production_rate_units = 'atoms/cm^2/s')
    sf = fitting.SingleFitter(cbm, cbm_model=model, hemisphere=snakemake.params.hemisphere)
    sf.load_data(snakemake.input[0])
    with open(snakemake.input[i+1], 'rb') as f:
        soln = pickle.load(f)
    sf.compile_production_model(model="control_points")
    ax1.plot(sf.time_data_fine, sf.dc14_fine(soln.x), color=colors[i])
    ax1.set_ylabel("$\Delta^{14}$C (‰)")
    fig.subplots_adjust(hspace=0.05)
    ax2.plot(sf.control_points_time, soln.x, "o", color=colors[i])
    ax2.plot(sf.control_points_time, soln.x, color=colors[i])

path = "data/" + event
file_names = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
for file in file_names:
    sf.load_data(path + '/' + file)
    ax1.errorbar(sf.time_data, sf.d14c_data, fmt="o", color="gray", yerr=sf.d14c_data_error, capsize=3, alpha=0.2)
    ax1.plot(sf.time_data, sf.d14c_data, "o",  color="gray",  alpha=0.2)
sf.load_data(snakemake.input[0])
ax1.errorbar(sf.time_data, sf.d14c_data, yerr=sf.d14c_data_error, fmt="ok", capsize=3,
             markersize=6.5, elinewidth=3, label="average $\Delta^{14}$C")
custom_lines = [Line2D([0], [0], color=colors[0], lw=1.5, label=models[0]),
                Line2D([0], [0], color=colors[1], lw=1.5, label=models[1]),
                Line2D([0], [0], color=colors[2], lw=1.5, label=models[2]),
                Line2D([0], [0], color=colors[3], lw=1.5, label=models[3]),
                Line2D([0], [0], color="k", marker="o", lw=1.5, label="average $\Delta^{14}$C")]
ax1.legend(handles=custom_lines)
ax2.set_xlabel("Calendar Year (CE)");
ax2.set_xlim(sf.start-0.2, sf.end+0.2);
ax2.set_ylabel("Production rate (atoms/cm$^2$/s)");
ax2.legend(loc="upper left")
plt.suptitle(event);
plt.tight_layout();
plt.savefig(snakemake.output[0])
