import matplotlib.pyplot as plt
import numpy as np
from ticktack import fitting

chains = []
for i in range(len(snakemake.input)):
    chains.append(np.load(snakemake.input[i]))
if chains[-1].shape[1] == 4:
    labels =["Start Date (yr)", "Duration (yr)", "phi (yr)", "Spike production"]
else:
    labels = ["Start Date (yr)", "Duration (yr)", "phi (yr)", "Spike production", 'Solar amplitude']
cf = fitting.CarbonFitter()
fig = cf.plot_multiple_chains(chains, 8,
                        params_names=labels,
                        labels = snakemake.params.cbm_model
                        )
fig.savefig(snakemake.output[0])
