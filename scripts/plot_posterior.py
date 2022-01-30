import matplotlib.pyplot as plt
import numpy as np
from ticktack import fitting
import matplotlib as mpl
mpl.style.use('seaborn-colorblind')

if snakemake.params.year > 0:
    start_date = "start date ({})".format("yr - " + str(int(snakemake.params.year) - 5) + "CE")
else:
    start_date = "start date ({})".format(str(-int(snakemake.params.year) + 5) + "BCE - yr")
chains = []
for i in range(len(snakemake.input)):
    chain = np.load(snakemake.input[i])
    chains.append(chain)
if snakemake.params.production_model == "simple_sinusoid":
    labels = [start_date, "duration (yr)", "$\phi$ (yr)", "spike production (atoms/cm$^2$ yr/s)"]
    idx = 0
elif snakemake.params.production_model == "flexible_sinusoid":
    labels = [start_date, "duration (yr)", "$\phi$ (yr)", "spike production (atoms/cm$^2$ yr/s)", "solar amplitude (atoms/cm$^2$/s)"]
    idx = 0
elif snakemake.params.production_model == "flexible_sinusoid_affine_variant":
    labels = ["gradient (atoms/cm$^2$/year$^2$)", start_date, "duration (yr)", "$\phi$ (yr)", "spike production (atoms/cm$^2$ yr/s)", "solar amplitude (atoms/cm$^2$/s)"]
    idx = 1
else:
    labels = None

for chain in chains:
    chain[:, idx] = chain[:, idx] - (snakemake.params.year - 5)

colors = mpl.rcParams['axes.prop_cycle'].by_key()['color']
cf = fitting.CarbonFitter()
fig = cf.plot_multiple_chains(chains, chain.shape[1] * 2,
                        params_labels=labels,
                        labels = snakemake.params.cbm_label,
                        label_font_size=6.5,
                        tick_font_size=8, colors=colors
                        )
plt.suptitle(snakemake.params.event_label, fontsize=25)
fig.savefig(snakemake.output[0])
