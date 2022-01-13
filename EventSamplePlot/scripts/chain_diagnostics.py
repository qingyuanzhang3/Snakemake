import matplotlib.pyplot as plt
import numpy as np

chain = np.load(snakemake.input[0])
if chain.shape[0] == 4:
    labels = ["Start Date (yr)", "Duration (yr)", "phi (yr)", "Spike production"]
else:
    labels = ["Start Date (yr)", "Duration (yr)", "phi (yr)", "Spike production", "Solar amplitude"]

fig, axs = plt.subplots(2, 3, figsize=(18, 8), sharex=True)
axs = axs.flatten()
for i in range(chain.shape[1]):
    axs[i].plot(chain[:, i], 'b.', markersize=1, alpha=0.5)
    axs[i].set_title(labels[i])
    axs[i].get_xaxis().set_visible(False)
for i in range(chain.shape[1], 6):
    axs[i].set_axis_off()
plt.suptitle(snakemake.params.event + ' -- ' + snakemake.params.cbm_model, fontsize=25)
fig.savefig('diagnostics/' + snakemake.params.event + '_' + snakemake.params.cbm_model + '.jpg')
