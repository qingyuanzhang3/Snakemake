import numpy as np
from ticktack import fitting

sf = fitting.SingleFitter(snakemake.params.cbm_model, snakemake.params.cbm_model)
sf.load_data(snakemake.input[0])
result = sf.MC_mean_std(iters=1000)
np.save(snakemake.output[0], result)
