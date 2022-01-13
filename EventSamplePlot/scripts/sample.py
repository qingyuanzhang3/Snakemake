import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
from ticktack import fitting

mf, sampler = fitting.fit_event(snakemake.params.year,
                                path=snakemake.input[0],
                                cbm_model=snakemake.params.cbm_model,
                                production_model="flexible_sinusoid", hemisphere='north',
                                sampler="MCMC", burnin=2000, production=2000,
                                oversample=108)
np.save(snakemake.output[0], sampler)
