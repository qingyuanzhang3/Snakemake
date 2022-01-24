import numpy as np
from ticktack import fitting
from astropy.table import Table

mf, sampler = fitting.fit_event(snakemake.params.year,
                                path=snakemake.input[0],
                                cbm_model=snakemake.params.cbm_model,
                                production_model=snakemake.params.production_model,
                                hemisphere=snakemake.params.hemisphere,
                                sampler="MCMC", burnin=2000, production=1000,
                                oversample=1008, burnin_time=2000)
np.save(snakemake.output[0], sampler)
