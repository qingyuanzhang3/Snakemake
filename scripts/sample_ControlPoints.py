import numpy as np
from ticktack import fitting

mf, sampler = fitting.fit_event(None,
                                path=snakemake.input[0],
                                cbm_model=snakemake.params.cbm_model,
                                production_model="control_points",
                                hemisphere=snakemake.params.hemisphere,
                                sampler="MCMC", burnin=1000, production=1000,
                                oversample=1008, burnin_time=2000)
np.save(snakemake.output[0], sampler)
