import numpy as np
from ticktack import fitting

soln = np.load(snakemake.input[1], allow_pickle=True)
mf, sampler = fitting.fit_event(None,
                                path=snakemake.input[0],
                                cbm_model=snakemake.params.cbm_model,
                                production_model="control_points",
                                params=soln,
                                hemisphere=snakemake.params.hemisphere,
                                sampler="MCMC", burnin=100, production=500,
                                oversample=1008, burnin_time=2000)
np.save(snakemake.output[0], sampler)
