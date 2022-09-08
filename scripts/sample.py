import numpy as np
from ticktack import fitting

if snakemake.params.event == "5259BCE":
    params = np.array([-8.28506168e-03, -5.25965399e+03,  np.log10(4.66141251e-01),
    3.09932956e+00, np.log10(6.24373446e+00),  np.log10(2.51351680e-02)])
else:
    params = None
mf, sampler = fitting.fit_event(snakemake.params.year,
                                path=snakemake.input[0],
                                cbm_model=snakemake.params.cbm_model,
                                production_model=snakemake.params.production_model,
                                hemisphere=snakemake.params.hemisphere,
                                params=params,
                                sampler="MCMC", burnin=1000, production=1000,
                                oversample=1008, burnin_time=2000)
np.save(snakemake.output[0], sampler)
