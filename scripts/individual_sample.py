import numpy as np
from ticktack import fitting

sf = fitting.SingleFitter(snakemake.params.cbm_model, snakemake.params.cbm_model, hemisphere=snakemake.params.hemisphere)
sf.load_data(snakemake.input[0])
sf.compile_production_model(model=snakemake.params.production_model)
params = np.array([775, 1. / 12, 3., 81. / 12, 0.18])
low_bounds = np.array([775 - 5, 1 / 52., 0, 0., 0.])
up_bounds = np.array([775 + 5, 5., 11, 15., 2.])
chain = sf.MarkovChainSampler(params, sf.log_joint_likelihood, burnin=2000,
                                production=1000, args=(low_bounds, up_bounds))
np.save(snakemake.output[0], chain)
