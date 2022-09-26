import numpy as np
from ticktack import fitting

if snakemake.params.event == "5259BCE":
    params = np.array([-9.87302625e-03, -5.25961987e+03, -7.53614127e-01,  1.43895132e+00,
        8.53135328e-01, -1.22434240e+00])
elif snakemake.params.event == "775AD-Prolonged":
    params = np.array([-1.74123467e-03, 7.73721759e+02, 2.52047467e-01, 2.27947073e+00,
                       7.87020467e-01, -1.51734656e+00])
elif snakemake.params.event == "993AD":
    params = np.array([-5.32783967e-04,  9.92948483e+02,  np.log10(1.12824990e+00),
    6.04790666e+00,  np.log10(3.10722315e+00),  np.log10(3.82991199e-02)])
else:
    params = None
mf, sampler = fitting.fit_event(snakemake.params.year,
                                path=snakemake.input[0],
                                cbm_model=snakemake.params.cbm_model,
                                production_model=snakemake.params.production_model,
                                hemisphere=snakemake.params.hemisphere,
                                params=params,
                                sampler="MCMC", burnin=1000, production=1000,
                                oversample=1008, burnin_time=2000, verbose=True)
np.save(snakemake.output[0], sampler)
