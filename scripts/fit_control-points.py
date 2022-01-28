from astropy.table import Table
from ticktack import fitting
import pickle


mf, soln = fitting.fit_event(None, path=snakemake.input[0], production_model='control_points', sampler="optimisation",
                            cbm_model=snakemake.params.cbm_model,
                            hemisphere=snakemake.params.hemisphere,)
with open(snakemake.output[0], 'wb') as file:
    pickle.dump(soln, file, protocol=pickle.HIGHEST_PROTOCOL)
