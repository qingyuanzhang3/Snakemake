from ticktack import fitting
import numpy as np


mf, soln = fitting.fit_event(None, path=snakemake.input[0], production_model='control_points', sampler="optimisation",
                            cbm_model=snakemake.params.cbm_model,
                            hemisphere=snakemake.params.hemisphere,)
np.save(snakemake.output[0], soln[1])
