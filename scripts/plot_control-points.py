from ticktack import fitting

fitting.plot_ControlPoints(average_path=snakemake.input[0], soln_path=snakemake.input[1:],
                            cbm_models=snakemake.params.cbm_model,
                            hemisphere=snakemake.params.hemisphere,
                            directory_path="data/" + snakemake.params.event,
                            savefig_path=snakemake.output[0], title=snakemake.params.event)
