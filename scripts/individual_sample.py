import numpy as np
from ticktack import fitting
import jax.numpy as jnp
from jax import jit

sf = fitting.SingleFitter(snakemake.params.cbm_model, snakemake.params.cbm_model, hemisphere=snakemake.params.hemisphere)
sf.load_data(snakemake.input[0])

if snakemake.params.production_model == "simple_sinusoid_sharp":
    chain = np.load("chain/775AD-Sharp_{}.npy".format(snakemake.params.cbm_model))
    mle = sf.chain_summary(chain, 12, mle=True)
    @jit
    def simple_sinusoid_sharp(t, *args):
        start_time, log_duration, log_area = jnp.array(list(args)).reshape(-1)
        duration, area = 10**log_duration, 10**log_area
        height = sf.super_gaussian(t, start_time, duration, area)
        production = sf.steady_state_production + mle[-1] * sf.steady_state_production * jnp.sin(
            2 * np.pi / 11 * t + mle[3] * 2 * np.pi / 11) + height
        return production
    model = simple_sinusoid_sharp
else:
    chain = np.load("chain/775AD-Prolonged_{}.npy".format(snakemake.params.cbm_model))
    mle = sf.chain_summary(chain, 12, mle=True)
    @jit
    def simple_sinusoid_prolonged(t, *args):
        start_time, log_duration, log_area = jnp.array(list(args)).reshape(-1)
        duration, area = 10 ** log_duration, 10 ** log_area
        height = sf.super_gaussian(t, start_time, duration, area)
        production = sf.steady_state_production + mle[-1] * sf.steady_state_production * jnp.sin(
            2 * np.pi / 11 * t + mle[3] * 2 * np.pi / 11) + height
        return production
    model = simple_sinusoid_prolonged

sf.compile_production_model(model=model)
params = np.array([775, np.log10(1. / 12), np.log10(81./ 12)])
low_bounds = np.array([775 - 5, np.log10(1 / 52.), -2])
up_bounds = np.array([775 + 5, np.log10(5.), 1.5])
chain = sf.MarkovChainSampler(params, sf.log_joint_likelihood, burnin=1000,
                                production=1000, args=(low_bounds, up_bounds))
np.save(snakemake.output[0], chain)
