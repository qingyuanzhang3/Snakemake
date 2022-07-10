configfile: "config.yaml"
import os

def get_param_year(wildcards):
    return float(config["event"][wildcards.event])

def get_param_event_label(wildcards):
    return config["event_label"][wildcards.event]

def get_param_hem(wildcards):
    return config["hemisphere"][wildcards.event]

def get_sample_directory(wildcards):
    return "data/" + wildcards.event

def get_cp_directory(wildcards):
    return "data-CP/" + wildcards.event

def get_sample_path_prolonged(wildcards):
    return "data/775AD-Prolonged/" + wildcards.event + ".csv"

def get_sample_path_sharp(wildcards):
    return "data/775AD-Sharp/" + wildcards.event + ".csv"

production_model = "flexible_sinusoid_affine_variant"

rule all:
    input:
        expand("plots/diagnostics/{event}.jpg", event=config["event"]), # continuous sample plot
        expand("data/means/{averages}.csv", averages=config["averages"]), # supplementary mean csv
        expand("data-CP/means/{averages}.csv", averages=config["averages"]), # supplementary mean csv
        expand("non-parametric/chain/{event}_{cbm_model}.npy", event=config["event"], cbm_model=config["cbm_model"]), # control-point chain
        expand("non-parametric/solutions/{event}_{cbm_model}.npy", event=config["event"], cbm_model=config["cbm_model"]), # control-point solution
        expand("non-parametric/solver/{event}_{cbm_model}.npy", event=config["event"], cbm_model=config["cbm_model"]), # inverse solver
        expand("individual_chain/Prolonged_chain/{event}_{cbm_model}.npy", event=config["Prolonged"], cbm_model=config["cbm_model"]),
        expand("individual_chain/Sharp_chain/{event}_{cbm_model}.npy", event=config["Sharp"], cbm_model=config["cbm_model"]),
        "plots/f4/event_duration.pdf",
        "plots/f5/solar_density.pdf",
        "plots/f3/production_density.pdf",
        "plots/f2/non-parametric.pdf",
        "plots/f1/parametric-fit.pdf",
        expand("plots/f10-16/{event}.pdf", event=config["event"]), # posterior
        # ## expand("plots/diagnostics/{event}_{cbm_model}.jpg", event=config["event"], cbm_model=config["cbm_model"]), # chain plot
        # ## expand("plots/control-points/{event}.jpg", event=config["event"]), # control-points plot
        ## expand("data/means/{event}.csv", event=config["event"]), # supplementary mean csv
        ## expand("data-CP/means/{event}.csv", event=config["event"]), # supplementary mean csv

rule sample_sharp:
    input:
        get_sample_path_sharp,
        expand("chain/775AD-Sharp_{cbm_model}.npy", cbm_model=config["cbm_model"])
    output:
        "individual_chain/Sharp_chain/{event}_{cbm_model}.npy"
    params:
        year = 775,
        cbm_model = "{cbm_model}",
        production_model = "simple_sinusoid_sharp",
        hemisphere = "north"
    script:
        "scripts/individual_sample.py"

rule sample_prolonged:
    input:
        get_sample_path_prolonged,
        expand("chain/775AD-Prolonged_{cbm_model}.npy", cbm_model=config["cbm_model"])
    output:
        "individual_chain/Prolonged_chain/{event}_{cbm_model}.npy"
    params:
        year = 775,
        cbm_model = "{cbm_model}",
        production_model = "simple_sinusoid_prolonged",
        hemisphere = "north"
    script:
        "scripts/individual_sample.py"

rule sample:
    input:
        get_sample_directory
    output:
        "chain/{event}_{cbm_model}.npy"
    params:
        event = "{event}",
        year = get_param_year,
        cbm_model = "{cbm_model}",
        hemisphere = get_param_hem,
        production_model = production_model
    script:
        "scripts/sample.py"

rule plot_posterior:
    input:
        expand("chain/{event}_{cbm_model}.npy", event="{event}", cbm_model=config["cbm_model"])
    output:
        "plots/f10-16/{event}.pdf"
    params:
        year = get_param_year,
        cbm_label = expand("{cbm_label}", cbm_label=config["cbm_label"].values()),
        event_label = get_param_event_label,
        production_model = production_model,
    script:
        "scripts/plot_posterior.py"

rule plot_duration:
    input:
        expand("chain/{event}_{cbm_model}.npy", event=config["event"],
        cbm_model=config["cbm_model"])
    output:
        "plots/f4/event_duration.pdf"
    script:
        "scripts/plot_duration.py"

rule plot_production:
    input:
        expand("chain/{event}_{cbm_model}.npy", event=config["event"],
        cbm_model=config["cbm_model"])
    output:
        "plots/f3/production_density.pdf"
    script:
        "scripts/plot_production.py"

rule plot_non_parametric_fit:
    input:
        expand("non-parametric/solutions/{event}_{cbm_model}.npy",
        event=config["event"], cbm_model=config["cbm_model"]),
        expand("non-parametric/solver/{event}_{cbm_model}.npy",
        event=config["event"], cbm_model=config["cbm_model"]),
        expand("non-parametric/chain/{event}_{cbm_model}.npy",
        event=config["event"], cbm_model=config["cbm_model"]),
    output:
        "plots/f2/non-parametric.pdf"
    script:
        "scripts/plot_non-param-fit.py"

rule plot_parametric_fit:
    input:
        expand("chain/{event}_{cbm_model}.npy", event=config["event"],
        cbm_model=config["cbm_model"]),
        expand("data/means/{event}.csv", event=config["event"])
    output:
        "plots/f1/parametric-fit.pdf"
    script:
        "scripts/plot_param-fit.py"


rule plot_solar:
    input:
        expand("chain/{event}_{cbm_model}.npy", event=config["event"],
        cbm_model=config["cbm_model"])
    output:
        "plots/f5/solar_density.pdf"
    script:
        "scripts/plot_solar.py"

rule plot_diagnostics:
    input:
        "chain/{event}_{cbm_model}.npy"
    output:
        "plots/diagnostics/{event}_{cbm_model}.jpg"
    params:
        event = "{event}",
        cbm_model = "{cbm_model}",
        production_model = production_model,
    script:
        "scripts/chain_diagnostics.py"

rule get_event_average:
    input:
        "data/{event}"
    output:
        "data/means/{event}.csv"
    script:
        "scripts/event_average.py"

rule get_event_cp_average:
    input:
        "data-CP/{event}"
    output:
        "data-CP/means/{event}.csv"
    script:
        "scripts/event_average.py"

rule plot_continuous_d14c:
    input:
        "data/means/{event}.csv",
        expand("chain/{event}_{cbm_model}.npy", event="{event}", cbm_model=config["cbm_model"])
    output:
        "plots/diagnostics/{event}.jpg"
    params:
        event = "{event}",
        cbm_model = expand("{cbm_model}", cbm_model=config["cbm_model"]),
        cbm_label = expand("{cbm_label}", cbm_label=config["cbm_label"]),
        production_model = production_model,
        event_label = get_param_event_label,
        hemisphere = get_param_hem
    script:
        "scripts/plot_continuous_samples.py"

rule fit_control_points:
    input:
        get_cp_directory
    output:
        "non-parametric/solutions/{event}_{cbm_model}.npy"
    params:
        cbm_model = "{cbm_model}",
        hemisphere = get_param_hem,
    script:
        "scripts/fit_control-points.py"

rule plot_control_points:
    input:
        "data-CP/means/{event}.csv",
        expand("non-parametric/solutions/{event}_{cbm_model}.npy", event="{event}", cbm_model=config["cbm_model"])
    output:
        "plots/control-points/{event}.jpg"
    params:
        event = "{event}",
        cbm_model = expand("{cbm_model}", cbm_model=config["cbm_model"]),
        cbm_label = expand("{cbm_label}", cbm_label=config["cbm_label"]),
        hemisphere = get_param_hem,
    script:
        "scripts/plot_control-points.py"

rule sample_ControlPoints_uncertainty:
    input:
        get_cp_directory,
        "non-parametric/solutions/{event}_{cbm_model}.npy"
    output:
        "non-parametric/chain/{event}_{cbm_model}.npy"
    params:
        cbm_model = "{cbm_model}",
        hemisphere = get_param_hem,
    script:
        "scripts/sample_ControlPoints.py"

rule sample_inverse_solver:
    input:
        "data-CP/means/{event}.csv"
    output:
        "non-parametric/solver/{event}_{cbm_model}.npy"
    params:
        cbm_model = "{cbm_model}",
    script:
        "scripts/sample_inverse_solver.py"
