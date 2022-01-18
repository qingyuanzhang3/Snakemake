configfile: "config.yaml"
import os

def get_param_year(wildcards):
    return float(config["event"][wildcards.event])

def get_param_hem(wildcards):
    return config["hemisphere"][wildcards.event]

def get_sample_directory(wildcards):
    return "data/" + wildcards.event

def get_plot_dataset_input_directory(wildcards):
    return config["dataset"][wildcards.dataset]

production_model = "flexible_sinusoid"

rule all:
    input:
        expand("plots/posterior/{event}.jpg", event=config["event"]),
        expand("plots/diagnostics/{event}_{cbm_model}.jpg", event=config["event"], cbm_model=config["cbm_model"]),
        # expand("plots/diagnostics/{event}.jpg", event=config["event"]),
        # expand("plots/datasets/{dataset}.jpg", dataset=config["dataset"]),

rule sample:
    input:
        get_sample_directory
    output:
        "chain/{event}_{cbm_model}.npy"
    params:
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
        "plots/posterior/{event}.jpg"
    params:
        cbm_model = expand("{cbm_model}", cbm_model=config["cbm_model"]),
        event = "{event}"
    script:
        "scripts/plot_posterior.py"

rule plot_diagnostics:
    input:
        "chain/{event}_{cbm_model}.npy"
    output:
        "plots/diagnostics/{event}_{cbm_model}.jpg"
    params:
        event = "{event}",
        cbm_model = "{cbm_model}"
    script:
        "scripts/chain_diagnostics.py"

rule get_event_average:
    input:
        "data/{event}"
    output:
        "data/means/{event}.csv"
    params:
        event = "{event}"
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
        production_model = production_model
    script:
        "scripts/continuous_d14c.py"

rule plot_dataset:
    input:
        get_plot_dataset_input_directory
    output:
        "plots/datasets/{dataset}.jpg"
    script:
        "scripts/plot_dataset.py"
