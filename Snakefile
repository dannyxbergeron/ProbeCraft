configfile: "config.yaml"

include: "rules/fetch_data.smk"
include: "rules/probe_design.smk"
include: "rules/report.smk"

rule all:
    input:
        expand("results/{gene}/{gene}_report.html", gene=config["genes"]),
        expand("results/{gene}/probes.tsv", gene=config["genes"]),
        expand("results/{gene}/{gene}_mutations.tsv", gene=config["genes"]),
