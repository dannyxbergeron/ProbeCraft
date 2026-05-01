rule design_probes:
    input:
        mutations="data/{gene}/clinvar_mutations.tsv",
        transcript="data/{gene}/transcript.tsv",
        junctions="data/{gene}/junctions.tsv",
    output:
        "data/{gene}/probes.tsv",
    params:
        gene="{gene}",
        target_length=config["probe_params"]["target_length"],
        min_length=config["probe_params"]["min_length"],
        max_length=config["probe_params"]["max_length"],
        junction_buffer=config["probe_params"]["junction_buffer"],
        scan_offset=config["probe_params"]["scan_offset"],
        region=config.get("region", "CDS"),
    log: "logs/design_probes/{gene}.log"
    conda: "../envs/design.yaml"
    script: "../scripts/design_probes.py"
