rule generate_report:
    input:
        mutations="data/{gene}/clinvar_mutations.tsv",
        transcript="data/{gene}/transcript.tsv",
        junctions="data/{gene}/junctions.tsv",
        probes="data/{gene}/probes.tsv",
    output:
        report="results/{gene}/{gene}_report.html",
        mutations_file="results/{gene}/{gene}_mutations.tsv",
    params:
        gene="{gene}",
    log: "logs/generate_report/{gene}.log"
    conda: "../envs/report.yaml"
    script: "../scripts/generate_report.py"


rule export_probes:
    input:
        "data/{gene}/probes.tsv",
    output:
        "results/{gene}/probes.tsv",
    shell:
        "cp {input} {output}"
