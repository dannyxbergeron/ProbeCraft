rule fetch_clinvar:
    output:
        mutations="data/{gene}/clinvar_mutations.tsv",
        transcript_id="data/{gene}/principal_transcript.txt",
    params:
        gene="{gene}",
        api_key=config.get("ncbi_api_key", ""),
        molecular_consequences=config.get("molecular_consequences", ["missense variant"]),
        region=config.get("region", "CDS"),
    log: "logs/fetch_clinvar/{gene}.log"
    conda: "../envs/fetch.yaml"
    script: "../scripts/fetch_clinvar.py"


rule fetch_transcript:
    input:
        transcript_id="data/{gene}/principal_transcript.txt",
    output:
        transcript="data/{gene}/transcript.tsv",
        junctions="data/{gene}/junctions.tsv",
    params:
        api_key=config.get("ncbi_api_key", ""),
    log: "logs/fetch_transcript/{gene}.log"
    conda: "../envs/fetch.yaml"
    script: "../scripts/fetch_transcript.py"
