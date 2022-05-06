rule deseq2:
    input:
        table=config["table"],
        metadata=config["metadata"]
    output:
        "results/deseq2/ranks.tsv",
        "results/deseq2/results.rds"
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/deseq2.R"


rule ancombc:
    input:
        table=config["table"],
        metadata=config["metadata"]
    output:
        "results/ancombc/ranks.tsv",
        "results/ancombc/results.rds"
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/ancombc.R"
