rule deseq2:
    input:
        table=config["table"],
        metadata=config["metadata"]
    output:
        "results/deseq2.tsv"
    conda:
        "../envs/qadabra-da-R.yml"
    script:
        "../scripts/deseq2.R"


rule ancombc:
    input:
        table=config["table"],
        metadata=config["metadata"]
    output:
        "results/ancombc.tsv"
    conda:
        "../envs/qadabra-da-R.yml"
    script:
        "../scripts/ancombc.R"
