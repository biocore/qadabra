stylesheet = "config/qadabra.mplstyle"


rule plot_differentials:
    input:
        "results/{tool}/differentials.processed.tsv"
    output:
        "figures/{tool}_differentials.pdf"
    params:
        stylesheet
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot/plot_differentials.py"


rule plot_rank_correlation:
    input:
        "results/concatenated_differentials.tsv"
    output:
        "figures/spearman_heatmap.pdf"
    params:
        stylesheet
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/plot/plot_rank_correlations.py"


rule interactive:
    input:
        "results/concatenated_differentials.tsv"
    output:
        "results/qadabra.html"
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/interactive_app.py"
