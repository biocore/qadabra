stylesheet = "config/qadabra.mplstyle"


rule rank_correlation:
    input:
        "results/concatenated_differentials.tsv"
    output:
        "results/figures/spearman_heatmap.pdf"
    params:
        stylesheet
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/rank_correlations.py"
