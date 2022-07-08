library(biomformat)
library(dplyr)
library(Maaslin2)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

table <- biomformat::read_biom(snakemake@input[["table"]])
table <- as.matrix(biomformat::biom_data(table))

metadata <- read.table(snakemake@input[["metadata"]], sep="\t", header=T,
                       row.names=1)

covariate <- snakemake@config[["model"]][["covariate"]]
target <- snakemake@config[["model"]][["target"]]
reference <- snakemake@config[["model"]][["reference"]]
confounders <- snakemake@config[["model"]][["confounders"]]

samples <- colnames(table)
metadata <- subset(metadata, rownames(metadata) %in% samples)
metadata[[covariate]] <- as.factor(metadata[[covariate]])
metadata[[covariate]] <- relevel(metadata[[covariate]], reference)
sample_order <- row.names(metadata)
table <- t(table[, sample_order])

fixed.effects <- c(covariate, confounders)

fit.data <- Maaslin2::Maaslin2(
    input_data=table,
    input_metadata=metadata,
    output=snakemake@output[["out_dir"]],
    fixed_effects=fixed.effects,
    min_prevalence=0,
    min_abundance=0,
    plot_scatter=F,
    plot_heatmap=F
)
results <- fit.data$results %>% dplyr::filter(metadata==covariate)
row.names(results) <- results$feature
results <- results %>% select(-c("feature"))

write.table(results, file=snakemake@output[["diff_file"]], sep="\t")
