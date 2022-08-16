library(biomformat)
library(dplyr)
library(Maaslin2)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

print("Loading table...")
table <- biomformat::read_biom(snakemake@input[["table"]])
table <- as.matrix(biomformat::biom_data(table))

print("Loading metadata")
metadata <- read.table(snakemake@input[["metadata"]], sep="\t", header=T,
                       row.names=1)

covariate <- snakemake@params[[1]][["factor_name"]]
target <- snakemake@params[[1]][["target_level"]]
reference <- snakemake@params[[1]][["reference_level"]]
confounders <- snakemake@params[[1]][["confounders"]]

print("Harmonizing table and metadata samples...")
samples <- colnames(table)
metadata <- subset(metadata, rownames(metadata) %in% samples)
metadata[[covariate]] <- as.factor(metadata[[covariate]])
metadata[[covariate]] <- relevel(metadata[[covariate]], reference)
sample_order <- row.names(metadata)

# Append F_ to features to avoid R renaming
row.names(table) <- paste0("F_", row.names(table))
table <- t(table[, sample_order])

fixed.effects <- c(covariate, confounders)
print(fixed.effects)

print("Running MaAsLin...")
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
row.names(results) <- gsub("^F_", "", results$feature)
results <- results %>% select(-c("feature"))

write.table(results, file=snakemake@output[["diff_file"]], sep="\t")
print("Saved differentials!")
