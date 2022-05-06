library(biomformat)
library(ANCOMBC)
library(phyloseq)

table <- biomformat::read_biom(snakemake@input[["table"]])
table <- as.matrix(biomformat::biom_data(table))

metadata <- read.table(snakemake@input[["metadata"]], sep="\t", header=T,
                       row.names=1)

covariate <- snakemake@config[["model"]][["covariate"]]
target <- snakemake@config[["model"]][["target"]]
reference <- snakemake@config[["model"]][["reference"]]

samples <- colnames(table)
metadata <- subset(metadata, rownames(metadata) %in% samples)
metadata[[covariate]] <- as.factor(metadata[[covariate]])
metadata[[covariate]] <- relevel(metadata[[covariate]], reference)
sample_order <- row.names(metadata)
table <- table[, sample_order]

taxa <- phyloseq::otu_table(table, taxa_are_rows=T)
meta <- phyloseq::sample_data(metadata)
physeq <- phyloseq::phyloseq(taxa, meta)

design.formula <- snakemake@config[["model"]][["covariate"]]
ancombc.results <- ANCOMBC::ancombc(phyloseq=physeq, formula=design.formula,
                                    zero_cut=1.0)
results <- ancombc.results$res$beta

write.table(results, file=snakemake@output[[1]], sep="\t")
