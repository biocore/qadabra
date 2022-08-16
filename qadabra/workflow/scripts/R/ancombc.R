library(biomformat)
library(ANCOMBC)
library(phyloseq)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

print("Loading table...")
table <- biomformat::read_biom(snakemake@input[["table"]])
table <- as.matrix(biomformat::biom_data(table))

print("Loading metadata...")
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
table <- table[, sample_order]
# Append F_ to features to avoid R renaming
row.names(table) <- paste0("F_", row.names(table))

print("Converting to phyloseq...")
taxa <- phyloseq::otu_table(table, taxa_are_rows=T)
meta <- phyloseq::sample_data(metadata)
physeq <- phyloseq::phyloseq(taxa, meta)

print("Creating design formula...")
design.formula <- covariate
if (length(confounders) != 0) {
    confounders_form = paste(confounders, collapse=" + ")
    design.formula <- paste0(design.formula, " + ", confounders_form)
}
print(design.formula)

print("Running ANCOMBC...")
ancombc.results <- ANCOMBC::ancombc(phyloseq=physeq, formula=design.formula,
                                    zero_cut=1.0)
saveRDS(ancombc.results, snakemake@output[[2]])
print("Saved RDS!")
results <- ancombc.results$res$beta
row.names(results) <- gsub("^F_", "", row.names(results))

write.table(results, file=snakemake@output[[1]], sep="\t")
print("Saved differentials!")
