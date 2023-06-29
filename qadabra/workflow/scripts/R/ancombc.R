library(biomformat)
library(ANCOMBC)
library(phyloseq)

# Set logging information
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# Load the input table
print("Loading table...")
table <- biomformat::read_biom(snakemake@input[["table"]])
table <- as.matrix(biomformat::biom_data(table))

# Load the metadata
print("Loading metadata...")
metadata <- read.table(snakemake@input[["metadata"]], sep="\t", header=T,
                       row.names=1)

covariate <- snakemake@params[[1]][["factor_name"]]
target <- snakemake@params[[1]][["target_level"]]
reference <- snakemake@params[[1]][["reference_level"]]
confounders <- snakemake@params[[1]][["confounders"]]

# Harmonize table and metadata samples
print("Harmonizing table and metadata samples...")
samples <- colnames(table)
metadata <- subset(metadata, rownames(metadata) %in% samples)
metadata[[covariate]] <- as.factor(metadata[[covariate]])
metadata[[covariate]] <- relevel(metadata[[covariate]], reference)
sample_order <- row.names(metadata)
table <- table[, sample_order]
row.names(table) <- paste0("F_", row.names(table)) # Append F_ to features to avoid R renaming

# Convert to phyloseq object
print("Converting to phyloseq...")
taxa <- phyloseq::otu_table(table, taxa_are_rows=T)
meta <- phyloseq::sample_data(metadata)
physeq <- phyloseq::phyloseq(taxa, meta)

# Create the design formula
print("Creating design formula...")
design.formula <- covariate
if (length(confounders) != 0) {
    confounders_form = paste(confounders, collapse=" + ")
    design.formula <- paste0(design.formula, " + ", confounders_form)
}

# Run ANCOMBC
print("Running ANCOMBC...")
ancombc.results <- ANCOMBC::ancombc(phyloseq=physeq, formula=design.formula, zero_cut=1.0, p_adj_method = "BH",)
saveRDS(ancombc.results, snakemake@output[[2]])
print("Saved RDS!")

# Access coefficients and p-values from the ANCOMBC results
coef_col <- paste("coefs", covariate, target, sep = ".")
pval_col <- paste("pvals", covariate, target, sep = ".")
coefs <- ancombc.results$res$beta
pvals <- ancombc.results$res$p_val
# qvals <- ancombc.results$res$q_val

# Modify result names
row.names(coefs) <- gsub("^F_", "", row.names(coefs))
row.names(pvals) <- gsub("^F_", "", row.names(pvals))
# row.names(qvals) <- gsub("^F_", "", row.names(qvals))

results_all <- data.frame(coefs=coefs, pvals=pvals)
colnames(results_all) <- c("coefs", "pvals")

# Save results to output file
write.table(results_all, file=snakemake@output[[1]], sep="\t")
print("Saved differentials!")
