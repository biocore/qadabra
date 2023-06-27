library(biomformat)
library(corncob)
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
design.formula <- paste0("~", covariate)
if (length(confounders) != 0) {
    confounders_form = paste(confounders, collapse=" + ")
    design.formula <- paste0(design.formula, " + ", confounders_form)
}
design.formula <- as.formula(design.formula)
print(design.formula)

# Run Corncob
print("Running corncob...")
fit <- corncob::differentialTest(
    formula=design.formula,
    formula_null=~1,
    phi.formula=design.formula,
    phi.formula_null=design.formula,
    test="Wald",
    data=physeq,
    boot=F,
    filter_discriminant=F,
    full_output=T
)
saveRDS(fit, snakemake@output[[2]])
print("Saved RDS!")

# Create results table
taxa_names <- names(fit$p)
colname <- paste0("mu.", covariate, target)

print("Aggregating models...")
all_models <- fit$all_models
coefs <- c()
for (mod in all_models) {
    coef <- mod$coefficients[c(colname), "Estimate"]
    coefs <- c(coefs, coef)
}
coefs <- as.data.frame(coefs)
row.names(coefs) <- taxa_names
row.names(coefs) <- gsub("^F_", "", row.names(coefs))

pvals <- as.data.frame(fit$p)
adjusted_p_values <- p.adjust(fit$p, method = "BH")
pval_BH_adj <- as.data.frame(adjusted_p_values)

results_all <- data.frame(coefs=coefs, pvals=pvals, pval_BH_adj=pval_BH_adj)

# Save results to output file
write.table(results_all, snakemake@output[[1]], sep="\t")
print("Saved differentials!")
