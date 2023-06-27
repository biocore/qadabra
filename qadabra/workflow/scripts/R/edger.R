library(biomformat)
library(edgeR)

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

# Create the design formula
print("Creating design formula...")
design.formula <- paste0("~", covariate)
if (length(confounders) != 0) {
    confounders_form = paste(confounders, collapse=" + ")
    design.formula <- paste0(design.formula, " + ", confounders_form)
}
design.formula <- as.formula(design.formula)
print(design.formula)
mm <- model.matrix(design.formula, metadata)

# Run edgeR
print("Running edgeR...")
d <- edgeR::DGEList(counts=table, group=metadata[[covariate]])
d <- edgeR::calcNormFactors(d)
d <- edgeR::estimateDisp(d, design=mm)
fit <- edgeR::glmQLFit(d, mm, robust=T)
saveRDS(fit, snakemake@output[[2]])
print("Saved RDS!")

# Obtain coefficients using glmQLFit
coeff_results <- fit$coefficients
row.names(coeff_results) <- gsub("^F_", "", row.names(coeff_results))

# Obtain p-values and corrected p-values using glmQLFTest
res <- edgeR::glmQLFTest(fit, coef=2)
pval_results <- res$table
adjusted_p_values <- p.adjust(pval_results$PValue, method = "BH")
pval_results$PValue_BH_adj <- adjusted_p_values
row.names(pval_results) <- gsub("^F_", "", row.names(pval_results))

# Create results table
results_all <- cbind(coeff_results, pval_results[match(rownames(coeff_results), rownames(pval_results)),])

# Save results to output file
write.table(results_all, file=snakemake@output[[1]], sep="\t")
print("Saved differentials!")
