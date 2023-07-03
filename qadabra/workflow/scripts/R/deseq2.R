library(biomformat)
library(DESeq2)

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

# Run DESeq2
print("Running DESeq2...")
dds <- DESeq2::DESeqDataSetFromMatrix(
    countData=table,
    colData=metadata,
    design=design.formula
)
dds.results <- DESeq2::DESeq(dds, sfType="poscounts")
saveRDS(dds.results, snakemake@output[[2]])
print("Saved RDS!")

# Create results table and modify result names
results <- DESeq2::results(
    dds.results,
    format="DataFrame",
    tidy=TRUE,
    cooksCutoff=FALSE,
    pAdjustMethod="BH",
    contrast=c(covariate, target, reference)
)
row.names(results) <- gsub("^F_", "", row.names(table))

# Save results to output file
write.table(results, file=snakemake@output[[1]], sep="\t")
print("Saved differentials!")
