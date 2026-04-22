# ==============================================================================
# Script Name: 06_wgcna.R
#
# Description:
#   This script performs Weighted Gene Co-expression Network Analysis (WGCNA).
#   It identifies clusters (modules) of highly correlated genes and relates 
#   these modules to external sample traits (e.g., Control vs Treatment, clinical data).
#
#   Tasks performed:
#   1. Data Preparation: VST normalization of raw counts and variance filtering 
#      (removes low-variance genes to reduce computational load).
#   2. Outlier Detection: Sample clustering to check for outliers.
#   3. Soft Thresholding: Calculates the optimal power for scale-free topology.
#   4. Network Construction: One-step module detection using blockwiseModules.
#   5. Visualization: Dendrogram of genes and their assigned module colors.
#
# Inputs:
#   - Raw count matrix (e.g., counts.txt from step 03).
#   - Trait/Phenotype data file (CSV or TXT) matching your samples.
#
# Outputs:
#   - Soft threshold evaluation plots.
#   - Module assignment dendrogram plot.
#   - Data frames containing module assignments and eigengenes.
#
# Usage:
#   1. Create a trait data file (sample names as row names, traits as columns).
#   2. Modify the "USER CONFIGURATION" section below.
#   3. Run the script: Rscript scripts/06_wgcna.R
# ==============================================================================

# --- USER CONFIGURATION START: Modify these parameters for your specific data ---
COUNTS_FILE <- "./results/quantification/counts.txt"

# You MUST create a CSV file containing your sample traits.
# Row names must exactly match the sample names in your count matrix.
# Columns should be numerical traits or binary indicators (e.g., 0 for Control, 1 for Treatment).
TRAIT_FILE <- "./data/sample_traits.csv" 

OUTPUT_DIR <- "./results/wgcna"

# WGCNA Parameters
SOFT_POWER <- 15          # Adjust this based on the Soft Threshold plot output if necessary
MIN_MODULE_SIZE <- 30
MERGE_CUT_HEIGHT <- 0.25  # Threshold to merge similar modules
# --- USER CONFIGURATION END ---

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

suppressPackageStartupMessages({
  library(WGCNA)
  library(DESeq2)
  library(ggplot2)
})

# Allow multi-threading within WGCNA
allowWGCNAThreads()

cat("-> Loading and preprocessing data...\n")
counts <- read.table(COUNTS_FILE, header = TRUE, row.names = 1)
counts <- counts[, -c(1:5)] # Remove featureCounts annotation columns

# We need a dummy colData to run VST using DESeq2
dummy_colData <- data.frame(row.names = colnames(counts), condition = rep("A", ncol(counts)))
dds <- DESeqDataSetFromMatrix(countData = counts, colData = dummy_colData, design = ~ 1)

cat("-> Performing VST transformation...\n")
vsd <- vst(dds, blind = FALSE)
expr_data_unfiltered <- assay(vsd)

cat("-> Filtering genes by variance (Keeping top 25%)...\n")
gene_variance <- apply(expr_data_unfiltered, 1, var)
variance_cutoff <- quantile(gene_variance, probs = 0.75) 
keep_genes <- gene_variance > variance_cutoff
expr_data <- t(expr_data_unfiltered[keep_genes, ]) # WGCNA requires genes as columns, samples as rows

cat("-> Loading trait data...\n")
# Note: Ensure your trait file has rownames matching colnames of counts
if (file.exists(TRAIT_FILE)) {
  trait_data <- read.csv(TRAIT_FILE, row.names = 1)
  # Ensure sample alignment
  common_samples <- intersect(rownames(expr_data), rownames(trait_data))
  expr_data <- expr_data[common_samples, ]
  trait_data <- trait_data[common_samples, , drop = FALSE]
} else {
  cat("   WARNING: Trait file not found. Proceeding with network construction only.\n")
}

cat("-> Step 1: Checking for sample outliers...\n")
sample_tree <- hclust(dist(expr_data), method = "average")
pdf(file.path(OUTPUT_DIR, "Sample_Clustering.pdf"), width = 8, height = 6)
plot(sample_tree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

cat("-> Step 2: Evaluating Soft Thresholding Powers...\n")
powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))
sft <- pickSoftThreshold(expr_data, powerVector = powers, verbose = 5, networkType = "signed")

pdf(file.path(OUTPUT_DIR, "Soft_Threshold_Plots.pdf"), width = 12, height = 5)
par(mfrow = c(1, 2))
# Scale-free topology fit index
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.90, col = "red")

# Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", 
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

cat("-> Step 3: Constructing Network and Identifying Modules...\n")
# Using the pre-defined SOFT_POWER. In a real analysis, you might pause and look at the plots first.
net <- blockwiseModules(expr_data, 
                        power = SOFT_POWER,
                        TOMType = "signed", 
                        minModuleSize = MIN_MODULE_SIZE,
                        reassignThreshold = 0, 
                        mergeCutHeight = MERGE_CUT_HEIGHT,
                        numericLabels = TRUE, 
                        saveTOMs = TRUE,
                        saveTOMFileBase = file.path(OUTPUT_DIR, "WGCNA_TOM"),
                        verbose = 3)

cat("-> Step 4: Visualizing Module Dendrogram...\n")
moduleColors <- labels2colors(net$colors)

pdf(file.path(OUTPUT_DIR, "Module_Dendrogram.pdf"), width = 10, height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Cluster Dendrogram and Module Colors")
dev.off()

# Save the module assignments
gene_module_info <- data.frame(Gene = colnames(expr_data), ModuleColor = moduleColors)
write.csv(gene_module_info, file.path(OUTPUT_DIR, "Gene_Module_Assignments.csv"), row.names = FALSE)

cat("=== WGCNA Analysis Completed! Results saved in", OUTPUT_DIR, "===\n")
