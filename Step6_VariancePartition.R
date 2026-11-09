# Step6_VariancePartition.R
# Microglia pseudobulk (one profile per sample): fraction of gene expression variance explained by
# genotype, sex and genotype x sex, after removing sequencing run as a technical batch (Fig 2G).
# Input:  output/step3_clustering/Plcg2_annotated.rds
# Output: output/step6_variancePartition/variancePartition_Microglia.csv
# Run:    Rscript Step6_VariancePartition.R

suppressPackageStartupMessages({
  library(Seurat)
  library(variancePartition)
  library(edgeR)
  library(limma)
})

# ---- Config ------------------------------------------------------------------
input_rds       <- "output/step3_clustering/Plcg2_annotated.rds"
out_dir         <- "output/step6_variancePartition"

celltype        <- "Microglia"            # value in the celltype column
genotype_levels <- c("WT", "HET", "KO")
min_cpm         <- 1      # keep genes with CPM > min_cpm ...
min_sample_frac <- 0.5    # ... in at least this fraction of samples

# ---- Setup -------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
seu <- readRDS(input_rds)

# ---- Pseudobulk: sum raw RNA counts per sample --------------------------------
sub <- subset(seu, cells = colnames(seu)[seu$celltype == celltype])
pb  <- AggregateExpression(sub, assays = "RNA", group.by = "sample")$RNA

# Sample metadata (one row per sample). AggregateExpression turns "_" into "-" in names.
# The genotype column is named Plcg2 so the variance components are labeled Plcg2 and Plcg2:sex.
info <- unique(sub[[]][, c("sample", "genotype", "sex", "run")])
rownames(info) <- gsub("_", "-", info$sample)
info <- info[colnames(pb), ]
info$Plcg2 <- factor(info$genotype, levels = genotype_levels)
info$sex   <- factor(info$sex)

# ---- voom, remove run effect, variance partition ------------------------------
keep <- rowSums(cpm(pb) > min_cpm) >= min_sample_frac * ncol(pb)
dge  <- calcNormFactors(DGEList(counts = pb[keep, ]))
vobj <- voom(dge, model.matrix(~ Plcg2 + sex + Plcg2:sex + run, data = info))

expr <- removeBatchEffect(vobj$E, batch = info$run, design = model.matrix(~ Plcg2 * sex, data = info))
vp   <- sortCols(fitExtractVarPartModel(expr, ~ (1 | Plcg2) + (1 | sex) + (1 | Plcg2:sex), info))

write.csv(vp, file.path(out_dir, paste0("variancePartition_", celltype, ".csv")))
