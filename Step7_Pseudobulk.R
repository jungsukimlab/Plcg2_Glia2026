# Step7_Pseudobulk.R
# Pseudobulk DESeq2 per cell type (one profile per sample): genotype effects in all samples,
# females only and males only, plus the genotype x sex interaction.
# Input:  output/step3_clustering/Plcg2_annotated.rds
# Output: output/step7_pseudobulk/<celltype>_*.csv
# Run:    Rscript Step7_Pseudobulk.R

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(edgeR)          # cpm() for the low-count filter
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

# ---- Config ------------------------------------------------------------------
input_rds       <- "output/step3_clustering/Plcg2_annotated.rds"
out_dir         <- "output/step7_pseudobulk"

celltypes       <- c(Astrocytes = "celltype", Microglia = "celltype")   # cell type = metadata column
genotype_levels <- c("WT", "HET", "KO")                                # first level is the reference
contrasts       <- list(c("KO", "WT"), c("HET", "WT"))
min_cpm         <- 5      # keep genes with CPM > min_cpm ...
min_sample_frac <- 0.75   # ... in at least this fraction of samples

# ---- Setup -------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
seu <- readRDS(input_rds)

# Protein-coding genes
gene_types <- AnnotationDbi::select(org.Mm.eg.db, keys = rownames(seu[["RNA"]]), keytype = "SYMBOL",
                                    columns = c("SYMBOL", "GENETYPE"))
pc_genes <- unique(na.omit(gene_types$SYMBOL[gene_types$GENETYPE == "protein-coding"]))

for (celltype in names(celltypes)) {
  prefix <- file.path(out_dir, celltype)

  # ---- Pseudobulk: sum raw RNA counts per sample ------------------------------
  sub <- subset(seu, cells = colnames(seu)[seu[[celltypes[[celltype]]]][, 1] == celltype])
  pb  <- round(AggregateExpression(sub, assays = "RNA", group.by = "sample")$RNA)
  pb  <- pb[rownames(pb) %in% pc_genes, ]
  pb  <- pb[rowSums(cpm(pb) > min_cpm) >= min_sample_frac * ncol(pb), ]

  # Sample metadata (one row per sample). AggregateExpression turns "_" into "-" in names.
  info <- unique(sub[[]][, c("sample", "genotype", "sex", "run")])
  rownames(info) <- gsub("_", "-", info$sample)
  info <- info[colnames(pb), ]
  info$genotype <- factor(info$genotype, levels = genotype_levels)
  info$sex      <- factor(info$sex, levels = c("F", "M"))
  info$run      <- factor(info$run)

  # ---- Genotype effects (Wald): all samples, females, males -------------------
  scopes <- list(AllSex = rep(TRUE, nrow(info)), F = info$sex == "F", M = info$sex == "M")
  for (scope in names(scopes)) {
    keep   <- scopes[[scope]]
    design <- if (scope == "AllSex") ~ sex + run + genotype else ~ run + genotype
    dds <- DESeq(DESeqDataSetFromMatrix(pb[, keep], droplevels(info[keep, ]), design), fitType = "local")
    for (cn in contrasts) {
      res <- results(dds, contrast = c("genotype", cn[1], cn[2]))
      write.csv(res[order(res$pvalue), ],
                sprintf("%s_GenotypeEffect_%s_vs_%s_%s.csv", prefix, cn[1], cn[2], scope))
    }
  }

  # ---- Genotype x sex interaction (Wald) --------------------------------------
  # Coefficient = genotype effect in males minus genotype effect in females
  dds <- DESeq(DESeqDataSetFromMatrix(pb, info, ~ genotype * sex + run), fitType = "local")
  for (g in setdiff(genotype_levels, genotype_levels[1])) {
    res <- results(dds, name = paste0("genotype", g, ".sexM"))
    write.csv(res[order(res$pvalue), ], sprintf("%s_InteractionEffect_%s.csv", prefix, g))
  }
}
