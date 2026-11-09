# Step8_MASTDEGs.R
# Single-cell DE (MAST) between genotypes, run separately in females and males, per cell type and
# (KO vs WT only) per microglial cluster.
# Input:  output/step3_clustering/Plcg2_annotated.rds
# Output: output/step8_MAST/<celltype>_MAST_<g1>_vs_<g2>_<sex>.csv,
#         output/step8_MAST/Microglia_cluster_<cluster>_MAST_KO_vs_WT_<sex>.csv
# Run:    Rscript Step8_MASTDEGs.R

suppressPackageStartupMessages({
  library(Seurat)
  library(MAST)
})

# ---- Config ------------------------------------------------------------------
input_rds       <- "output/step3_clustering/Plcg2_annotated.rds"
out_dir         <- "output/step8_MAST"

sexes           <- c("F", "M")
contrasts       <- list(c("KO", "WT"), c("HET", "WT"), c("KO", "HET"))   # ident.1 vs ident.2
microglia_label <- "Microglia"   # value in the celltype column
latent_vars     <- c("percent.mt", "percent.rb", "nCount_SCT_data")      # from steps 1 and 3
min_cells       <- 3             # skip a comparison if either genotype has fewer cells
seed            <- 42


# ---- Setup -------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seu <- readRDS(input_rds)
Idents(seu) <- "genotype"

# ---- DEG ---------------------------------------------------------------------
for (sex in sexes) {
  seu_sex <- subset(seu, cells = colnames(seu)[seu$sex == sex])
  
  # 1. Each cell type, all contrasts
  for (ct in sort(unique(seu_sex$celltype))) {
    x <- subset(seu_sex, cells = colnames(seu_sex)[seu_sex$celltype == ct])
    
    for (cn in contrasts) {
      file <- file.path(out_dir, sprintf("%s_MAST_%s_vs_%s_%s.csv",
                                         gsub("[^A-Za-z0-9]+", "_", ct), cn[1], cn[2], sex))
      n <- table(factor(Idents(x), levels = cn))
      if (any(n < min_cells)) { message("Skipping ", basename(file)); next }
      
      message("Running ", basename(file))
      write.csv(FindMarkers(x, assay = "SCT", ident.1 = cn[1], ident.2 = cn[2], recorrect_umi = FALSE,
                            latent.vars = latent_vars, test.use = "MAST"), file)
    }
  }
  
  # 2. Each microglial cluster, KO vs WT only
  microglia <- subset(seu_sex, cells = colnames(seu_sex)[seu_sex$celltype == microglia_label])
  
  for (cl in sort(unique(microglia$cluster_label))) {
    x <- subset(microglia, cells = colnames(microglia)[microglia$cluster_label == cl])
    
    file <- file.path(out_dir, sprintf("Microglia_cluster_%s_MAST_KO_vs_WT_%s.csv", cl, sex))
    n <- table(factor(Idents(x), levels = c("KO", "WT")))
    if (any(n < min_cells)) { message("Skipping ", basename(file)); next }
    
    message("Running ", basename(file))
    write.csv(FindMarkers(x, assay = "SCT", ident.1 = "KO", ident.2 = "WT", recorrect_umi = FALSE,
                          latent.vars = latent_vars, test.use = "MAST"), file)
  }
}
