# Step3_clustering.R
# Cluster the integrated object, find cluster markers, annotate clusters, and test cell-type and
# cluster proportions between genotypes (speckle propeller).
# Input:  output/step2_integration/Plcg2_integrated.rds
#         config/Annotated_clusters_0.6.xlsx, one row per cluster:
#         column 1 = cluster number, 2 = cell type, 3 = subtype, 4 = cluster label
# Output: output/step3_clustering/Plcg2_annotated.rds, cluster marker tables, cell counts and
#         proportion tests
# Run:    Rscript Step3_clustering.R

suppressPackageStartupMessages({
  library(Seurat)
  library(readxl)
  library(speckle)
})

# ---- Config ------------------------------------------------------------------
input_rds       <- "output/step2_integration/Plcg2_integrated.rds"
annotation_file <- "config/Annotated_clusters_0.6.xlsx"
out_dir         <- "output/step3_clustering"

n_dims          <- 50        # PCs from integrated.dr for neighbors and UMAP
cluster_res     <- 0.6      
rb_pattern      <- "^Rp"     # ribosomal genes 
microglia_label <- "Microglia"

# ---- Setup -------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

# ---- Clustering --------------------------------------------------------------
seu <- readRDS(input_rds)

seu <- FindNeighbors(seu, reduction = "integrated.dr", dims = 1:n_dims)
seu <- RunUMAP(seu, reduction = "integrated.dr", dims = 1:n_dims, reduction.name = "umap.dr")
seu <- FindClusters(seu, resolution = cluster_res)

seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
seu[["RNA"]] <- NormalizeData(seu[["RNA"]])
seu[["percent.rb"]] <- PercentageFeatureSet(seu, pattern = rb_pattern)

# Prepare the SCT assay for marker/DE testing
seu <- PrepSCTFindMarkers(seu)
DefaultAssay(seu) <- "SCT"
seu$nCount_SCT_data <- colSums(GetAssayData(seu, assay = "SCT", layer = "data"))

# ---- Cluster markers ---------------------------------------------------------
Idents(seu) <- "seurat_clusters"
write.csv(FindAllMarkers(seu, assay = "SCT", only.pos = TRUE),
          file.path(out_dir, paste0("FindAllMarkers_SCT_wilcox_res", cluster_res, ".csv")))

# ---- Annotate clusters -------------------------------------------------------
annot <- as.data.frame(read_excel(annotation_file))

label_by_cluster <- function(col_idx) {
  map <- setNames(as.character(annot[[col_idx]]), as.character(annot[[1]]))
  unname(map[as.character(seu$seurat_clusters)])
}
seu$celltype      <- label_by_cluster(2)
seu$subtype       <- label_by_cluster(3)
seu$cluster_label <- label_by_cluster(4)

# Microglia-only markers between cluster labels (SCT re-prepared on the subset)
microglia <- subset(seu, cells = colnames(seu)[seu$subtype == microglia_label])
microglia <- PrepSCTFindMarkers(microglia)
Idents(microglia) <- "cluster_label"
write.csv(FindAllMarkers(microglia, assay = "SCT", only.pos = TRUE),
          file.path(out_dir, paste0("FindAllMarkers_Microglia_SCT_wilcox_res", cluster_res, ".csv")))

# ---- Cell counts and proportions by genotype (speckle propeller) -------------
write.csv(table(cluster = seu$cluster_label, genotype = seu$genotype, sample = seu$sample),
          file.path(out_dir, "cells_per_cluster.csv"))
write.csv(table(subtype = seu$subtype, genotype = seu$genotype, sample = seu$sample),
          file.path(out_dir, "cells_per_subtype.csv"))

write.csv(propeller(clusters = seu$seurat_clusters, sample = seu$sample, group = seu$genotype),
          file.path(out_dir, "propeller_clusters.csv"))
write.csv(propeller(clusters = seu$subtype, sample = seu$sample, group = seu$genotype),
          file.path(out_dir, "propeller_subtype.csv"))

# Microglial clusters within each sex
for (s in c("F", "M")) {
  x <- subset(microglia, cells = colnames(microglia)[microglia$sex == s])
  write.csv(propeller(clusters = x$seurat_clusters, sample = x$sample, group = x$genotype),
            file.path(out_dir, paste0("propeller_Microglia_clusters_", s, ".csv")))
}

saveRDS(seu, file.path(out_dir, "Plcg2_annotated.rds"))
