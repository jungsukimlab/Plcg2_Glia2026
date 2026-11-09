# Step1_preprocessing.R
# Per-sample QC: SoupX ambient RNA removal -> Seurat QC + SCTransform -> DoubletFinder.
# Input:  config/samples.csv, one row per sample. cellranger_dir is a Cell Ranger "outs" folder
#         with raw_feature_bc_matrix/, filtered_feature_bc_matrix/ and analysis/ (SoupX needs the clusters).
# Output: output/step1_preprocessing/<ID>_postSoupX_DF.rds
# Run:    Rscript Step1_preprocessing.R

suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
  library(DoubletFinder)
})

# ---- Config ------------------------------------------------------------------
sample_sheet       <- "config/samples.csv"   # columns: ID, sample, genotype, sex, run, cellranger_dir
out_dir            <- "output/step1_preprocessing"

min_cells          <- 3
min_features       <- 200
max_percent_mt     <- 10
mt_pattern         <- "^mt-"    
n_dims             <- 50        # PCs for neighbors
cluster_res        <- 0.2       # low on purpose: homotypic proportion should reflect broad cell types
n_pcs_doublet      <- 10
pN                 <- 0.25
multiplet_per_cell <- 4.6e-6    # multiplet rate per loaded cell (10x Genomics multiplet-rate table)
recovery_rate      <- 0.57      # recovered / loaded cells (same table)

# ---- Setup -------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
samples <- read.csv(sample_sheet, stringsAsFactors = FALSE)

# ---- Per sample --------------------------------------------------------------
for (i in seq_len(nrow(samples))) {
  info <- samples[i, ]
  out  <- function(suffix) file.path(out_dir, paste0(info$ID, "_", suffix))
  message("Processing ", info$ID, " (", i, "/", nrow(samples), ")")

  # 1. SoupX ambient RNA correction
  soup <- load10X(info$cellranger_dir)
  soup <- autoEstCont(soup, forceAccept = TRUE, doPlot = FALSE)
  seu <- CreateSeuratObject(adjustCounts(soup), project = paste0(info$sample, "_", info$ID),
                            min.cells = min_cells, min.features = min_features)

  # 2. QC and metadata
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = mt_pattern)
  for (col in c("sample", "ID", "genotype", "sex", "run")) seu[[col]] <- info[[col]]
  seu <- subset(seu, subset = nFeature_RNA > min_features & percent.mt < max_percent_mt)

  # 3. Normalization and clustering
  seu <- SCTransform(seu)
  seu <- RunPCA(seu)
  seu <- FindNeighbors(seu, dims = 1:n_dims)
  seu <- FindClusters(seu, resolution = cluster_res)

  # 4. DoubletFinder
  sweep <- summarizeSweep(paramSweep(seu, PCs = 1:n_pcs_doublet, sct = TRUE), GT = FALSE)
  bcmvn <- find.pK(sweep)
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  n_cells  <- ncol(seu)
  nExp     <- round(n_cells / recovery_rate * multiplet_per_cell * n_cells)
  nExp_adj <- round(nExp * (1 - modelHomotypic(seu$seurat_clusters)))
  seu <- doubletFinder(seu, PCs = 1:n_pcs_doublet, pN = pN, pK = pK,
                       nExp = nExp_adj, sct = TRUE)

  # Rename DoubletFinder's parameter-named columns to pANN / DF
  md_cols <- colnames(seu@meta.data)
  colnames(seu@meta.data)[grep("^pANN_", md_cols)] <- "pANN"
  colnames(seu@meta.data)[grep("^DF\\.classifications", md_cols)] <- "DF"

  saveRDS(seu, out("postSoupX_DF.rds"))
}

