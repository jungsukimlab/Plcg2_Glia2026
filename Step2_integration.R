# Step2_integration.R
# Merge per-sample objects, remove doublets, SCTransform, and integrate with CCA (Seurat v5).
# Input:  config/samples.csv, output/step1_preprocessing/<ID>_postSoupX_DF.rds
# Output: output/step2_integration/Plcg2_integrated.rds
# Run:    Rscript Step2_integration.R

suppressPackageStartupMessages(library(Seurat))

# ---- Config ------------------------------------------------------------------
sample_sheet <- "config/samples.csv"
in_dir       <- "output/step1_preprocessing"
out_dir      <- "output/step2_integration"

# ---- Setup -------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load, merge, filter -----------------------------------------------------
ids  <- read.csv(sample_sheet, stringsAsFactors = FALSE)$ID
objs <- lapply(ids, function(id) readRDS(file.path(in_dir, paste0(id, "_postSoupX_DF.rds"))))

# Prefix barcodes with the sample ID so identical barcodes across samples can't collide
seu <- merge(objs[[1]], y = objs[-1], add.cell.ids = ids)
seu <- subset(seu, subset = DF == "Singlet")   # percent.mt was already filtered in step 1

# ---- Normalize and integrate -------------------------------------------------
# Each sample stays a separate layer, so SCTransform runs per sample
seu <- SCTransform(seu)
seu <- RunPCA(seu)
seu <- IntegrateLayers(seu, method = CCAIntegration, normalization.method = "SCT", verbose = FALSE)

saveRDS(seu, file.path(out_dir, "Plcg2_integrated.rds"))

