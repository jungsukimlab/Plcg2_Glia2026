# Step12_CellChat.R
# CellChat cell-cell communication: one CellChat object per genotype (WT, KO) x sex, then
# differential ligand-receptor signaling KO vs WT within each sex.
# Input:  output/step3_clustering/Plcg2_annotated.rds
# Output: output/step12_CellChat/{objects,tables}/
# Run:    Rscript Step12_CellChat.R

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(future)
})

# ---- Config ------------------------------------------------------------------
input_rds       <- "output/step3_clustering/Plcg2_annotated.rds"
out_dir         <- "output/step12_CellChat"

group_col       <- "celltype"                    # cell groups used for communication
sexes           <- c("F", "M")
ref_genotype    <- "WT"                          # one CellChat object per genotype x sex;
test_genotype   <- "KO"                          # compared test vs reference within each sex

min_cells       <- 10      # drop communication from groups with fewer cells
n_workers       <- 4       # parallel workers (1 = sequential)
options(future.globals.maxSize = 4e9)

# Differential expression between datasets (ligand-receptor level)
thresh_pc       <- 0.1
thresh_fc       <- 0.05
thresh_p        <- 0.05
ligand_logFC    <- 0.05

# ---- Setup -------------------------------------------------------------------
for (d in c("objects", "tables")) dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
out <- function(d, f) file.path(out_dir, d, f)

seu <- readRDS(input_rds)
seu$samples <- seu$sample   # CellChat expects a "samples" column
md <- seu[[]]

# ---- 1. One CellChat object per genotype x sex -------------------------------
if (n_workers > 1) plan("multisession", workers = n_workers)
cc <- list()

for (g in c(ref_genotype, test_genotype)) {
  for (s in sexes) {
    nm <- paste(g, s, sep = "_")
    message("CellChat: ", nm)

    sub <- subset(seu, cells = rownames(md)[md$genotype == g & md$sex == s])
    Idents(sub) <- group_col

    x <- createCellChat(sub, group.by = "ident", assay = "RNA")
    x@DB <- CellChatDB.mouse
    x <- subsetData(x)
    x <- identifyOverExpressedGenes(x)
    x <- identifyOverExpressedInteractions(x)
    x <- computeCommunProb(x)
    x <- filterCommunication(x, min.cells = min_cells)
    x <- computeCommunProbPathway(x)
    x <- aggregateNet(x)
    x <- netAnalysis_computeCentrality(x, slot.name = "netP")   # centrality of the pathway-level network

    saveRDS(x, out("objects", paste0("cellchat_", nm, ".rds")))
    cc[[nm]] <- x
  }
}
plan("sequential")

# ---- 2. Test vs reference genotype within each sex ---------------------------
for (s in sexes) {
  tag <- sprintf("%s_vs_%s_%s", test_genotype, ref_genotype, s)
  object_list <- setNames(list(cc[[paste(ref_genotype, s, sep = "_")]], cc[[paste(test_genotype, s, sep = "_")]]),
                          c(ref_genotype, test_genotype))
  merged <- mergeCellChat(object_list, add.names = names(object_list), cell.prefix = TRUE)

  # Differentially expressed ligands/receptors in test vs reference, mapped to the L-R network
  features_name <- paste0(test_genotype, ".merged")
  merged <- identifyOverExpressedGenes(merged, group.dataset = "datasets", pos.dataset = test_genotype,
                                       features.name = features_name, only.pos = FALSE,
                                       thresh.pc = thresh_pc, thresh.fc = thresh_fc, thresh.p = thresh_p,
                                       group.DE.combined = FALSE)
  net      <- netMappingDEG(merged, features.name = features_name, variable.all = TRUE)
  net_up   <- subsetCommunication(merged, net = net, datasets = test_genotype,
                                  ligand.logFC = ligand_logFC, receptor.logFC = NULL)
  net_down <- subsetCommunication(merged, net = net, datasets = ref_genotype,
                                  ligand.logFC = -ligand_logFC, receptor.logFC = NULL)

  write.csv(net, out("tables", paste0("Net_DEA_CellChat_", tag, ".csv")))
  saveRDS(list(object_list = object_list, merged = merged, net = net, net_up = net_up, net_down = net_down),
          out("objects", paste0("comparison_", tag, ".rds")))
}
