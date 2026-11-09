# Step4_UCell.R
# UCell scores for homeostatic and inflammatory microglia signatures across microglial clusters
# (Fig 1D, Fig S2A-B). Inflammatory = union of the IRM, CRM, DAM and HLA signatures.
# Input:  output/step3_clustering/Plcg2_annotated.rds
# Output: output/step4_UCell/UCell_scores_Microglia.csv (one row per microglial cell)
# Run:    Rscript Step4_UCell.R

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
})

# ---- Config ------------------------------------------------------------------
input_rds       <- "output/step3_clustering/Plcg2_annotated.rds"
out_dir         <- "output/step4_UCell"
microglia_label <- "Microglia"   # value in the subtype column

# Microglia state signatures (gene lists as used in the paper)
IRM <- c(  # interferon response microglia
  "Adar", "Apoe", "Axl", "B2m", "Bcl2a1d", "Bst2", "Ccl9", "Ccl2", "Ccl3", "Ccl4",
  "Cd83", "Clec7a", "Cmpk2", "Csf1", "Cst7", "Ctsl", "Cxcl10", "Gnas", "Ifi213", "Ifi35",
  "Ifi44", "Ifit2", "Ifit3", "Ifitm3", "Igf1", "Irf7", "Isg15", "Itgax", "Lgals3bp",
  "Ly6e", "Mx1", "Oas1a", "Oas2", "Oas3", "Oasl1", "Parp14", "Parp9", "Rnf213", "Rsad2",
  "Samd9l", "Serpine2", "Spp1", "Stat1", "Stat2", "Trim5", "Usp18", "Xaf1")
CRM <- c(  # cytokine response microglia
  "Ank", "Apoe", "Axl", "B2m", "Baiap2l2", "Cadm1", "Capg", "Ccl9", "Ccl3", "Ccl4",
  "Cd52", "Cd68", "Cd83", "Cebpb", "Clec7a", "Csf1", "Cst7", "Ctsa", "Ctsz",
  "Fth1", "Gusb", "H2-D1", "H2-K1", "Hexa", "Hspa1a", "Igf1", "Itgax", "Lgals3bp", "Lpl",
  "Lyz1", "Npc2", "Or5v1", "Pld3", "Plek", "Rgs1", "Rpl18a", "Rpl32", "Rplp1", "Rps14",
  "Rps28", "Serpine2", "Spp1", "Tyrobp")
DAM <- c(  # disease-associated microglia
  "Cd9", "Spp1", "Cd83", "Ramp1", "Ctsd", "Ccl3", "Trem2", "Cpnmb", "Itgax", "Lpl",
  "Cd63", "Slc11a1", "Fabp3", "Myo1e", "Apoe", "Npc2", "Ctsb", "Ctsz", "Axl", "Rpl37a",
  "Rpl19", "Rplp1", "Cadm1", "Ank", "Cst7", "Clec7a", "H2-D1", "Fth1", "Lyz1",
  "Serpine2", "Ctsl", "Lgals3bp", "Cd52", "Ccl9", "H2-K1", "Cd68", "Csf1", "Cd63", "Or5v1")
HLA <- c(  # HLA (antigen-presenting) microglia
  "Cd74", "H2-Ab1", "Apoe", "Cd9", "B2m", "Rps19", "Clec7a", "Rps28", "Rps14", "Rps23",
  "Rpl32", "Rpl41", "Rps16", "Rpl13", "Rpl19", "Lyz1", "Rplp1", "Rpl13a", "Npc2", "Rplp2",
  "Itgax", "Rpl37a", "Rpl23", "Ftl1", "Tyrobp", "Trem2", "Lpl", "Cadm1", "Ctsz", "Cts7",
  "H2-D1", "Fth1", "Serpine2", "Ank", "Ctsl", "Lgals3bp", "Cd52", "Ccl9", "H2-K1",
  "Cd68", "Csf1", "Cd63", "Or5v1", "Hexa", "Ctsa", "Gusb", "Baiap2l2", "Ifg1")
signatures <- list(
  HM = c(  # homeostatic microglia
    "Adgrg1", "Adrb2", "Alox5ap", "Bin1", "Bin1", "Ccnd1", "Ccr5", "Chd9", "Clec4b1",
    "Cmtm6", "Csf1r", "Cx3cr1", "Cysltr1", "Epb41l2", "Otulinl", "Frmd4a", "Fscn1",
    "Gal3st4", "Gpr34", "Adgrg1", "Ifngr1", "Lpar6", "Lrba", "Ltc4s", "Maf", "Malat1",
    "Marcks", "Olfml3", "P2ry12", "P2ry13", "Pmepa1", "Ptgs1", "Rgs2", "Ear3",
    "Selplg", "Slc2a5", "Slc7a8", "Slco2b1", "Sparc", "Srgap2", "Ssh2", "Stab1",
    "Tmem119", "Sting1", "Txnip", "Zfhx3"),
  Inflammatory = sort(unique(c(IRM, CRM, DAM, HLA)))
)

# ---- Setup -------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
seu <- readRDS(input_rds)

# ---- UCell scores in microglia -----------------------------------------------
microglia <- subset(seu, cells = colnames(seu)[seu$subtype == microglia_label])
microglia <- AddModuleScore_UCell(microglia, features = signatures, assay = "RNA")
score_cols <- paste0(names(signatures), "_UCell")

write.csv(microglia[[]][, c("sample", "genotype", "sex", "cluster_label", score_cols)],
          file.path(out_dir, "UCell_scores_Microglia.csv"))
