# Step5_ADRiskGenes.R
# AD GWAS risk genes: select genes enriched in microglia (mean log-normalized expression in microglia
# above every other cell type), then summarize their expression across microglial clusters
# (Table S5) and as cluster averages (for Fig 1E). Plcg2 is always included.
# Input:  output/step3_clustering/Plcg2_annotated.rds
# Output: output/step5_ADRiskGenes/
# Run:    Rscript Step5_ADRiskGenes.R

suppressPackageStartupMessages(library(Seurat))

# ---- Config ------------------------------------------------------------------
input_rds       <- "output/step3_clustering/Plcg2_annotated.rds"
out_dir         <- "output/step5_ADRiskGenes"
microglia_label <- "Microglia"   # value in the celltype (gene selection) and subtype (clusters) columns

ad_genes <- c(  # AD GWAS risk genes (mouse orthologs)
  "Sort1", "Cr1l", "Adam17", "Prkd3", "Nck2", "Bin1", "Wdr12", "Inpp5d", "Mme", "Idua",
  "Clnk", "Rhoh", "Ank", "Cox7c", "Tnip1", "Rasgef1c", "H2-Aa", "Unc5cl", "Trem2",
  "Treml2", "Cd2ap", "Hs3st5", "Umad1", "Ica1", "Tmem106b", "Jazf1", "Epdr1",
  "Sec61g", "Epha1", "Ctsb", "Ptk2b", "Clu", "Sharpin", "Abca1", "Usp6nl",
  "Ank3", "Tspan14", "Blnk", "Plekha1", "Spi1", "Ms4a4a", "Eed", "Sorl1", "Tpcn1",
  "Fermt2", "Slc24a4", "Sppl2a", "Mindy2", "Aph1b", "Snx1", "Ctsh", "Doc2a", "Bckdk",
  "Il34", "Maf", "Plcg2", "Foxf1", "Prdm9", "Wdr81", "Scimp", "Myo15a", "Grn",
  "Wnt3", "Abi3", "Tspoap1", "Ace", "Abca7", "Klf16", "Siglecg", "Pirb",
  "Rbck1", "Cass4", "App", "Adamts1")

# ---- Setup -------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
seu <- readRDS(input_rds)
DefaultAssay(seu) <- "RNA"

# ---- 1. Microglia-enriched genes (mean RNA log-normalized expression per cell type) ----
expr <- FetchData(seu, vars = ad_genes, layer = "data")
avg  <- sapply(split(expr, seu$celltype), colMeans)   # genes x cell types

selection <- data.frame(
  Gene          = rownames(avg),
  Microglia_Avg = avg[, microglia_label],
  Max_Other_Avg = apply(avg[, colnames(avg) != microglia_label, drop = FALSE], 1, max)
)
selection$Difference <- selection$Microglia_Avg - selection$Max_Other_Avg
selection$Selected   <- selection$Difference > 0
selection <- selection[order(selection$Gene), ]

genes <- unique(c(selection$Gene[selection$Selected], "Plcg2"))

# ---- 2. Expression per microglial cluster (Table S5) --------------------------
microglia <- subset(seu, cells = colnames(seu)[seu$subtype == microglia_label])
mexpr <- FetchData(microglia, vars = genes, layer = "data")

summary_tab <- do.call(rbind, lapply(genes, function(g) {
  x <- split(mexpr[[g]], microglia$cluster_label)
  data.frame(cluster_label       = names(x),
             Gene                = g,
             Number_of_cells     = lengths(x),
             Number_expressing   = sapply(x, function(v) sum(v > 0)),
             Percent_expressing  = 100 * sapply(x, function(v) mean(v > 0)),
             Mean_log_expression = sapply(x, mean))
}))
summary_tab <- summary_tab[order(summary_tab$cluster_label, summary_tab$Gene), ]
write.csv(summary_tab, file.path(out_dir, "AD_genes_Plcg2_microglia_cluster_expression.csv"), row.names = FALSE)

# ---- 3. Average expression per microglial cluster (Fig 1E) -------------------
avg_mg <- as.matrix(AverageExpression(microglia, assays = "RNA", layer = "data",
                                      group.by = "cluster_label", features = genes)$RNA)
avg_mg <- avg_mg[order(rownames(avg_mg)), order(colnames(avg_mg)), drop = FALSE]
write.csv(avg_mg, file.path(out_dir, "AD_genes_microglia_cluster_average.csv"))
