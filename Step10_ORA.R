# Step10_ORA.R
# Per cell type, using the MAST DEGs from step 8:
#   1. Shared DEGs between two comparisons: direction, correlation/regression of log2FC, monotonic genes
#   2. KEGG and MSigDB Hallmark over-representation analysis (clusterProfiler)
# Then KO vs WT log2FC of KEGG mmu05022 genes across microglial clusters (Fig 5E).
# Needs internet access (KEGG).
# Input:  output/step3_clustering/Plcg2_annotated.rds, output/step8_MAST/*.csv,
#         config/mh.all.v2026.1.Mm.entrez.gmt
# Output: output/step10_ORA/<celltype>/{gene_sets,KEGG,MSigDB}/, output/step10_ORA/shared_DEG_statistics.csv,
#         output/step10_ORA/Microglia/neurodegeneration_mmu05022_cluster_log2FC.csv
# Run:    Rscript Step10_ORA.R

suppressPackageStartupMessages({
  library(Seurat)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(KEGGREST)
})

# ---- Config ------------------------------------------------------------------
input_rds     <- "output/step3_clustering/Plcg2_annotated.rds"
deg_dir       <- "output/step8_MAST"
gmt_file      <- "config/mh.all.v2026.1.Mm.entrez.gmt"   # MSigDB Hallmark, mouse, Entrez IDs
out_dir       <- "output/step10_ORA"

celltypes     <- c("Astrocytes", "Microglia")   # values in the celltype column
sexes         <- c("F", "M")
comparisons   <- c(KO = "KO_vs_WT", HET = "HET_vs_WT", KOvHET = "KO_vs_HET")
padj_cutoff   <- 0.05

# Pairs of DEG lists to compare as c(x, y); log2FC of y is regressed on x.
#   Female vs male pairs: their same- and opposite-direction genes also go to ORA.
#   KO vs HET pairs (same sex): monotonic genes are also written.
shared_pairs  <- list(c("F_KO", "M_KO"), c("F_HET", "M_HET"), c("F_KO", "F_HET"), c("M_KO", "M_HET"))
fig5e_pathway <- "mmu05022"   # KEGG Pathways of neurodegeneration - multiple diseases

# ---- Setup -------------------------------------------------------------------
seu  <- readRDS(input_rds)
m_df <- read.gmt(gmt_file)
to_entrez <- function(x) bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
stats <- list()
sig_by_celltype <- list()

for (ct in celltypes) {
  message("==== ", ct, " ====")
  ct_dir <- file.path(out_dir, ct)
  for (d in c("gene_sets", "KEGG", "MSigDB")) dir.create(file.path(ct_dir, d), recursive = TRUE, showWarnings = FALSE)
  
  # ---- Background (universe): genes detected in this cell type ---------------
  counts   <- GetAssayData(seu, assay = "RNA", layer = "counts")[, seu$celltype == ct]
  universe <- to_entrez(rownames(counts)[rowSums(counts > 0) > 0])
  
  # ---- Significant DEGs, named like F_KO, M_HET, F_KOvHET ---------------------
  deg <- list()
  for (s in sexes) {
    for (cmp in names(comparisons)) {
      file <- file.path(deg_dir, sprintf("%s_MAST_%s_%s.csv", ct, comparisons[[cmp]], s))
      if (!file.exists(file)) { warning("Missing (skipped in step 8?): ", file); next }
      res <- read.csv(file, row.names = 1)
      deg[[paste0(s, "_", cmp)]] <- res[res$p_val_adj < padj_cutoff, ]
    }
  }
  sig  <- lapply(deg, rownames)
  up   <- lapply(deg, function(d) rownames(d)[d$avg_log2FC > 0])
  down <- lapply(deg, function(d) rownames(d)[d$avg_log2FC < 0])
  saveRDS(list(sig = sig, up = up, down = down), file.path(ct_dir, "gene_sets.rds"))
  sig_by_celltype[[ct]] <- sig
  
  # ---- 1. Shared DEGs between two comparisons ---------------------------------
  ora_sets <- sig   # every DEG list goes to ORA; same/opposite-direction sets are added below
  for (pr in shared_pairs) {
    x <- pr[1]; y <- pr[2]
    shared <- intersect(sig[[x]], sig[[y]])
    lfc_x  <- deg[[x]][shared, "avg_log2FC"]
    lfc_y  <- deg[[y]][shared, "avg_log2FC"]
    same   <- sign(lfc_x) == sign(lfc_y)
    
    # Each shared gene with both log2FCs and its direction (up_both, down_both, up_<x>_down_<y>, ...)
    direction <- ifelse(same, ifelse(lfc_x > 0, "up_both", "down_both"),
                        ifelse(lfc_x > 0, paste0("up_", x, "_down_", y), paste0("down_", x, "_up_", y)))
    tab <- setNames(data.frame(shared, lfc_x, lfc_y, direction), c("gene", x, y, "direction"))
    write.csv(tab, file.path(ct_dir, "gene_sets", sprintf("%s_shared_DEGs_%s_%s.csv", ct, x, y)), row.names = FALSE)
    
    # Correlation and linear fit of log2FC (y ~ x)
    r   <- cor.test(lfc_x, lfc_y, method = "pearson")
    fit <- lm(lfc_y ~ lfc_x)
    ci  <- confint(fit)["lfc_x", ]
    stats[[length(stats) + 1]] <- data.frame(
      celltype = ct, x = x, y = y, n_shared = length(shared),
      n_concordant = sum(same), n_discordant = sum(!same),
      pearson_r = unname(r$estimate), p_value = r$p.value, R2 = summary(fit)$r.squared,
      slope = coef(fit)[["lfc_x"]], slope_CI_low = ci[[1]], slope_CI_high = ci[[2]])
    
    if (substr(x, 1, 1) != substr(y, 1, 1)) {
      # Female vs male: genes changing in the same / opposite direction in both sexes go to ORA
      g <- sub("^F_", "", x)
      ora_sets[[paste0("samebysex_", g)]]     <- shared[same]
      ora_sets[[paste0("oppositebysex_", g)]] <- shared[!same]
    } else {
      # KO vs HET: monotonic genes change in the same direction, more in KO
      write.csv(tab[same & abs(lfc_x) > abs(lfc_y), ],
                file.path(ct_dir, "gene_sets", sprintf("%s_monotonic_%s_%s.csv", ct, x, y)), row.names = FALSE)
    }
  }
  
  # ---- 2. Over-representation analysis: KEGG (all sets), MSigDB (DEG lists only) ----
  for (nm in names(ora_sets)) {
    if (length(ora_sets[[nm]]) == 0) next
    message("ORA: ", ct, " ", nm)
    genes <- to_entrez(ora_sets[[nm]])
    
    kegg <- enrichKEGG(genes, organism = "mmu", keyType = "ncbi-geneid", universe = universe,
                       pvalueCutoff = 1, qvalueCutoff = 1)
    write.csv(as.data.frame(kegg), file.path(ct_dir, "KEGG", sprintf("%s_KEGG_%s.csv", ct, nm)), row.names = FALSE)
    
    if (!nm %in% names(sig)) next
    msig <- enricher(genes, TERM2GENE = m_df, universe = universe, pvalueCutoff = 1, qvalueCutoff = 1)
    write.csv(as.data.frame(msig), file.path(ct_dir, "MSigDB", sprintf("%s_MSigDB_%s.csv", ct, nm)), row.names = FALSE)
  }
}

write.csv(do.call(rbind, stats), file.path(out_dir, "shared_DEG_statistics.csv"), row.names = FALSE)

# ---- Fig 5E: KO vs WT log2FC of mmu05022 genes across microglial clusters -----
# Genes: pathway members that are microglia KO vs WT DEGs in either sex
kegg_ids      <- keggGet(fig5e_pathway)[[1]]$GENE[c(TRUE, FALSE)]   # odd elements are Entrez IDs
pathway_genes <- unname(mapIds(org.Mm.eg.db, keys = kegg_ids, column = "SYMBOL", keytype = "ENTREZID"))
mg_sig        <- sig_by_celltype[["Microglia"]]
genes         <- intersect(pathway_genes, union(mg_sig$F_KO, mg_sig$M_KO))

# One column per cluster and sex from step 8, named like F_5_M (genes MAST did not report are 0)
files  <- list.files(deg_dir, pattern = "^Microglia_cluster_.*_MAST_KO_vs_WT_[FM]\\.csv$")
labels <- sub("^Microglia_cluster_(.*)_MAST_KO_vs_WT_([FM])\\.csv$", "\\2_\\1", files)
lfc <- sapply(file.path(deg_dir, files), function(f) read.csv(f, row.names = 1)[genes, "avg_log2FC"])
dimnames(lfc) <- list(genes, labels)
lfc[is.na(lfc)] <- 0
lfc <- lfc[rowSums(lfc != 0) > 0, , drop = FALSE]

# Columns: females then males, clusters in numeric order
lfc <- lfc[, order(substr(labels, 1, 1), as.numeric(sub("^[FM]_(\\d+)_.*", "\\1", labels)))]
write.csv(cbind(gene = rownames(lfc), as.data.frame(lfc)),
          file.path(out_dir, "Microglia", "neurodegeneration_mmu05022_cluster_log2FC.csv"), row.names = FALSE)
