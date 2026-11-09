# Step9_GSEA.R
# KEGG and MSigDB Hallmark GSEA (clusterProfiler) on every pseudobulk DESeq2 result,
# with genes ranked by the DESeq2 Wald statistic. Needs internet access (KEGG).
# Input:  output/step7_pseudobulk/*.csv, config/mh.all.v2026.1.Mm.entrez.gmt
# Output: output/step9_GSEA/{KEGG,MSigDB}/<result>_*_GSEA.csv
# Run:    Rscript Step9_GSEA.R

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

# ---- Config ------------------------------------------------------------------
in_dir   <- "output/step7_pseudobulk"
gmt_file <- "config/mh.all.v2026.1.Mm.entrez.gmt"   # MSigDB Hallmark, mouse, Entrez IDs
out_dir  <- "output/step9_GSEA"

min_gs_size <- 10
max_gs_size <- 500

# ---- Setup -------------------------------------------------------------------
for (d in c("KEGG", "MSigDB")) dir.create(file.path(out_dir, d), recursive = TRUE, showWarnings = FALSE)
m_df <- read.gmt(gmt_file)

# ---- GSEA per DESeq2 result --------------------------------------------------
for (file in list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)) {
  name <- tools::file_path_sans_ext(basename(file))
  message("Running ", name)

  res <- read.csv(file, row.names = 1)   # gene symbols are the row names
  res <- res[is.finite(res$stat), ]

  # Symbol -> Entrez ID; if several symbols map to one ID, keep the largest |stat|
  ids  <- bitr(rownames(res), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  rank <- merge(data.frame(SYMBOL = rownames(res), stat = res$stat), ids, by = "SYMBOL")
  rank <- rank[order(-abs(rank$stat)), ]
  rank <- rank[!duplicated(rank$ENTREZID), ]
  gene_list <- sort(setNames(rank$stat, rank$ENTREZID), decreasing = TRUE)

  kegg <- gseKEGG(gene_list, organism = "mmu", keyType = "ncbi-geneid",
                  minGSSize = min_gs_size, maxGSSize = max_gs_size, pvalueCutoff = 1, eps = 0,
                  by = "fgsea", verbose = FALSE)
  write.csv(as.data.frame(kegg), file.path(out_dir, "KEGG", paste0(name, "_KEGG_GSEA.csv")), row.names = FALSE)

  msig <- GSEA(gene_list, TERM2GENE = m_df,
               minGSSize = min_gs_size, maxGSSize = max_gs_size, pvalueCutoff = 1, eps = 0,
               by = "fgsea", verbose = FALSE)
  write.csv(as.data.frame(msig), file.path(out_dir, "MSigDB", paste0(name, "_MSigDB_GSEA.csv")), row.names = FALSE)
}
