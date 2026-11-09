# Step11_hdWGCNA.R
# hdWGCNA co-expression network for microglia: metacells, modules, eigengenes, hub genes,
# differential module eigengenes (DMEs) by genotype within each sex, GO enrichment (Enrichr)
# of each module's top hub genes, and module-trait correlations per cell type (Fig S12).
# Needs internet access (Enrichr).
# Input:  output/step3_clustering/Plcg2_annotated.rds
# Output: output/step11_hdWGCNA/
# Run:    Rscript Step11_hdWGCNA.R

suppressPackageStartupMessages({
  library(Seurat)
  library(WGCNA)
  library(hdWGCNA)
  library(UCell)
  library(enrichR)
})

# ---- Config ------------------------------------------------------------------
input_rds          <- "output/step3_clustering/Plcg2_annotated.rds"
out_dir            <- "output/step11_hdWGCNA"

wgcna_name         <- "Microglia"   # name of the hdWGCNA experiment
group_name         <- "Microglia"   # cell group the network is built on ...
group_col          <- "subtype"     # ... found in this metadata column

gene_fraction      <- 0.05          # gene must be expressed in this fraction of cells
metacell_group_by  <- c("sample", "subtype", "cluster_label")
k                  <- 25            # cells per metacell
max_share          <- 10            # max cells shared between metacells
network_type       <- "signed"
tom_name           <- "Micro"       # topological overlap matrix written to disk (TOM/)
n_score_genes      <- 25            # top hub genes per module for UCell scores
n_hubs             <- 100

sexes              <- c("F", "M")
contrasts          <- list(c("KO", "WT"))   # DMEs are tested in all cells

enrichr_dbs        <- c("GO_Biological_Process_2023", "GO_Molecular_Function_2023")
enrichr_n_genes    <- 50            # top hub genes per module (by kME) used for enrichment

# ---- Setup -------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Network -----------------------------------------------------------------
seu <- readRDS(input_rds)
DefaultAssay(seu) <- "RNA"

seu <- SetupForWGCNA(seu, gene_select = "fraction", fraction = gene_fraction, wgcna_name = wgcna_name)
seu <- MetacellsByGroups(seu, group.by = metacell_group_by, ident.group = group_col,
                         k = k, max_share = max_share)
seu <- NormalizeMetacells(seu)

seu <- SetDatExpr(seu, group_name = group_name, group.by = group_col, assay = "RNA")
seu <- TestSoftPowers(seu, networkType = network_type)
seu <- ConstructNetwork(seu, tom_name = tom_name, overwrite_tom = TRUE)

# ---- Module eigengenes, connectivity, scores ---------------------------------
seu <- ModuleEigengenes(seu)
seu <- AddMetaData(seu, GetMEs(seu, harmonized = FALSE))   # MEs as cell metadata (one column per module)
seu <- ModuleConnectivity(seu, group.by = group_col, group_name = group_name)
seu <- ModuleExprScore(seu, n_genes = n_score_genes, method = "UCell")

write.csv(GetHubGenes(seu, n_hubs = n_hubs), file.path(out_dir, paste0("hubgenes_top", n_hubs, ".csv")))
modules <- GetModules(seu)
write.csv(modules[modules$module != "grey", ], file.path(out_dir, "modules.csv"))

# ---- Differential module eigengenes (genotype within each sex) ---------------
md <- seu[[]]

dmes <- list()
for (s in sexes) {
  for (cn in contrasts) {
    nm <- sprintf("%s_vs_%s_%s", cn[1], cn[2], s)
    dmes[[nm]] <- FindDMEs(seu, test.use = "wilcox", wgcna_name = wgcna_name,
                           barcodes1 = rownames(md)[md$sex == s & md$genotype == cn[1]],
                           barcodes2 = rownames(md)[md$sex == s & md$genotype == cn[2]])
    write.csv(dmes[[nm]], file.path(out_dir, paste0("DMEs_", nm, ".csv")))
  }
}
saveRDS(dmes, file.path(out_dir, "DMEs.rds"))

# ---- GO enrichment of hub genes (Enrichr) ------------------------------------
seu <- RunEnrichr(seu, dbs = enrichr_dbs, max_genes = enrichr_n_genes)
write.csv(GetEnrichrTable(seu), file.path(out_dir, "enrichr_GO_hubgenes.csv"), row.names = FALSE)

# ---- Module-trait correlations per cell type (Fig S12) -----------------------
# Model-style indicators: female and WT are the reference levels
seu$sex_M     <- as.numeric(seu$sex == "M")
seu$Plcg2_HET <- as.numeric(seu$genotype == "HET")
seu$Plcg2_KO  <- as.numeric(seu$genotype == "KO")
seu$sex_HET   <- seu$sex_M * seu$Plcg2_HET   # sex x genotype interactions
seu$sex_KO    <- seu$sex_M * seu$Plcg2_KO

seu <- ModuleTraitCorrelation(seu, traits = c("sex_M", "Plcg2_HET", "Plcg2_KO", "sex_HET", "sex_KO"),
                              cor_method = "pearson", group.by = "celltype")

# One row per cell group x trait x module
mt <- GetModuleTraitCorrelation(seu)
mt_df <- do.call(rbind, lapply(names(mt$cor), function(g) {
  r <- mt$cor[[g]]
  data.frame(group = g, setNames(as.data.frame(as.table(r)), c("trait", "module", "cor")),
             pval = c(mt$pval[[g]][rownames(r), colnames(r)]),
             fdr  = c(as.matrix(mt$fdr[[g]])[rownames(r), colnames(r)]))
}))
write.csv(mt_df, file.path(out_dir, "module_trait_correlation.csv"), row.names = FALSE)

saveRDS(seu, file.path(out_dir, paste0("hdWGCNA_", wgcna_name, ".rds")))
