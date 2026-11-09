# Plcg2_Glia2026

Code associated with Hartigan et al., 2026 Glia

| Step | Script | What it does |
|------|--------|--------------|
| 1 | `Step1_preprocessing.R` | Per sample: SoupX ambient RNA removal, QC filtering, SCTransform, DoubletFinder |
| 2 | `Step2_integration.R` | Merge samples, keep singlets, SCTransform, CCA integration |
| 3 | `Step3_clustering.R` | Clustering, cluster markers, cluster annotation, cell counts and speckle proportion tests |
| 4 | `Step4_UCell.R` | Homeostatic and inflammatory UCell scores across microglial clusters |
| 5 | `Step5_ADRiskGenes.R` | Microglia-enriched AD GWAS risk genes and their expression across microglial clusters |
| 6 | `Step6_VariancePartition.R` | Microglia pseudobulk variance partition (genotype, sex, genotype x sex; run removed) |
| 7 | `Step7_Pseudobulk.R` | Pseudobulk DESeq2 per cell type: genotype effects (all, F, M) and genotype x sex interaction |
| 8 | `Step8_MASTDEGs.R` | Single-cell MAST DE between genotypes within each sex, per cell type and (KO vs WT) per microglial cluster |
| 9 | `Step9_GSEA.R` | KEGG and Hallmark GSEA on the step 7 results |
| 10 | `Step10_ORA.R` | DEG overlaps, shared-DEG correlation/regression and monotonic genes, KEGG/Hallmark ORA on the step 8 results, mmu05022 log2FC across microglial clusters |
| 11 | `Step11_hdWGCNA.R` | Microglia co-expression modules, hub genes, DMEs, hub-gene GO enrichment (Enrichr), module-trait correlations |
| 12 | `Step12_CellChat.R` | Cell-cell communication for WT and KO in each sex; KO vs WT within each sex |

Steps 4-12 all start from the step 3 object (`output/step3_clustering/Plcg2_annotated.rds`);
9 needs step 7 and 10 needs step 8. Each script writes tables and R objects to its output folder.

```sh
Rscript Step1_preprocessing.R
Rscript Step2_integration.R
Rscript Step3_clustering.R
Rscript Step4_UCell.R
Rscript Step5_ADRiskGenes.R
Rscript Step6_VariancePartition.R
Rscript Step7_Pseudobulk.R
Rscript Step8_MASTDEGs.R
Rscript Step9_GSEA.R
Rscript Step10_ORA.R
Rscript Step11_hdWGCNA.R
Rscript Step12_CellChat.R
```

Parameters (paths, cutoffs, contrasts) are set in the `Config` section at the top of each script.

## Inputs (`config/`)

- `samples.csv`: one row per sample, with columns `ID, sample, genotype, sex, run, cellranger_dir`.
  `genotype` is `WT`/`HET`/`KO` and `sex` is `F`/`M`. `cellranger_dir` is the Cell Ranger `outs`
  folder containing `raw_feature_bc_matrix/`, `filtered_feature_bc_matrix/` and `analysis/`.
- `Annotated_clusters_0.6.xlsx`: one row per cluster at resolution 0.6; columns are
  cluster number, cell type, subtype, cluster label.
- `mh.all.v2026.1.Mm.entrez.gmt`: MSigDB mouse Hallmark gene sets (Entrez IDs), from
  [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/mouse/collections.jsp).

## Requirements

R (>= 4.3) with:

- CRAN: `Seurat` (v5), `readxl`, `WGCNA`, `enrichR`, `future`
- Bioconductor: `glmGamPoi`, `DESeq2`, `edgeR`, `limma`, `variancePartition`, `AnnotationDbi`, `org.Mm.eg.db`, `MAST`,
  `clusterProfiler`, `KEGGREST`, `UCell`, `speckle`
- GitHub: `SoupX` (constantAmateur/SoupX or CRAN), `DoubletFinder` (chris-mcginnis-ucsf/DoubletFinder),
  `hdWGCNA` (smorabit/hdWGCNA), `CellChat` (jinworks/CellChat)

Steps 9-11 need internet access (KEGG and Enrichr are queried at run time).
