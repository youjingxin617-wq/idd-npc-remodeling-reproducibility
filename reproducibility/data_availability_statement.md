# Data Availability Statement

## Draft Statement For Manuscript

The datasets analyzed in this study are publicly available from the Gene Expression Omnibus. The main single-cell RNA-seq dataset was GSE165722. External bulk transcriptomic datasets included GSE70362 and GSE147383. Ferroptosis-related driver and suppressor gene sets were obtained from FerrDb. No new human specimens were collected for this study.

Analysis scripts, processed result tables, and figure-generation code should be deposited in a public repository before submission or upon acceptance. Until a public repository is created, this statement should not claim that code is already publicly available.

## Public Data Sources

| Dataset/resource | Role in study | Access route | Notes |
|---|---|---|---|
| GSE165722 | Main human NP scRNA-seq dataset | GEO | Used for single-cell atlas, NPC subclustering, ferroptosis scoring, DEG/enrichment, pseudotime, and CellChat |
| GSE70362 | Main external bulk NP support dataset | GEO | Used for low-grade vs high-grade NP module score comparison |
| GSE147383 | Small supportive bulk NP trend dataset | GEO | Used only as direction-level supportive evidence because sample size is small |
| FerrDb | Ferroptosis driver/suppressor gene sets | FerrDb database | Used for transcriptomic gene-set scoring, not for proving ferroptotic cell death |

## Files Not Recommended For Journal Upload

The following are large intermediate files and should normally not be uploaded as supplementary files:

- Raw GEO archives.
- Processed Seurat RDS objects.
- CellChat RDS objects.
- Slingshot SingleCellExperiment RDS objects.

Instead, provide accession numbers, scripts, and processed result tables.

## Files Recommended For Supplementary Upload

- `08_BMC_submission/supplementary_tables/Supplementary_Table_Index.md`
- `08_BMC_submission/supplementary_tables/Supplementary_Table_S1_NPC_subtype_annotation.csv`
- `08_BMC_submission/supplementary_tables/Supplementary_Table_S2_NPC_subtype_counts_and_proportions.csv`
- `08_BMC_submission/supplementary_tables/Supplementary_Table_S3_FerrDb_gene_sets_used.csv`
- `08_BMC_submission/supplementary_tables/Supplementary_Table_S4_ferroptosis_score_summary_by_subtype_group.csv`
- `08_BMC_submission/supplementary_tables/Supplementary_Table_S5_ferroptosis_score_wilcox_stats.csv`
- `08_BMC_submission/supplementary_tables/Supplementary_Table_S6_DEG_count_summary_by_subtype.csv`
- `08_BMC_submission/supplementary_tables/Supplementary_Table_S7_top30_DEGs_by_subtype_direction.csv`
- `08_BMC_submission/supplementary_tables/Supplementary_Table_S8_top_enrichment_terms.csv`
- `08_BMC_submission/supplementary_tables/Supplementary_Table_S9_GSE70362_module_score_stats.csv`
- `08_BMC_submission/supplementary_tables/Supplementary_Table_S10_GSE147383_supportive_module_score_stats.csv`
- `08_BMC_submission/supplementary_tables/Supplementary_Table_S11_Slingshot_lineages.csv`
- `08_BMC_submission/supplementary_tables/Supplementary_Table_S12_Slingshot_global_pseudotime_stats.csv`
- `08_BMC_submission/supplementary_tables/Supplementary_Table_S13_CellChat_global_network_summary.csv`
- `08_BMC_submission/supplementary_tables/Supplementary_Table_S14_CellChat_centrality_delta_by_subtype.csv`
- `08_BMC_submission/supplementary_tables/Supplementary_Table_S15_CellChat_top_severe_enhanced_pathways.csv`
- `08_BMC_submission/supplementary_tables/Supplementary_Table_S16_CellChat_top_severe_enhanced_LR_pairs.csv`
