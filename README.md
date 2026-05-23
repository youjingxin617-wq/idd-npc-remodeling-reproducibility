# NPC state remodeling in intervertebral disc degeneration

This repository contains analysis scripts, processed result tables, supplementary tables, and reproducibility documentation for a public-data transcriptomic study of nucleus pulposus cell state remodeling in intervertebral disc degeneration.

## Study Summary

The analysis reuses public transcriptomic datasets:

- GSE165722: main human nucleus pulposus single-cell RNA-seq dataset.
- GSE70362: external bulk transcriptomic support dataset.
- GSE147383: small-cohort supportive bulk transcriptomic dataset.
- FerrDb: ferroptosis driver and suppressor gene sets.

The study focuses on NPC subtype redistribution, ferroptosis-related transcriptomic scores, severe-versus-mild differential expression, enrichment analysis, exploratory external bulk support, pseudotime inference, and CellChat-based communication inference.

## Repository Structure

```text
scripts/
results/
  tables/
  tables_external_bulk/
  tables_slingshot/
  tables_cellchat/
supplementary_tables/
reproducibility/
software_versions_sessionInfo.txt
```

## How To Reproduce

1. Install R 4.4.1 or a compatible recent R version.
2. Install the R packages listed in `software_versions_sessionInfo.txt`.
3. Download public GEO data using `scripts/01_download_GSE165722.ps1`.
4. Run scripts in the order described in `reproducibility/analysis_workflow.md`.
5. Compare regenerated tables with the files in `results/` and `supplementary_tables/`.

## Main Workflow

The main script order is:

1. `scripts/01_download_GSE165722.ps1`
2. `scripts/02_reproduce_GSE165722_scRNA.R`
3. `scripts/05_apply_manual_celltype_labels.R`
4. `scripts/06_reproduce_article_npc_subset.R`
5. `scripts/07_clean_npc_subset_plot.R`
6. `scripts/08_finalize_npc_names_and_ferroptosis_plots.R`
7. `scripts/09_ferroptosis_driver_suppressor_scores.R`
8. `scripts/10_deg_and_enrichment_fig2f.R`
9. `scripts/13_external_bulk_validation_step1_download_inspect.R`
10. `scripts/14_external_bulk_GSE70362_annotate_validate.R`
11. `scripts/15_external_bulk_GSE70362_exploratory_plots.R`
12. `scripts/16_external_bulk_GSE147383_inspect_validate.R`
13. `scripts/17_external_bulk_GSE147383_supportive_plots.R`
14. `scripts/18_external_bulk_combined_figure5.R`
15. `scripts/21_exploratory_pseudotime_mst.R`
16. `scripts/22_slingshot_pseudotime.R`
17. `scripts/23_cellchat_npc_communication.R`

Manuscript-generation scripts are retained for traceability but are not required to reproduce the biological analyses.

## Data Availability

Raw datasets are publicly available from GEO. Large intermediate objects such as Seurat RDS files are not included in this repository because they can be regenerated from public data using the scripts above.

## Code Availability

All analysis scripts used for the current manuscript figures and processed result tables are included in this repository. If the repository is archived with Zenodo, cite the Zenodo DOI in the manuscript reference list.

## License

Add a license before public release. For code, MIT is a common permissive option. For documentation and processed tables, CC BY 4.0 can be considered if acceptable to all authors.
