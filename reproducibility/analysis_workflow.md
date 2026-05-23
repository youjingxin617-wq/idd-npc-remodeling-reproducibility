# Analysis Workflow

This workflow records the script order used to generate the current manuscript figures, tables, and supplementary materials.

Project root:

Run all commands from the repository root after cloning or downloading this repository.

Rscript:

Use R 4.4.1 or a compatible recent R version. The original analysis used R 4.4.1 on Windows 11; selected package versions are listed in `software_versions_sessionInfo.txt`.

## Stage 1. Main GSE165722 scRNA-seq Processing

| Step | Script | Purpose | Main outputs |
|---|---|---|---|
| 1 | `scripts/01_download_GSE165722.ps1` | Download and extract GSE165722 raw GEO files | `data/raw/GSE165722_RAW`, `GSE165722_sample_metadata.csv` |
| 2 | `scripts/02_reproduce_GSE165722_scRNA.R` | Build Seurat object, QC, normalization, Harmony integration, clustering | Global UMAP/QC figures, clustered RDS |
| 3 | `scripts/05_apply_manual_celltype_labels.R` | Assign major cell-type labels | Figure 1 source image; labeled Seurat object |
| 4 | `scripts/06_reproduce_article_npc_subset.R` | Extract and recluster NPC-like/stromal compartment | NPC reclustered object and figures |
| 5 | `scripts/07_clean_npc_subset_plot.R` | Clean NPC subset and assign preliminary NPC labels | Cleaned NPC subset and preliminary subtype figures |
| 6 | `scripts/08_finalize_npc_names_and_ferroptosis_plots.R` | Final NPC naming, including `Hom/Reg-like NPCs`; generate Figure 2 source | Final named NPC RDS; subtype annotation and proportion tables |

## Stage 2. Ferroptosis Scores, DEG, and Enrichment

| Step | Script | Purpose | Main outputs |
|---|---|---|---|
| 7 | `scripts/09_ferroptosis_driver_suppressor_scores.R` | Calculate FerrDb driver/suppressor module scores | Figure 3 source; FDS/FSS summary/stat tables |
| 8 | `scripts/10_deg_and_enrichment_fig2f.R` | Subtype-resolved severe vs mild DEG and Enrichr analysis | Figure 4 source; DEG and enrichment tables |

## Stage 3. External Bulk Transcriptomic Support

| Step | Script | Purpose | Main outputs |
|---|---|---|---|
| 9 | `scripts/13_external_bulk_validation_step1_download_inspect.R` | Download/inspect candidate bulk datasets | Raw and phenotype inspection tables |
| 10 | `scripts/14_external_bulk_GSE70362_annotate_validate.R` | Prepare GSE70362 NP samples and module scores | GSE70362 module statistics |
| 11 | `scripts/15_external_bulk_GSE70362_exploratory_plots.R` | Plot GSE70362 module support | GSE70362 module figures |
| 12 | `scripts/16_external_bulk_GSE147383_inspect_validate.R` | Inspect and prepare GSE147383 | GSE147383 processed expression and annotations |
| 13 | `scripts/17_external_bulk_GSE147383_supportive_plots.R` | Generate small-cohort supportive plots | GSE147383 support tables and figures |
| 14 | `scripts/18_external_bulk_combined_figure5.R` | Combine external bulk module plots | Figure 5 source |

## Stage 4. Pseudotime and CellChat

| Step | Script | Purpose | Main outputs |
|---|---|---|---|
| 15 | `scripts/21_exploratory_pseudotime_mst.R` | Exploratory MST pseudotime | MST pseudotime tables and figures |
| 16 | `scripts/22_slingshot_pseudotime.R` | Slingshot pseudotime and module trends | Figure 6 source; Slingshot tables |
| 17 | `scripts/23_cellchat_npc_communication.R` | Exploratory CellChat communication inference | Figure 7 source; CellChat tables |

## Stage 5. Manuscript and Supplementary Materials

| Step | Script or location | Purpose | Main outputs |
|---|---|---|---|
| 18 | `scripts/25_generate_bmc_musculoskeletal_version.py` | Generate BMC-oriented manuscript MD/DOCX/PDF | `08_BMC_submission` manuscript files |
| 19 | `08_BMC_submission/supplementary_tables` | Curated supplementary table package | Supplementary Tables S1-S16 |

## Notes

Some early draft-generation scripts are retained for provenance but are not part of the current BMC manuscript-generation chain:

- `scripts/11_generate_sci_manuscript_doc.py`
- `scripts/12_generate_integrated_sci_manuscript_doc.py`
- `scripts/20_generate_results_draft_doc_pdf.py`
- `scripts/24_generate_initial_full_manuscript_doc_pdf.py`
- `scripts/26_generate_mentor_reading_package.py`
