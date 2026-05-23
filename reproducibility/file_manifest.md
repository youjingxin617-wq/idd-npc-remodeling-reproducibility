# File Manifest

This manifest links the manuscript figures and tables to their source scripts and result files.

## Manuscript Files

| File | Purpose |
|---|---|
| `08_BMC_submission/BMC_Musculoskeletal_Disorders_oriented_manuscript.md` | Markdown source of current BMC-oriented draft |
| `08_BMC_submission/BMC_Musculoskeletal_Disorders_oriented_manuscript.docx` | Word version |
| `08_BMC_submission/BMC_Musculoskeletal_Disorders_oriented_manuscript.pdf` | PDF version |

## Main Figures

| Figure | Current file | Source script |
|---|---|---|
| Figure 1 | `03_图件总审查/main_figure_candidates/Figure1_global_scRNA_atlas_manual_celltypes_candidate.png` | `scripts/05_apply_manual_celltype_labels.R` |
| Figure 2 | `03_图件总审查/main_figure_candidates/Figure2_NPC_subtypes_mild_vs_severe_candidate.png` | `scripts/08_finalize_npc_names_and_ferroptosis_plots.R` |
| Figure 3 | `03_图件总审查/main_figure_candidates/Figure3_ferroptosis_scores_candidate.png` | `scripts/09_ferroptosis_driver_suppressor_scores.R` |
| Figure 4 | `03_图件总审查/main_figure_candidates/Figure4_DEG_enrichment_candidate.png` | `scripts/10_deg_and_enrichment_fig2f.R` |
| Figure 5 | `01_外部bulk验证/figures/Figure5_external_bulk_module_validation_main_candidate.png` | `scripts/18_external_bulk_combined_figure5.R` |
| Figure 6 | `05_伪时序分析/figures/56_npc_slingshot_main_Figure6_candidate.png` | `scripts/22_slingshot_pseudotime.R` |
| Figure 7 | `06_细胞通讯分析/figures/63_npc_cellchat_main_Figure7_candidate.png` | `scripts/23_cellchat_npc_communication.R` |

## Main Tables

| Main table | Source file | Source script |
|---|---|---|
| Table 1 | `results/tables/npc_final_subtype_counts_and_proportions.csv` | `scripts/08_finalize_npc_names_and_ferroptosis_plots.R` |
| Table 2 | `results/tables/clean_npc_deg_summary_by_subtype.csv` | `scripts/10_deg_and_enrichment_fig2f.R` |
| Table 3 | `results/tables_external_bulk/GSE70362_NP_low_high_module_score_stats.csv` | `scripts/14_external_bulk_GSE70362_annotate_validate.R` |

## Supplementary Tables

The curated supplementary table package is located at:

`08_BMC_submission/supplementary_tables`

See:

`08_BMC_submission/supplementary_tables/Supplementary_Table_Index.md`

## Large Intermediate Files

Large RDS files in `data/processed` are useful for local reproduction but should not be uploaded as ordinary supplementary files unless specifically requested.
