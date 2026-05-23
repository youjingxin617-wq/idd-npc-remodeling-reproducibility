# Public Repository Staging Manifest

Date prepared: 2026-05-23

Staging folder:

`10_重复性材料/public_repository_staging`

## Purpose

This folder is a clean pre-publication package for a future GitHub or Zenodo-linked repository. It is separate from the working analysis directory so that large raw data and intermediate RDS objects are not accidentally uploaded.

## Current Contents

| Folder/file | Contents | Intended use |
|---|---|---|
| `README.md` | Repository README | Public repository landing page |
| `.gitignore` | Exclusion rules for raw data, RDS files, caches, and local IDE files | Prevents accidental upload of large/private files |
| `software_versions_sessionInfo.txt` | R session and selected package versions | Computational reproducibility |
| `scripts/` | R and PowerShell analysis scripts | Reproduce data processing, scoring, figures, tables |
| `results/tables/` | Main GSE165722 scRNA-seq result tables | Figure/table traceability |
| `results/tables_external_bulk/` | GSE70362 and GSE147383 result tables | External bulk support traceability |
| `results/tables_slingshot/` | Slingshot pseudotime result tables | Exploratory pseudotime traceability |
| `results/tables_cellchat/` | CellChat communication result tables | Exploratory communication traceability |
| `supplementary_tables/` | Supplementary Tables S1-S16 as CSV and Excel workbook | Journal upload and repository copy |
| `reproducibility/` | Workflow, data availability, manifest, upload map, and checklist | Documentation |

## Current Package Check

The repository staging folder currently contains 119 files. Hard-coded local absolute paths were removed from the runnable scripts and repository-facing workflow document. The public repository URL has been added to the manuscript and reproducibility documentation. The main remaining pre-release item is the license decision.

## Exclusion Check

The staging folder was checked for common large intermediate file types:

- `.rds`
- `.RData`
- `.tar`
- `.gz`
- `.zip`

No such files were present at the time of staging.

## Before Public Release

Complete these steps before making the repository public:

1. Rename the repository to match the final manuscript title.
2. Review `README.md` for final title and author-approved wording.
3. Add a license file after author agreement.
4. Confirm that all scripts run from relative paths or clearly document required working directory.
5. Confirm the license with all authors.
6. Optionally create a Zenodo release DOI after author approval.
