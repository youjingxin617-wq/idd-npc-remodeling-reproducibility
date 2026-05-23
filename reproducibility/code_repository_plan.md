# Code Repository Plan

## Recommended Repository Contents

A public repository for this manuscript should include:

- `README.md`
- `scripts/`
- `results/tables/` or curated result tables
- `results/tables_external_bulk/`
- `results/tables_slingshot/`
- `results/tables_cellchat/`
- `supplementary_tables/`
- `software_versions_sessionInfo.txt`
- `LICENSE`

## Recommended Exclusions

Do not commit large raw or processed files directly to GitHub:

- `data/raw/GSE165722_RAW.tar`
- `data/raw/GSE165722_RAW/*.gz`
- `data/processed/*.rds`
- `data/external_bulk/processed/*.rds`
- large temporary or cache files

Instead, use:

- GEO accession numbers for raw data.
- Scripts to regenerate processed objects.
- Zenodo or OSF only if large processed objects must be shared.

## Suggested `.gitignore`

```gitignore
data/raw/
data/processed/
data/external_bulk/raw/
data/external_bulk/processed/
tmp/
.Rhistory
.RData
.Rproj.user/
*.rds
*.RData
*.tar
*.gz
*.zip
__pycache__/
```

## Suggested Repository Release Before Submission

Before submission, create either:

1. A GitHub repository with scripts and supplementary tables.
2. A Zenodo DOI linked to the GitHub release.

If this is not ready before submission, use conservative wording:

> Analysis scripts and processed result tables will be made available upon reasonable request and will be deposited in a public repository before publication.

Do not write:

> All code is publicly available at ...

until the link exists.

## Current Staging Plan

A clean staging folder should be created outside the analysis workspace. It should include scripts, selected result tables, supplementary tables, reproducibility documentation, `README.md`, `.gitignore`, and a license placeholder. It should exclude raw GEO archives, extracted GEO raw folders, RDS files, IDE folders, and cache files.
