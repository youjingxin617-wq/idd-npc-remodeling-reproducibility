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

## Current Repository Release

The public GitHub repository has been created:

https://github.com/youjingxin617-wq/idd-npc-remodeling-reproducibility

Current manuscript wording may state that analysis scripts, processed result tables, supplementary tables, and software/session information are available at this repository.

Optional next step: create a Zenodo DOI linked to a GitHub release after the author team confirms the repository contents.

## Current Staging Plan

A clean public repository folder has been created outside the analysis workspace. It includes scripts, selected result tables, supplementary tables, reproducibility documentation, `README.md`, `.gitignore`, and a license placeholder. It excludes raw GEO archives, extracted GEO raw folders, RDS files, IDE folders, and cache files.
