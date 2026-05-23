# Reproducibility README

Project: Public single-cell and bulk transcriptomic analysis of NPC state remodeling in intervertebral disc degeneration

Date prepared: 2026-05-23

## Purpose

This folder documents how the analysis underlying the BMC Musculoskeletal Disorders-oriented draft can be reproduced from public data and local scripts.

The package is intended for:

1. Internal lab review.
2. Supervisor review.
3. Future GitHub or Zenodo deposition.
4. Preparing the Data availability and Code availability sections for manuscript submission.

## Current Status

This is a reproducibility preparation package, not yet a public repository release.

The manuscript analysis is currently based on:

- GSE165722 single-cell RNA-seq data.
- GSE70362 public bulk NP transcriptomic data.
- GSE147383 public bulk NP transcriptomic data used as small-cohort supportive trend evidence.
- FerrDb-derived ferroptosis driver and suppressor gene sets.

## What Should Be Uploaded With Submission

Usually upload as supplementary files:

- Main supplementary table package: `08_BMC_submission/supplementary_tables`.
- Preferred single-file upload if the submission system allows Excel: `Supplementary_Tables_S1_S16.xlsx`.
- Any required reporting checklist if requested by the journal.

Usually provide through GitHub/Zenodo or upon request:

- Analysis scripts.
- Reproducibility README.
- File manifest.
- Software/session information.

Usually do not upload because the data are public and large:

- Raw GEO archives.
- Large processed Seurat RDS objects.
- Large CellChat or Slingshot RDS objects.

Instead, provide accession numbers and code to regenerate them.

## Minimal Reproducibility Claim

Recommended manuscript wording:

> The public datasets analyzed in this study are available from the Gene Expression Omnibus under accession numbers GSE165722, GSE70362, and GSE147383. Analysis scripts, processed result tables, and figure-generation code will be made available in a public repository upon acceptance or before submission.

Use stronger wording only after the repository is actually created.

## Current Practical Decision

For this manuscript, reproducibility materials should be prepared before submission. The journal-facing upload should include the supplementary tables and the Data availability statement. The code-facing materials should be deposited in a public repository or kept ready for immediate release during review.
