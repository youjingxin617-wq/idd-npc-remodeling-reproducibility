# Submission Reproducibility Checklist

## Ready

- [x] Public dataset accession numbers identified: GSE165722, GSE70362, GSE147383.
- [x] Main analysis scripts are present locally.
- [x] Main figures are generated and linked to source scripts.
- [x] Figure 7 has been simplified for the main manuscript, with detailed CellChat outputs moved to Supplementary Figure S1.
- [x] Supplementary Tables S1-S16 are assembled as CSV files.
- [x] Supplementary Tables S1-S16 have been combined into one Excel workbook for easier submission.
- [x] `software_versions_sessionInfo.txt` has been generated.
- [x] File manifest has been created.
- [x] Data availability draft has been created.
- [x] Submission upload map has been created.
- [x] Repository README and `.gitignore` templates have been created.

## Still Needed Before Submission

- [x] Public GitHub repository created: https://github.com/youjingxin617-wq/idd-npc-remodeling-reproducibility.
- [ ] Add a license file if scripts are shared publicly.
- [ ] Add a repository-level `README.md` that matches the final manuscript title.
- [ ] Remove or ignore large raw and processed data files before public upload.
- [ ] Confirm all script paths are relative or clearly documented.
- [ ] Add software versions to supplementary methods or repository.
- [ ] Confirm whether BMC requires separate source data files for figures.
- [ ] Confirm author list, affiliations, contributions, funding, and competing-interest statements.
- [ ] Confirm whether the journal wants editable figure files or high-resolution PNG/TIFF files.

## Practical Upload Strategy

For the first supervisor-facing package:

1. Manuscript PDF and DOCX.
2. Supplementary table folder.
3. Reproducibility materials folder: `10_重复性材料`.
4. Main figure PNG files.

For formal journal submission:

1. Manuscript DOCX.
2. Main figures as separate high-resolution files if required.
3. Supplementary tables as CSV or Excel. Current recommended file: `08_BMC_submission/supplementary_tables/Supplementary_Tables_S1_S16.xlsx`.
4. Supplementary Figure S1 if the journal allows separate supplementary figures.
5. Data availability statement in the manuscript.
6. Repository link: https://github.com/youjingxin617-wq/idd-npc-remodeling-reproducibility.

## Caution

Do not overstate reproducibility until the full code package has been cleaned and tested from a fresh directory. The current package is suitable for transparency and supervisor review; a public repository release should receive one more cleanup pass.
