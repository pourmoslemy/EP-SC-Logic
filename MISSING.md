# Exclusions, unverifiable items and release blockers

- Neither related_work_v11.tex nor requested related_work_v10.tex exists; no replacement was invented. The main manuscript has no input directive requiring either file.
- Excluded code EU/ep_experiments/results_v3/report_results.py: embeds identifiable members with individual voting-activity statements; unmodified-copy rule prevents redaction.
- Raw VoteWatch workbooks, export archives, aligned caches, every input NPZ/NPY/PKL, member/group tables, member metadata, and chronological vote-ID tables are excluded. No individual vote arrays are included.
- Corrected EP9 NPZ: 811 x 13459 states/Yes/No/Abstain/presence arrays plus linked member identifiers. Corrected EP6 NPZ: 940 x 6199 equivalent arrays plus linked member identifiers. Both fail the no-per-member-data rule and are excluded; results/input is intentionally absent.
- Availability-text change: the sentence claiming corrected input is in the archive was replaced by a sentence explicitly excluding that per-member input. All other supplied availability wording and placeholders are retained.
- The byte-identical manuscript still contains its original claim that corrected input is included. The author must revise that claim separately; packaging did not edit the manuscript.
- EP9 has no deviations.csv, comparison_old_new.csv, or execution_status.json in its source directory. EP6 versions are included only in results/EP6, never relabelled as EP9.
- Source output manifests are copied unchanged and describe the source runs, including files intentionally excluded here. COPY_MANIFEST.json is the separate package copy manifest.
- Requested Linux mount is unavailable in this Windows workspace; this package was built project-locally, not uploaded to the linked remote location.
- Restricted input, original caches, external original source fingerprints, and raw date metadata are unavailable in the public package. Full numerical reproduction is not certified; no analysis or manuscript compilation was run.
- Author must confirm MIT permission for code, rights for derived data/manuscript/figures, and publisher-class redistribution terms. Existing publisher notices are preserved.
- Open metadata: permanent repository URL, real DOI, full source citation, workbook access dates, author ORCID, and final code-licence choice. No external identifier was created or verified.
- 7 requested input-layout pattern matches remain in unmodified relative paths or provenance text. They are already relative, not machine-absolute paths; original code/results are not rewritten. Additional drive-letter-pattern matches in TeX macros are false positives. See VERIFICATION.json.
- All 45 credential-pattern lexical matches were reviewed: they are publisher-class TeX control sequences, not credential material. Locations are recorded in VERIFICATION.json.
