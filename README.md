# Beyond the Grand Coalition: Simplicial Heyting Semantics and Higher-Order Coalition Capacity in the European Parliament

Author: Alireza Pourmoslemi; Payame Noor University.
Version: 0.1.0. Release status: preparation draft; see MISSING.md.

## Abstract (verbatim from the v11 manuscript)

Pairwise party separation and higher-order legislative capacity need not move together. We apply simplicial Heyting semantics (SC-Logic) to European Parliament roll calls, tracking weighted Newman--Girvan party modularity and clique-based depth over thresholded co-Yes graphs. In the revised EP9 case study (2019--June 2022), modularity rises as the maximum face dimension falls at the selected local step. This is an exploratory finding: it does not identify the cause of the decline, and the location of the drop depends on the threshold mesh. Supplementary EP6 results from the 2004--2009 term, which use a legacy normalization and dimension limit, show the same direction without a comparable null test.

## Data and code availability

Data and code availability. The raw VoteWatch Europe roll-call data used in this study are not redistributed. The derived EP9 records and analysis code are available at [permanent repository URL] (DOI: [real DOI]). The archive contains party-aware null outputs, but excludes the corrected per-member input and the raw VoteWatch files. The code is released under the [MIT License / chosen code licence]. The source data were obtained from [full VoteWatch source citation, workbook names, and access dates].

## Contents and provenance

paper/paper_EU_v11.tex and paper/references_v5.bib are byte-identical source copies.
paper/sc_logic_figures contains the four figures referenced by the manuscript.
code/ep_experiments preserves the original analysis directory structure and unmodified Python source.
results/observed.csv, results/main_steps.csv, results/density_at_050.csv and results/null_pvalues.csv are EP9 records.
results/mesh.csv and results/bootstrap_CI.csv contain the EP9 mesh and approximate-bootstrap records.
results/EP6 contains separately labelled EP6 outputs, including approved Stage-2 null results.
Source output manifests are not package completeness manifests; use COPY_MANIFEST.json.
No raw workbooks, individual vote arrays, or member-linked tables are distributed.

## Reproduction

Install dependencies from requirements.txt in an isolated Python environment.
Dependency versions describe the inspected environment, not a verified original-run lockfile.
The supplied derived results can be inspected without acquiring restricted input.
To rerun numerical analyses, first obtain permission and the source workbooks described in data/README.md.
Rebuild corrected input privately using the supplied input-construction code; never add that input to this public package.
The unchanged historical entry points use their original relative input/cache conventions and source-fingerprint checks. They are not turnkey public-package commands.
Review code/ep_experiments/results_v3/common.py and its imports, then privately restore the original input/cache structure before running observed, refinement, null and optional bootstrap stages.
The EP6 code and data are supplementary and use their own output folder; preserve the recorded refinement guard failure.
Analysis scripts were not executed for packaging. No numerical values or copied source files were edited.
Compile the manuscript from paper/ using a suitable TeX installation; bibliography style is cas-model2-names.

## Limitations and release blockers

See MISSING.md and VERIFICATION.json. Individual corrected input is deliberately absent.
Depth timeout bounds depend on runtime hardware/scheduling. Bootstrap depth is approximate.
The manuscript availability sentence still claims corrected input is included; it requires an author edit before release.
The manuscript has no separate supplied related-work file, and public-input reconstruction is not certified.
MIT licensing is provisional pending author confirmation; it does not grant rights to VoteWatch source data or third-party publisher assets.

## Citation

Use CITATION.cff after confirming author ORCID, [permanent repository URL], and [real DOI].
No repository, deposit, or DOI has been created by this packaging operation.
