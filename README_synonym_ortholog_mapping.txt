
Gene Synonym and Ortholog Mapping Reference (DepMap-anchored)
============================================================

Overview
--------

This folder contains a curated, versioned gene synonym and ortholog reference
table built to support robust mapping of user-provided gene lists (human or mouse)
onto DepMap human genes.

The core design goal is:

    Maximally complete but biologically honest gene name mapping,
    with explicit control over confidence and error propagation.


Backbone principle
------------------

- One row per DepMap gene
- Rows are anchored by `depmap_symbol`
- Rows are never added, removed, or duplicated
- All enrichment happens by adding columns only

This guarantees stability and auditability across versions.


Canonical reference
-------------------

Canonical file:

    depmap_gene_reference_clean_v6.tsv

(and derived v8 with risk tiers)

Number of rows: 44,083
This row count is invariant across all versions.


Core identifier columns
-----------------------

- depmap_symbol
  Canonical DepMap human gene symbol (row anchor)

- locus_group
  Gene biotype (protein-coding, lncRNA, pseudogene, etc.)

- ensembl_gene_id
  Human Ensembl gene ID


Source synonym / ortholog columns
---------------------------------

Each of the following columns represents one independent source of gene naming
or orthology evidence. They are intentionally kept separate to prevent
uncontrolled error propagation.

- all_names_depmap
  DepMap / HGNC gene symbols and aliases (high confidence)

- all_names_MGI_HMD
  Human–mouse orthology derived from MGI Human Phenotype data

- all_names_MRK
  Mouse marker symbols from MGI (historical and alternative mouse symbols)

- all_names_HCOP
  Human–mouse orthologs from HCOP (HGNC Comparison of Orthology Predictions)

- all_names_HCOP_conservative
  Conservative subset of HCOP orthologs

- all_names_human_mouse_error_prone
  Broad human+mouse synonym source known to introduce false positives.
  Retained for transparency and debugging, not for default use.


Risk-tiered aggregate columns (v8)
----------------------------------

These columns combine names already present on the same row.
They do NOT perform recursive expansion or graph closure.

Duplicates are removed string-wise, but case is preserved
(e.g. A1BG and A1bg are both kept).

- all_names_low_risk
  Union of:
    all_names_depmap
    all_names_MGI_HMD
    all_names_MRK

  Intended as the default, high-precision matching set.

- all_names_mid_risk
  Union of:
    all_names_low_risk
    all_names_HCOP_conservative

  Balanced precision and recall.

- all_names_high_risk
  Union of:
    all_names_mid_risk
    all_names_HCOP
    all_names_human_mouse_error_prone

  Maximum recall, known to be error-prone.
  Intended for exploratory analysis and debugging only.


Why risk tiers exist
--------------------

Gene synonym expansion is not monotonic in quality.

Adding more sources increases recall, but also increases false positives,
especially across species.

Rather than collapsing everything into a single column, this reference
preserves provenance and exposes confidence explicitly.


File versioning philosophy
--------------------------

- Each transformation produces a new vX file
- Older versions are kept until explicitly retired
- No version overwrites a previous one
- All changes are additive and traceable

This allows reproducibility, debugging, and safe rollback.


Recommended usage
-----------------

- Default gene matching:
    all_names_low_risk

- Cross-species support enabled:
    all_names_mid_risk

- Exploratory / debug mode:
    all_names_high_risk

Applications should not silently switch between tiers.


Non-goals
---------

This reference does NOT:
- infer paralogy
- resolve gene families
- rank ortholog confidence numerically
- collapse symbols across species automatically

Those decisions are left to downstream logic.


Summary
-------

This folder contains a carefully engineered gene synonym and ortholog mapping
layer that is DepMap-anchored, cross-species aware, and explicitly confidence-
controlled.

The design favors clarity and correctness over convenience.

