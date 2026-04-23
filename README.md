# Biol7180 Code Review 
### Kira Noordwijk

## Scripts to set up for and down stream analysis for BayesPrism Workflow created by Danko-Lab (https://github.com/Danko-Lab/BayesPrism) to perform RNA bulk sequencing deconvolution.

### Using remotes:
```r
remotes::install_github("kjn0012/Biol7180_Class_Project")
```

This repository contains a complete pipeline to run **BayesPrism** on your own bulk RNA-seq and scRNA-seq data, including:

* Data loading
* QC of cell type/state labels
* Gene filtering
* Prism object construction
* Running BayesPrism
* Downstream analysis

---

# Package Requirements

### R packages

Install required packages in R:

```r
install.packages("devtools")
install.packages("BiocManager")

# BayesPrism
devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")

# Downstream analysis
BiocManager::install("DESeq2")
```

---

## Input files

You need **3 input files**:

### 1. Bulk RNA-seq (`bulk_counts.tsv`) - a minimum of 4 sample sets is recommended

* Rows = samples
* Columns = genes
* Values = raw counts (preferred) or TPM/RPKM

### 2. Single-cell RNA-seq (`sc_counts.tsv`) - this is the reference sample (further recommendations can be found in Dank_Lab files vignette)

* Rows = cells
* Columns = genes
* Values = raw counts

### 3. Metadata (`sc_metadata.tsv`) - also part of the reference file

* Rows = cells (must match `sc_counts.tsv`)
* Required columns:

  * `cell_type` (broad types)
  * `cell_state` (subclusters or finer labels)

---

### Important requirements

* Gene names must match between bulk and scRNA-seq
* Use **raw counts if possible** (avoid log-transformed data)
* `nrow(sc.dat)` must equal:

  * `length(cell.type.labels)`
  * `length(cell.state.labels)`

---

# Pipeline overview

Run scripts in this order:

```bash
Rscript load_data.R
Rscript qc_labels.R
Rscript filter_genes.R
Rscript construct_prism.R
Rscript run_bayesprism.R
Rscript downstream_analysis.R
```

---

# Step-by-step description

## 1. Load data (`load_data.R`)

* Reads input files
* Aligns metadata
* Performs basic checks

---

## 2. QC of labels (`qc_labels.R`)

* Plots correlation heatmaps:

  * cell types
  * cell states
* Flags low-cell-count groups

Purpose: Detect poorly defined or redundant clusters

---

## 3. Filter outlier genes (`filter_genes.R`)

Removes problematic genes:

* Ribosomal genes (`Rb`, `Mrp`)
* Mitochondrial genes (`chrM`)
* Sex chromosome genes (`chrX`, `chrY`)
* MALAT1
* Low-expression genes

Outputs:

* `sc.dat.filtered`
* `sc.dat.filtered.pc` (protein-coding only)
* optional `sc.dat.filtered.pc.sig` (marker genes)

---

## 4. Construct prism object (`construct_prism.R`)

Creates:

```r
myPrism <- new.prism(...)
```

Key parameters:

* `reference` → filtered scRNA matrix
* `mixture` → bulk matrix
* `key` → malignant cell type (e.g. `"tumor"` or `NULL`)

---

## 5. Run BayesPrism (`run_bayesprism.R`)

```r
bp.res <- run.prism(myPrism)
```

Output:

* `bp.res.rds`

Default parameters are recommended but can be altered to fit needs

---

## 6. Downstream analysis (`downstream_analysis.R`)

### Outputs:

#### Cell fractions

* `theta.tsv`
* clustering: `theta_clustering.pdf`

#### Tumor expression (Z)

* `Z_raw.rds`
* `Z_vst.tsv`
* clustering: `Z_clustering.pdf`

#### Signature scores

* `signature_scores.tsv`

#### Correlations

* `cor_Z_vs_<celltype>.tsv`

#### Clinical-ready table

* `theta_for_clinical_analysis.tsv`

---

# Downstream analysis options

## 1. Cluster samples

* By **cell fractions (`theta`)**
* By **tumor expression (`Z`)**

## 2. Signature scoring

Compute z-scores for gene sets:

```r
signature.score <- rowMeans(scale(Z_subset))
```

---

## 3. Correlation analysis

Correlate tumor expression with microenvironment:

* tumor gene expression (Z)
* vs non-malignant cell fractions (theta)

Useful for:

* pathway discovery
* microenvironment interactions

---

## 4. Survival / clinical analysis

Use:

```
theta_for_clinical_analysis.tsv
```

Merge with clinical metadata and analyze in R:

* Cox regression
* Kaplan-Meier curves

---

# ⚙️ Key parameters explained

### `input.type`

* `"count.matrix"` → recommended (scRNA counts)
* `"GEP"` → for precomputed expression profiles

### `key`

* malignant cell type (e.g. `"tumor"`)
* use `NULL` if:

  * no malignant cells
  * matched tumor reference

### `outlier.cut` / `outlier.fraction`

* filters extreme genes in bulk
* default values usually fine

---

# Troubleshooting

### “no package called BayesPrism”

→ install with `devtools::install_github()`

### dimension mismatch

→ metadata not aligned with scRNA matrix

### few shared genes

→ inconsistent gene naming (Ensembl vs symbols)

### bad clustering / noisy results

→ likely:

* poor cell labels
* unfiltered genes
* batch effects

---

# Tips

* Start simple:

  * use `cell.state = cell.type`
* Only add states when biologically meaningful
* Use **filtered + protein-coding genes** for best results
* Marker selection is optional but helpful in noisy datasets

---

# Minimal working example

```r
myPrism <- new.prism(
  reference = sc.dat.filtered.pc,
  mixture = bk.dat,
  input.type = "count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key = "tumor"
)

bp.res <- run.prism(myPrism)
```

---


## Project Proposal
The purpose of this study is to perform Bayesian deconvolution methods on bulk RNA sequencing data from equine sarcoid tumor samples.

### Aim
Adapt exisitng code for the application of Bayesian modelling for RNA sequencing in oncology to bulk RNA sequencing data in equine sarcoid tumors (and hopefully sn-RNA-seq if we get it to work asap)

### Source
Chu, T., Wang, Z., Pe’er, D. et al. Cell type and gene expression deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and single-cell RNA sequencing in oncology. Nat Cancer 3, 505–517 (2022). https://doi.org/10.1038/s43018-022-00356-3

GitHub: https://github.com/Danko-Lab/BayesPrism.git

## Data
The data was collect retrospectively from the hospital electronic database and electronic surgical records. Additional data was collect from published material to increase data volume.
