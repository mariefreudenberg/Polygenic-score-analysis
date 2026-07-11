<<<<<<< HEAD
# Polygenic Risk Score Analysis

This repository contains the R scripts and R Markdown notebooks used for my master's thesis on polygenic risk scores.
=======
# Polygenic Risk Score Construction and Evaluation Across Ancestries

## Overview

This repository contains the R scripts and analysis notebooks used for my master's thesis investigating the construction, evaluation, and transferability of polygenic risk scores (PRS) across diverse ancestry groups.

The project compares several commonly used PRS construction methods using publicly available GWAS summary statistics, evaluates published PRS from the PGS Catalog, and investigates ancestry-related differences using the 1000 Genomes Project and the LIFE-Adult cohort.

The repository contains only scripts and documentation. Large datasets, genotype files, summary statistics, and intermediate results are intentionally excluded.

---

# Datasets

1000 Genomes Project Phase 3 (https://www.internationalgenome.org/data) used as reference panel for

- PRS calculation
- allele frequency estimation
- ancestry comparisons
- FST calculation

Individual-level genotype and phenotype data from the LIFE-Adult cohort (not publicly available) was used for

- PRS evaluation
- recalibration methods

---

## GWAS summary statistics

Trait-specific GWAS summary statistics were downloaded from the original publications or associated repositories.

Examples include

- GIANT Consortium
- DIAMANTE Consortium
- UK Biobank based GWAS

---

## PGS Catalog

Published polygenic scores were downloaded from

https://www.pgscatalog.org/

---

# PRS construction methods

## LDpred2

Implementation based on the **bigsnpr** R package.

Reference

Privé F, Arbel J, Vilhjálmsson BJ (2020)
LDpred2: better, faster, stronger.

Repository

https://github.com/privefl/bigsnpr

---

## PRS-CS

Implementation of PRS-CS using the official Python software.

Repository

https://github.com/getian107/PRScs

---

## PRSice

PRS generated using PRSice-2.

Repository

https://github.com/choishingwan/PRSice

---

## lassosum

Implementation using the lassosum R package.

Repository

https://github.com/tshmak/lassosum

---

## snpboost

Implementation of component-wise boosting using individual-level genotype data.

Repository

https://github.com/biometrische-gesellschaft/snpboost

Unlike the other methods, snpboost does **not** rely on GWAS summary statistics but instead trains directly on genotype and phenotype data from the LIFE-Adult cohort.

---

# Additional analyses

## PGS Catalog evaluation

Scripts for

- downloading score files
- harmonizing genome builds
- calculating PRS
- ancestry comparisons

---

## FST analyses

Using Hudson's FST to quantify genetic differentiation between ancestry groups in the 1000 Genomes Project.

Analyses include

- pairwise FST
- superpopulation comparisons
- correlation between FST and PRS differences

---

## Manhattan plots

Utilities for visualizing GWAS summary statistics.

---
> Year.
>>>>>>> 7644bbcd680f5f5f402aa543c3a53af6c46657a8
