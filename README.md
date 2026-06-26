# MAVIS — Mitochondrial Analysis, Variant Inspection & Summary

An integrated pipeline for mitochondrial DNA (mtDNA) variant calling, cohort-level quality control reporting, and curation from whole-genome sequencing (WGS) data. Designed for mitochondrial disease and rare disease research.


---

## Overview

The pipeline has three main modules:

| Module | Purpose |
|--------|---------|
| [`mtDNA_calling`](#1-mtdna_calling) | Batch variant calling via MitoHPC (Mutect2) on HPC |
| [`mtDNA_report`](#2-mtdna_report) | Cohort-level QC and summary HTML report |
| [`mtDNA_curation`](#3-mtdna_curation) | Interactive Shiny app for multi-sample variant review |

**Reference databases used:** gnomAD v3.1 (mtDNA), MITOMAP, mitoTIP

---

## Requirements

- R ≥ 4.2
- SLURM-based HPC (configurable for SGE/bash)
- MitoHPC (installed via `mtDNA_calling/install.sh`)
- Reference genome: GRCh38/hs38DH (with decoy)

### R packages

```r
# Core
install.packages(c("data.table", "tidyverse", "ggplot2", "plotly", "patchwork",
                   "DT", "shiny", "openxlsx", "Rsamtools"))

# Bioconductor
BiocManager::install("VariantAnnotation")
```

---

## Quick Start

### 1. Install MitoHPC

```bash
cd mtDNA_calling
bash install.sh
```

### 2. Generate input manifest

```r
# Edit paths in generate_input_list.R, then:
Rscript mtDNA_calling/generate_input_list.R
```

### 3. Submit jobs

```r
Rscript mtDNA_calling/run_multisample.R
```

### 4. Check job completion

```r
Rscript mtDNA_calling/check_jobs.R
```

### 5. Generate cohort report

```r
# Edit mtDNA_report/config.R, then:
Rscript mtDNA_report/render.R
```

### 6. Launch curation app

```r
shiny::runApp("mtDNA_curation/app.R")
```

---

## Module Details

### 1. `mtDNA_calling`

Wraps [MitoHPC](https://github.com/dpaulson45/mitoHPC) to run Mutect2 in mtDNA mode across a cohort.

**Scripts:**

| Script | Description |
|--------|-------------|
| `install.sh` | Clones and installs MitoHPC |
| `generate_input_list.R` | Scans BAM directories, validates files, writes input manifest |
| `run_multisample.R` | Generates and submits per-sample SLURM jobs |
| `check_jobs.R` | Identifies failed/incomplete samples, writes rerun script |

**Per-sample outputs** (in `out/<sample>/`):

| File | Description |
|------|-------------|
| `.mutect2.00.vcf` | Mutect2 variant calls |
| `.mutect2.haplogroup` | Haplogroup assignment |
| `.mutect2.haplocheck` | Contamination estimate |
| `.mutect2.cvg.stat` | Coverage statistics (median, mean, min, max) |
| `.mutect2.cvg` | Per-position depth |
| `count.tab` | mtDNA copy number |

**Key parameters** (set in `init.sh`):

```bash
HETEROPLASMY_FILTER=0.01     # minimum heteroplasmy threshold
REFERENCE=hs38DH             # GRCh38 with decoy
CALLER=mutect2               # mutect2 | freebayes | mutserve
```

Filtered flags: `strand_bias`, `weak_evidence`, `strict_strand`, `base_qual`, `position`

---

### 2. `mtDNA_report`

Generates a self-contained interactive HTML cohort report.

**Configuration:** Edit `mtDNA_report/config.R` before rendering.

**Report sections:**

| Section | Tabs |
|---------|------|
| Coverage | Depth per sample, mtDNA-CN, base coverage, cumulative coverage, CV distribution |
| Contamination | Level per sample, level vs depth |
| Haplogroup | Macro-haplogroup distribution, assignment quality |
| PCA | PC1 vs PC2, PC2 vs PC3, scree plot |
| Variants | Summary, VAF distribution, homoplasmy by clade, heteroplasmy by clade, MSS, consequence breakdown, MLC score, frequency spectrum |
| Relatedness | Pairwise Jaccard similarity, flagged pairs (conditional) |
| QC Summary | Flag overview, per-sample table |
| Glossary | Definitions for VAF, CV, haplogroup, heteroplasmy, Jaccard, MLC, MSS, mtDNA-CN, NUMT |

**QC thresholds** (configurable in `config.R`):

```r
MIN_MEDIAN_DEPTH    <- 100    # flag samples below this
CONTAMINATION_CUTOFF <- 0.02  # flag samples at or above 2%
HOMOPLASMY_VAF      <- 0.95   # VAF threshold for homoplasmic calls
```

**Outputs:**
- `mtDNA_report_<date>.html` — self-contained interactive report
- `qc_summary.tsv` — per-sample QC metrics

---

### 3. `mtDNA_curation`

Interactive Shiny app for reviewing, filtering, and annotating variants across a cohort.

**Setup:**

```r
# 1. Build gnomAD haplogroup-stratified allele frequency table (once only)
Rscript mtDNA_curation/prepare_databases.R

# 2. Generate sample manifest
Rscript mtDNA_curation/generate_manifest.R

# 3. Launch app
shiny::runApp("mtDNA_curation/app.R")
```

**Genotype classification:**

| Class | VAF range |
|-------|-----------|
| Homoplasmic | ≥ 0.95 |
| Heteroplasmic | 0.03 – 0.95 |
| Low heteroplasmy | < 0.03 |

**Filter options:**

- VAF threshold (default minimum 3%)
- gnomAD allele frequency cutoff
- PASS-only filter
- Strand bias artefact blacklist
- Hypervariable region exclusion
- Synonymous variant exclusion
- Region / gene / sample ID filter
- Maternal inheritance (optional FAM file)

**Blacklisted positions** (known artefacts):  
301, 302, 310, 316, 3107, 5894, 10933, 16179, 16181–16183, 16188, 16189, 16192

**Annotations per variant:**

- VAF, depth, genotype class
- Haplogroup assignment (per sample)
- gnomAD v3.1 allele frequency (global and per-haplogroup)
- MITOMAP disease associations (confirmed and reported)
- mitoTIP pathogenicity score

---

## Reference Databases

Place files under `database/`:

```
database/
├── gnomad/
│   ├── gnomad.genomes.v3.1.sites.chrM.vcf.bgz        # download from gnomAD
│   ├── gnomad.genomes.v3.1.sites.chrM.vcf.bgz.tbi
│   └── gnomad_chrM_hap_AF.tsv.gz                      # generated by prepare_databases.R
├── MITOMAP/
│   ├── MutationsCodingControl_*_VAR.csv               # processed with add_VAR_MitoMAP.R
│   └── MutationstRNA_*_VAR.csv
└── MitoTIP/
    └── mitotip_scores_27_04_2020_VAR.txt              # processed with add_VAR_MitoTIP.R
```

**Sources:**
- gnomAD mtDNA v3.1: https://gnomad.broadinstitute.org/downloads#v3-mitochondrial-dna
- MITOMAP: https://www.mitomap.org/foswiki/bin/view/MITOMAP/Resources
- mitoTIP: https://github.com/sonnymijael/mitoTIP

After downloading MITOMAP and mitoTIP files, run the helper scripts to add the `VAR` join key (format: `position_REF_ALT`):

```r
Rscript database/MITOMAP/add_VAR_MitoMAP.R
Rscript database/MitoTIP/add_VAR_MitoTIP.R
```

---

## Variant Representation

All databases and outputs use a `VAR` column as the join key:

```
position_REF_ALT    e.g., 8573_G_A, 3243_A_G, 16189_T_C
```

Based on the revised Cambridge Reference Sequence (rCRS, NC_012920.1).

---

## Glossary

| Term | Definition |
|------|-----------|
| mtDNA | Mitochondrial DNA (~16.5 kb, circular, maternally inherited) |
| rCRS | Revised Cambridge Reference Sequence (NC_012920.1) |
| Heteroplasmy | Mixture of mtDNA variants within a cell/tissue |
| Homoplasmy | Single mtDNA sequence at high frequency (≥ 95%) |
| Haplogroup | Phylogenetic mtDNA lineage defined by characteristic variants |
| NUMT | Nuclear mitochondrial DNA insert |
| MitoHPC | Mutect2-based mtDNA variant calling pipeline |
| Haplocheck | Contamination detection tool for mtDNA |

See [`disease_glossary.md`](disease_glossary.md) for mitochondrial disease abbreviations.

---

## Citation

If you use this pipeline, please cite:

- **MitoHPC:** Paulson et al. (cite relevant publication)
- **gnomAD mtDNA:** Karczewski et al., *Nature* 2020; Laricchia et al., *Genome Research* 2022
- **MITOMAP:** https://www.mitomap.org
- **mitoTIP:** Sonney et al., *PLOS Genetics* 2017
- **Haplocheck:** Weissensteiner et al., *Genome Research* 2021

---

## Contact

Wang Lo — wang.lo@wehi.edu.au  
Walter and Eliza Hall Institute of Medical Research
