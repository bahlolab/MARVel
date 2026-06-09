# mtDNA_report

Cohort-level mitochondrial DNA variant report. Takes [mitoHPC](https://github.com/bahlolab/mitoHPC) output and produces an interactive HTML report covering coverage QC, contamination, haplogroup assignment, PCA, variant calling, relatedness, and a summary QC table.

---

## Requirements

R ≥ 4.2 with the following packages:

```r
install.packages(c(
  "rmarkdown", "knitr", "ggplot2", "plotly", "patchwork",
  "DT", "htmltools", "scales", "data.table"
))

# BEDMatrix for PCA projection
install.packages("BEDMatrix")
```

---

## Directory structure

```
mtDNA_report/
├── config.R          # all user-facing settings — edit this before running
├── report.Rmd        # main report template
├── modules/
│   ├── 01_load_data.R
│   ├── 02_coverage.R
│   ├── 03_contamination.R
│   ├── 04_haplogroup.R
│   ├── 05_pca.R
│   ├── 06_variants.R
│   ├── 07_relatedness.R
│   └── 08_qc_summary.R
└── 1000G/            # PCA reference files (see below)
    ├── 1000G_mtPCApanel.bed/bim/fam
    ├── mtPCA.variant_ids.txt
    ├── mtPCA.mu.txt
    ├── mtPCA.loadings.PC1-20.txt
    └── mtPCA.eigenvalues.txt
```

The `1000G/` reference files are fixed and should not be modified.

---

## Input: mitoHPC output

Each sample must have the following files under `<results_root>/<ID>/mitoHPC/out/<inner_ID>.merged/`:

| File | Content |
|------|---------|
| `<ID>.merged.mutect2.cvg.stat` | Coverage statistics (median, mean, min, max) |
| `<ID>.merged.mutect2.cvg` | Per-position depth (3-column: chrom, pos, depth) |
| `<ID>.merged.mutect2.00.vcf` | Mutect2 VCF (mitochondrial mode) |
| `<ID>.merged.mutect2.haplogroup` | mitoHPC haplogroup assignment |
| `<ID>.merged.mutect2.haplocheck` | Haplocheck contamination output |
| `<inner_ID>/mutect2.haplocheck.tab` | Haplocheck extended table |

mtDNA copy number is read from `count.tab` (cohort-level tab file in the sample output directory).

Sample IDs are auto-discovered as all immediate subdirectories of `results_root`. Use `cohort_id_pattern` (a regex string) in `config.R` to filter by name, or supply `cohort_ids` explicitly.

---

## Configuration

Edit `config.R` before running. Key settings:

```r
# Paths
results_root <- "/path/to/mitoHPC/results"   # parent directory of per-sample folders
report_dir   <- "/path/to/output/directory"  # where the report and QC TSV are written
pca_dir      <- file.path(report_dir, "1000G")

# Sample discovery
cohort_ids        <- NULL   # NULL = auto-discover all subdirs of results_root
cohort_id_pattern <- NULL   # optional regex filter, e.g. "^SAMPLE" or "^AA"
                             # supply cohort_ids explicitly to bypass discovery

# QC thresholds
min_median_depth  <- 100    # flag samples with median depth below this
contam_threshold  <- 0.02   # flag samples with haplocheck level >= this
min_hom_vaf       <- 0.95   # VAF threshold for homoplasmic classification

# Cohort label (appears in the report title)
cohort_label <- "My cohort"
```

---

## Running

Use the provided `render.R` wrapper — it ensures the correct `config.R` is sourced and writes the HTML to `report_dir` as set in your config.

```bash
cd /path/to/mtDNA_report/
Rscript render.R
```

Or from R:

```r
setwd("/path/to/mtDNA_report/")
source("render.R")
```

> **Important:** always run from the `mtDNA_report/` directory. If you call `rmarkdown::render()` directly from a different working directory, R will pick up the wrong `config.R` and write output to the wrong location.

Runtime is roughly 5–15 minutes for ~500 samples, depending on VCF parsing speed.

---

## Outputs

| File | Description |
|------|-------------|
| `mtDNA_report_<date>.html` | Self-contained interactive HTML report |
| `qc_summary.tsv` | Per-sample QC table (written to `report_dir`) |

The HTML report is fully self-contained (no external dependencies) and can be shared directly.

### Report sections

| Section | Content |
|---------|---------|
| Coverage | Median depth per sample, mtDNA copy number, per-position base coverage, cumulative coverage, CV distribution |
| Contamination | Haplocheck contamination level per sample, level vs depth scatter |
| Haplogroup | Macro-haplogroup distribution, assignment quality |
| PCA | Cohort projected onto 1000 Genomes PCA space (PC1/2, PC2/3, scree) |
| Variants | VAF distribution, consequence breakdown, MLC score, MSS, homoplasmy/heteroplasmy by clade, frequency spectrum, summary pies |
| Relatedness | Pairwise Jaccard similarity distribution, flagged duplicate/relative pairs |
| QC Summary | Interactive per-sample table with all QC metrics; exports to CSV/Excel |

### QC flags

A sample is flagged (`QC_pass = FALSE`) if any of the following apply:

- Median depth < `min_median_depth`
- Haplocheck contamination level ≥ `contam_threshold`
- Haplogroup not assigned

The `related` column tags samples with Jaccard similarity = 1.0 to another sample (probable duplicate or sample swap) but does not affect `QC_pass`.

---

## Variant filtering

PASS variants are retained after:

1. Removing positions in the blacklist (`variant_blacklist` in `config.R`)
2. Excluding variants with `strand_bias` in the FILTER field
3. Minimum VAF ≥ 0.03

Classification:
- **Homoplasmic**: VAF ≥ `min_hom_vaf` (default 0.95)
- **Heteroplasmic**: 0.03 ≤ VAF < `min_hom_vaf`
