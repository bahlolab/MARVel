# Projection of External Cohorts into mtDNA PCA Reference Space

## Overview

This repository provides a standardized workflow for projecting external whole-genome sequencing (WGS) cohorts into a pre-computed mitochondrial DNA (mtDNA) principal component analysis (PCA) reference space derived from the 1000 Genomes Project.

External cohorts can be embedded into this reference space through a harmonized workflow to ensure consistency and reproducibility across independent datasets.

---

## Workflow

The projection procedure consists of five key steps:

### 1. Restrict to reference variant set
Subset your dataset to the same set of mtDNA SNPs used in the reference PCA: mtPCA.variant_ids.txt

### 2. Harmonize genotype encoding
Ensure genotype encodings match the reference (e.g., allele order and coding scheme).

### 3. Impute missing genotypes
Replace missing genotypes with the reference allele to ensure compatibility with the reference matrix.

### 4. Center genotype matrix
Center genotype vectors using the pre-computed reference means: mtPCA.mu.txt

### 5. Project into PCA space
Project samples into the reference PCA space by multiplying the centered genotype matrix by the pre-computed loadings: mtPCA.loadings.PC1-20.txt


---

## Output

- Principal component (PC) scores for each sample (PC1–PC20)  
- Coordinates directly comparable to the 1000 Genomes mtDNA reference space  

This approach ensures reproducible embedding across independent datasets.

---

## Important Considerations

### Variant overlap requirement

We recommend using the pre-computed PCA loadings only when variant overlap between the 1000 Genomes reference and the external dataset is high:

- **>95% overlap**: Use pre-computed loadings  
- **≤95% overlap**: Recompute PCA loadings  

---

## Alternative Workflow (Low Variant Overlap)

For datasets with ≤95% variant overlap:

1. Identify the intersection of high-quality variants between:
   - the external dataset  
   - the 1000 Genomes reference  
2. Use the shared variant set (n = 317)  
3. Recompute PCA loadings using this subset  

---

## Resources Provided

- 1000 Genomes mtDNA reference dataset  
- Pre-computed PCA parameters:
  - `mtPCA.variant_ids.txt`
  - `mtPCA.mu.txt`
  - `mtPCA.loadings.PC1-20.txt`
- R scripts implementing:
  - Projection using pre-computed loadings  
  - Re-computation of PCA loadings  

---

## Notes

- This workflow is designed for short-read WGS mtDNA data  
- Ensure consistent QC and variant filtering prior to projection  
- Interpretation of PCs should consider mtDNA haplogroup structure  
