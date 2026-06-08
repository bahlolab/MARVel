# MitoHPC : Mitochondrial High Performance Caller #

This tool has been installed at:
```
/stornext/Bioinf/data/lab_bahlo/software/apps/MitoHPC/
```
How it's installed: [install.sh](https://github.com/bahlolab/mtSwirl_HPC/blob/main/mitoHPC/install.sh)

## Citing ##

A bioinformatics pipeline for estimating mitochondrial DNA copy number and heteroplasmy levels from whole genome sequencing data, Battle et. al, NAR 2022
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9112767/ 

## PIPELINE USAGE (for multiple samples) ##

### Step 1: Set up input table ###
* Run `generate_input.R` to generate input file 

### Step 2: Run multiple samples ###
* Using the input file from Step 1, run `run_multisample.R` to generate a script file for each input sample
* Run combined shell file to run all samples

```bash
sbatch /path/to/combined/run_all.sh
```
