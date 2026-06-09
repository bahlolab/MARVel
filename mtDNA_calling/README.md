# MitoHPC: Mitochondrial High Performance Caller

Scripts for running [MitoHPC](https://github.com/dpuiu/MitoHPC) — a pipeline for estimating mitochondrial DNA copy number and heteroplasmy from whole-genome sequencing data.

## Citation

Battle et al. (2022). A bioinformatics pipeline for estimating mitochondrial DNA copy number and heteroplasmy levels from whole genome sequencing data. *Nucleic Acids Research*. https://doi.org/10.1093/nar/gkac290

---

## Installation

See [install.sh](install.sh) for a minimal setup. In brief:

```bash
git clone https://github.com/dpuiu/MitoHPC.git
cd MitoHPC/scripts
export HP_SDIR=$(pwd)
. ./init.sh   # or init.hs38DH.sh / init.hg19.sh / init.mm39.sh
$HP_SDIR/install_prerequisites.sh
$HP_SDIR/checkInstall.sh
```

---

## Single-sample usage

### 1. Set up working directory

```bash
mkdir -p /path/to/workdir && cd /path/to/workdir

HP_SDIR=/path/to/MitoHPC/scripts
cp $HP_SDIR/init.sh .
```

Edit `init.sh` to set `HP_ADIR` (directory containing your BAM) and any other parameters (e.g. set `HP_L=` to disable subsampling for more accurate heteroplasmy estimation).

```bash
. ./init.sh
printenv | grep '^HP_' | sort   # verify
```

### 2. Create input file

`in.txt` is a tab-separated file with three columns (no header): sample name, BAM path, output prefix.

```
SAMPLE001	/path/to/SAMPLE001.bam	out/SAMPLE001/SAMPLE001
```

### 3. Run

```bash
$HP_SDIR/run.sh > run.all.sh
bash ./run.all.sh > output.log 2>&1
```

---

## Multi-sample usage (SLURM)

### Step 1: Generate input table

Edit the paths in [generate_input_list.R](generate_input_list.R) then run:

```bash
Rscript generate_input_list.R
```

This scans BAM directories, checks for missing indices, warns on duplicate sample names, and writes `input_list_complete.txt`.

### Step 2: Generate per-sample SLURM scripts

Edit paths and SBATCH settings in [run_multisample.R](run_multisample.R) then run:

```bash
Rscript run_multisample.R
```

This writes one `.sh` script per sample under `sh_scripts/`, plus a combined `sh_scripts/run_all.sh`.

### Step 3: Submit

```bash
bash sh_scripts/run_all.sh
```

### Step 4: Check for incomplete jobs (optional)

After jobs finish, run [check_jobs.R](check_jobs.R) to identify any samples that did not complete:

```bash
Rscript check_jobs.R
```

This checks for expected mitoHPC output files per sample and writes `rerun_failed.sh` listing any incomplete jobs. Submit with:

```bash
bash rerun_failed.sh
```

---

## Key parameters in `init.sh`

| Variable  | Description                                              |
|-----------|----------------------------------------------------------|
| `HP_ADIR` | Directory containing input BAM files                     |
| `HP_IN`   | Path to input file (`in.txt`)                            |
| `HP_L`    | Subsampling depth (leave blank for no subsampling)       |
| `HP_M`    | Variant caller (`gatk` default; `mutect2` also supported)|
