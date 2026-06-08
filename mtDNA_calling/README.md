# MitoHPC : Mitochondrial High Performance Caller #

This tool has been installed at:
```
/stornext/Bioinf/data/lab_bahlo/software/apps/MitoHPC/
```
How it's installed: [install.sh](https://github.com/bahlolab/mtSwirl_HPC/blob/main/mitoHPC/install.sh)

## Citing ##

A bioinformatics pipeline for estimating mitochondrial DNA copy number and heteroplasmy levels from whole genome sequencing data, Battle et. al, NAR 2022
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9112767/ 

## PIPELINE USAGE ##

### SETUP ENVIRONMENT ###

```bash
# Move to your working directory
cd /vast/scratch/users/$USER/1000G/mitoHPC

# Define the path to the MitoHPC scripts
HP_SDIR=/stornext/Bioinf/data/lab_bahlo/software/apps/MitoHPC/scripts

# Copy the parameter initialization script (with prompt before overwriting)
cp -i $HP_SDIR/init.sh .

# View and edit the init file if needed (e.g., update paths)
cat init.sh
nano init.sh

# ⚠️ IMPORTANT: Update this line in init.sh if needed
# export HP_ADIR=/vast/scratch/users/wang.lo/1000G/data/WGS

# set no subsampling for more accurate heteroplasmy level
# export HP_L=

# Source the init file to load environment variables
. ./init.sh    # or: source ./init.sh

# Verify HP_ variables are correctly set
printenv | grep '^HP_' | sort

# View input file (if present)
nano "$HP_IN"

```

### RUN PIPELINE  ###

```bash
# Generate the command script from run.sh and save it as run.all.sh
$HP_SDIR/run.sh > run.all.sh

# Execute the generated command script
# Stdout and stderr will be logged in output.log
bash ./run.all.sh > output.log 2>&1
# /usr/bin/time -v bash ./run.all.sh > output.log 2>&1
```

### RUN PIPELINE IN PARALLEL ###

```bash
screen

# Generate the command script from run.sh and save it as run.all.sh
$HP_SDIR/run.sh > run.all.sh

module load parallel

# grep filter.sh ./run.all.sh | parallel && getSummary.sh
grep filter.sh ./run.all.sh | parallel
bash /stornext/Bioinf/data/lab_bahlo/software/apps/MitoHPC/scripts//getSummary.sh /vast/scratch/users/$USER/1000G/mitoHPC/out/

# This approach isn't scalable — we need to write scripts that submit jobs via SLURM.

```

### RE-RUN PIPELINE (optional) ###

```bash
cd /vast/scratch/users/$USER/1000G/mitoHPC

# Set the path to the MitoHPC script directory
HP_SDIR=/stornext/Bioinf/data/lab_bahlo/software/apps/MitoHPC/scripts

# Remove old input file if it exists
rm -f in.txt

# Edit parameters as needed (e.g., input/output directories, sample list)
nano init.sh

# Load updated configuration into the current shell session
. ./init.sh

# Regenerate the full command script based on current settings
$HP_SDIR/run.sh > run.all.sh

# Execute the pipeline and log all output
bash ./run.all.sh > output.log 2>&1

```

## PIPELINE USAGE (for multiple samples) ##

### Step 1: Set up input table ###
* Run `generate_input.R` to generate input file 

### Step 2: Run multiple samples ###
* Using the input file from Step 1, run `run_multisample.R` to generate a script file for each input sample
* Run combined shell file to run all samples

```bash
sbatch /path/to/combined/run_all.sh
```
