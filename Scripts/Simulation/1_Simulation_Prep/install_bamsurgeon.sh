# load conda
module load miniconda3
# create a new environment
conda create -y -p /vast/scratch/users/$USER/conda_envs/bamsurgeon python=3.10
#activate it
conda activate /vast/scratch/users/$USER/conda_envs/bamsurgeon
# install BamSurgeon
conda install -y -c bioconda -c conda-forge bamsurgeon
# check
addsnv.py --help
