### DOWNLOAD PIPELINE
git clone https://github.com/dpuiu/MitoHPC.git

### # SETUP ENVIRONMENT (important)

# Navigate to the script directory
cd MitoHPC/scripts

# Set HP_SDIR to current path
export HP_SDIR=$(pwd)

# Optional: persist this variable for future sessions
echo "export HP_SDIR=$(pwd)" >> ~/.bashrc

# Load environment variables and default settings
. ./init.sh
# Or choose a reference-specific init:
# . ./init.hs38DH.sh     # for human GRCh38
# . ./init.hg19.sh       # for human hg19
# . ./init.mm39.sh       # for mouse mm39

### INSTALL PREREQUISITES & CHECK INSTALL
# Install tools (bwa, samtools, bedtools, etc.) if not already available
$HP_SDIR/install_prerequisites.sh

# Optionally, force re-install latest versions:
# $HP_SDIR/install_prerequisites.sh -f

# Run installation check script, if successfull => "Success message!"
$HP_SDIR/checkInstall.sh

# View log of tool paths and versions
cat checkInstall.log
