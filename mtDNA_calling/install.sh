#!/usr/bin/env bash
set -euo pipefail

# Clone MitoHPC
git clone https://github.com/dpuiu/MitoHPC.git
cd MitoHPC/scripts

export HP_SDIR=$(pwd)

# Source default init (edit to use genome-specific init if needed)
# Options: init.hs38DH.sh (GRCh38), init.hg19.sh, init.mm39.sh
. ./init.sh

# Install prerequisites (add -f to force reinstall)
$HP_SDIR/install_prerequisites.sh

# Verify installation
$HP_SDIR/checkInstall.sh
cat checkInstall.log
