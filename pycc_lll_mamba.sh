#!/bin/bash
#SBATCH --job-name=pycc_lll_mamba
#SBATCH --output=pycc_lll_mamba_%j.log
#SBATCH -p interactive
#SBATCH --ntasks=1
#SBATCH -c 4
#SBATCH --mem=60000
#SBATCH --time=16:00:00

eval "$(mamba shell hook --shell bash)"
mamba activate pybioenv2

# Works with the old spack python 
#/ref/rmlab/software/spack/opt/spack/linux-rocky8-x86_64/gcc-8.5.0/miniconda3-4.10.3-btkingmhzajcvwsakbdkqvl23c3svufg/bin/conda
python -c "import sys; print(sys.executable)"

#export PATH=/home/l.llaci/.conda/envs/myrenv/bin/jupyter:$PATH

export XDG_RUNTIME_DIR=${TMPDIR}
 
port=$(shuf -i 1026-9998 -n 1)
 
echo -e "\nStarting Jupyter Lab on port ${port} on the $(hostname) server."
echo -e "\nSSH tunnel command: ssh -NL ${port}:$(hostname):${port} ${USER}@login.htcf.wustl.edu"
echo -e "\nLocal URI: http://localhost:${port}"

jupyter lab --no-browser --port=${port} --ip='*' #--NotebookApp.token='' --NotebookApp.password=''
 
date


echo -e "\nSSH tunnel command: ssh -NL ${port}:$(hostname):${port} ${USER}@login.htcf.wustl.edu"

set -x
