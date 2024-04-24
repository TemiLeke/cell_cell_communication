#!/bin/bash
#SBATCH --job-name=topic_modelling        # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --mem-per-cpu=16G
#SBATCH --ntasks=1
#SBATCH --time=96:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email on job start, end and fault
#SBATCH --mail-user=tadeoye@usf.edu

module load python/3.8.5
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip

pip install --no-index -r /work/t/tadeoye/scRNAseq_AD_meta_analysis/scripts/requirements.txt
python topic_model.py