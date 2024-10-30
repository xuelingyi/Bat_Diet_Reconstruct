#!/bin/sh

#SBATCH --job-name=mcmc
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10g
#SBATCH --array=1-10
#SBATCH --mail-user=xueling.yi@senckenberg.de
#SBATCH --mail-type=END,FAIL

# array index is used as the chain ID
module unload R/4.0.3
module load R/4.4.1

batch=arg1
index=arg2
nitt=2000000
thin=50
fix=1
diet_coding=arg7
## raw, conserve
rdata=arg8

Rscript ../../vertlife_MCMC_diet6.R ${batch} ${index} ${SLURM_ARRAY_TASK_ID} ${nitt} ${thin} ${fix} ${diet_coding} ${rdata}

