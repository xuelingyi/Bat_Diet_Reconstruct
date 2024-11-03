#!/bin/sh

#SBATCH --job-name=summary
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8g
#SBATCH --array=1-10
#SBATCH --mail-user=xueling.yi@senckenberg.de
#SBATCH --mail-type=END,FAIL

# array index is used as the tree ID
## needed larger mem because estimates from all 10 chains will be stored
module unload R/4.0.3
module load R/4.4.1

batch=1
tree_id=$(sed -n ${SLURM_ARRAY_TASK_ID}p tree100_batch${batch}.txt)
rdata=sif179_100tree_batch1_diet6

for code in raw conserve
do
cd trees_${code}/tree_${tree_id}
Rscript ../../vertlife_extract_threshold_probability.R ${batch} ${SLURM_ARRAY_TASK_ID} ${rdata}
cd ../../
done
