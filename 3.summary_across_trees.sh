#!/bin/sh

#SBATCH --job-name=summary
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2g
#SBATCH --mail-user=xueling.yi@senckenberg.de
#SBATCH --mail-type=END,FAIL

module unload R/4.0.3
module load R/4.4.1

rdata=sif179_100tree_batch1_diet6
id_start=1
id_end=10

for code in raw conserve
do
Rscript aggregate.R ${rdata} ${id_start} ${id_end} ${code} 
#Rscript aggregate_genus.R ${batch} ${id_start} ${id_end} ${code}
done
