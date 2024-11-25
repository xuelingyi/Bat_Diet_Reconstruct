#!/usr/bin/env Rscript
library(ape)

########################################################################
## This script aggregates results across trees using the estimates saved in summary.RData (i.e., only tips and the labeled nodes of interest). 
## run this script in the same dir as the folders trees_conserve or trees_raw
# batch=1
# rdata=sif176_100tree_batch1_diet5
# id_start=1
# id_end=100

# for code in raw conserve
# do
# Rscript aggregate.R ${batch} ${rdata} ${id_start} ${id_end} ${code}
# done
#######################################################################


args <- commandArgs(TRUE)
batch = as.numeric(args[1])
rdata = args[2]
id_start = as.numeric(args[3])
id_end = as.numeric(args[4])
diet.code=args[5]

trees = read.table(paste0("tree100_batch", batch, ".txt"))
## or tree23_batch2 

load(paste0(rdata, ".RData"))
## tree_100_sub, sif, diet.in, diet.states

setwd(paste0("trees_", diet.code))

## get the first tree
id = id_start
#tree = tree_100_sub[[id]]
load(paste0("tree_", trees[id, "V1"], "/summary.RData"))
# diet_summary
for(d in diet.in){
  assign(paste0(d, "_summary_aggregate"), get(paste0(d, "_summary")))
}

## get the remaining trees
for(id in (id_start+1):id_end){
  #tree = tree_100_sub[[id]]
  load(paste0("tree_", trees[id, "V1"], "/summary.RData"))
  # diet_summary
  for(d in diet.in){
    data = get(paste0(d, "_summary_aggregate"))
    data = rbind(data, get(paste0(d, "_summary")))
    assign(paste0(d, "_summary_aggregate"), data)
  }
}

setwd("../")

save(list = grep("aggregate", ls(), value=T), 
     file=paste0("mcmc_aggregate_diet6_", diet.code, "_", id_start, "_", id_end, ".RData"))
