#!/usr/bin/env Rscript
library(ape)
## run this in sif179
args <- commandArgs(TRUE)
batch = as.numeric(args[1])
rdata = args[2]
id_start = as.numeric(args[3])
id_end = as.numeric(args[4])
diet.code=args[5]

trees = read.table(paste0("tree100_batch", batch, ".txt"))
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
