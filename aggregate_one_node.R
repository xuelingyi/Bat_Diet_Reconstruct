#!/usr/bin/env Rscript
library(ape)

args <- commandArgs(TRUE)
batch = as.numeric(args[1])
rdata = args[2]
id_start = as.numeric(args[3])
id_end = as.numeric(args[4])
diet.code=args[5]
node_name = args[6]

trees = read.table(paste0("tree100_batch", batch, ".txt"))

load(paste0(rdata, ".RData"))
## tree_100_sub, sif, diet.in, diet.states
my.tips = read.table(paste0(node_name, "_tips"))

### check monophyly first
cat("check monophyly...")
cat("\n")
for(i in tree_100_sub){
  if(!is.monophyletic(i, my.tips$V1)){
    print("not monophyly!")
    break()
  }
}

setwd(paste0("trees_", diet.code))


## get the first tree
id = id_start
tree = tree_100_sub[[id]]
n.edge = length(tree$edge)/2
tree_id = trees[id, "V1"]
## get the raw estimates of the focal node in the first tree
my.node = tree$node.label[getMRCA(tree, my.tips$V1)-length(tree$tip.label)]
load(paste0("tree_", trees[id, "V1"], "/estimates.RData"))
for (d in diet.in){
  data = as.data.frame(get(paste0(d, "_estimates"))[[my.node]])
  data$Taxon = node_name
  data$tree = tree_id
  assign(paste0(d, "_summary_aggregate"), data)
}


## get the remaining trees
for(id in (id_start+1):id_end){
  tree = tree_100_sub[[id]]
  n.edge = length(tree$edge)/2
  tree_id = trees[id, "V1"]
  my.node = tree$node.label[getMRCA(tree, my.tips$V1)-length(tree$tip.label)]
  ## get the raw estimates of the focal node in the rest of trees
  load(paste0("tree_", trees[id, "V1"], "/estimates.RData"))
  for (d in diet.in){
    summary = get(paste0(d, "_summary_aggregate"))
    data = as.data.frame(get(paste0(d, "_estimates"))[[my.node]])
    data$Taxon = node_name
    data$tree = tree_id
    summary = rbind(summary, data)
    assign(paste0(d, "_summary_aggregate"), summary)
  }
}

setwd("../")

save(list = grep("aggregate", ls(), value=T),
     file=paste0(node_name, "_mcmc_aggregate_diet", length(diet.in), "_", diet.code, "_", id_start, "_", id_end, ".RData"))
