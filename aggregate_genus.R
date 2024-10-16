#!/usr/bin/env Rscript
library(ape)
## run this in sif179
args <- commandArgs(TRUE)
batch = as.numeric(args[1])
id_start = as.numeric(args[2])
id_end = as.numeric(args[3])
diet.code=args[4]

trees = read.table(paste0("tree100_batch", batch, ".txt"))
load(paste0("sif179_100tree_batch", batch, "_diet6.RData"))
## tree_100_sub, sif, diet.in, diet.states

setwd(paste0("trees_", diet.code))

## get the first tree
id = id_start
tree_id = trees[id, "V1"]
tree = tree_100_sub[[id]]
load(paste0("tree_", tree_id, "/estimates.RData"))
not_mono = NULL
for (d in diet.in){
  summary = NULL
  for(g in unique(sif$Genus)){
    if(length(sif[sif$Genus ==g, "vertlife_name"]) > 1){
      ## only using the genera having more than one species
      if(is.monophyletic(tree, sif[sif$Genus == g, "vertlife_name"])){
        n = tree$node.label[(getMRCA(tree, sif[sif$Genus == g, "vertlife_name"]) - 179)]
        data = as.data.frame(get(paste0(d, "_estimates"))[[n]])
        data$Taxon = n
        data$tree = tree_id
        summary = rbind(summary, data)
      } else {
        ## only do this once per tree
        if(d == "Arthropods"){not_mono = rbind(not_mono, c(tree_id, g))}
      }
    }
  }
  assign(paste0(d, "_summary_aggregate"), summary)
}

## get the remaining trees
for(id in (id_start+1):id_end){
  tree_id = trees[id, "V1"]
  tree = tree_100_sub[[id]]
  load(paste0("tree_", tree_id, "/estimates.RData"))
  
  for(d in diet.in){
    summary = NULL
    for(g in unique(sif$Genus)){
      if(length(sif[sif$Genus ==g, "vertlife_name"]) > 1){
        ## only using the genera having more than one species
        if(is.monophyletic(tree, sif[sif$Genus == g, "vertlife_name"])){
          n = tree$node.label[(getMRCA(tree, sif[sif$Genus == g, "vertlife_name"]) - 179)]
          data = as.data.frame(get(paste0(d, "_estimates"))[[n]])
          data$Taxon = n
          data$tree = tree_id
          summary = rbind(summary, data)
        } else {
          ## only do this once per tree
          if(d == "Arthropods"){not_mono = rbind(not_mono, c(tree_id, g))}
        }
      }
    }
    
    data = get(paste0(d, "_summary_aggregate"))
    data = rbind(data, summary)
    assign(paste0(d, "_summary_aggregate"), data)
  }
}

setwd("../")

not_mono = as.data.frame(not_mono)
names(not_mono) = c("tree", "Genus")

save(list = c("not_mono", grep("aggregate", ls(), value=T)), 
     file=paste0("mcmc_aggregate_diet6_", diet.code, "_GENUS_", id_start, "_", id_end, ".RData"))

