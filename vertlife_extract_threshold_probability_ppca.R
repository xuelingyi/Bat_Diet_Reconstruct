#!/usr/bin/env Rscript

## run this in the working dir that contains results 
# This script loads MCMCglmm chains and extracts the liability and threshold of each iteration to estimate the trait states of each iteration. Then it estimates the probabilities of each trait states at each node across chains. 

library(MCMCglmm, lib="~/R_4.4")
library(phytools, lib="~/R_4.4")
library(janitor, lib="~/R_4.4")

## NOTE: run this script within the dir vertlife/sif179/trees_${code}/tree_id

### load input files 
args <- commandArgs(TRUE)
batch=as.numeric(args[1])
batch_tree_index=as.numeric(args[2])
working_dir = "./"

load(paste0("../../inputs/sif176_ppca_topo11_index_", batch_tree_index, ".RData"))
# tree, scores, diet.in, diet.states, mono.nodes
n.edge = length(tree$edge)/2

trees = read.table(paste0("../../tree100_batch", batch, ".txt"))
tree_id = trees[batch_tree_index, "V1"]

# Look for model fit files in the working directory.
Rda_files_in_working_dir <- list.files(
  path = working_dir, pattern = '^\\d+\\.Rda$'
)

## total iterations
iterations = 0

# Load the first .Rda file
cat('Now loading ', Rda_files_in_working_dir[1], ' ...\n', sep ='')
load(paste(working_dir, Rda_files_in_working_dir[1], sep = ''))
iterations = sum(iterations, nrow(fit$Sol))

# Get the names of extant species and internal nodes.
sp_names <- gsub(
  'traitPC1.vertlife_name.', '', 
  grep('PC1', colnames(fit$Sol), value = TRUE)
)
sp_names[sp_names == 'traitPC1'] <- 'root'

estimates = as.data.frame(matrix(nrow = length(sp_names), ncol=2))
names(estimates) = c("PC1", "PC2")
row.names(estimates) = sp_names

## results of chain 1
# Add the samples for the root node.
estimates['root',] = c(sum(fit$Sol[, "traitPC1"]), sum(fit$Sol[, "traitPC2"]))
# for all remaining taxa; mean across iterations of the first chain
for(x in 2:length(sp_names)){
  estimates[sp_names[x], "PC1"] <- sum(fit$Sol[, "traitPC1"] + fit$Sol[, paste0("traitPC1.vertlife_name.", sp_names[x])])
  estimates[sp_names[x], "PC2"] <- sum(fit$Sol[, "traitPC2"] + fit$Sol[, paste0("traitPC2.vertlife_name.", sp_names[x])])
}

# Load all remaining .Rda files and summarize
for ( i in 2:length(Rda_files_in_working_dir) ){
  rm(fit)
  gc()
  
  cat('Now loading ', Rda_files_in_working_dir[i], ' ...\n', sep ='')
  load(paste(working_dir, Rda_files_in_working_dir[i], sep = ''))
  iterations = sum(iterations, nrow(fit$Sol))
  
  estimates['root', "PC1"] = sum(estimates['root', "PC1"] + sum(fit$Sol[, "traitPC1"])) 
  estimates['root', "PC2"] = sum(estimates['root', "PC2"] + sum(fit$Sol[, "traitPC2"]))
  
  # for all remaining taxa.
  for(x in 2:length(sp_names)){
    estimates[sp_names[x], "PC1"] <- sum(estimates[sp_names[x], "PC1"], 
                                         sum(fit$Sol[, "traitPC1"] + fit$Sol[, paste0("traitPC1.vertlife_name.", sp_names[x])]))
    estimates[sp_names[x], "PC2"] <- sum(estimates[sp_names[x], "PC2"],
                                         sum(fit$Sol[, "traitPC2"] + fit$Sol[, paste0("traitPC2.vertlife_name.", sp_names[x])]))
  }
}

if(iterations != 360000){print("wrong numbers of iterations!!")}

estimates$PC1 = estimates$PC1 / iterations
estimates$PC2 = estimates$PC2 / iterations

estimates$vertlife_name = row.names(estimates)
estimates$type = "output"
scores$type = "input"
estimates = rbind(estimates, scores[, c("PC1", "PC2", "vertlife_name", "type")])
estimates$tree_id = paste0("tree", tree_id)

write.csv(estimates, paste0(working_dir, "estimates.csv"), row.names = F, quote = F)
