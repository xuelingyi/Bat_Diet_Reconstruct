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
rdata = args[3]
working_dir = "./"

load(paste0("../../", rdata, ".RData"))
# tree_100_sub, sif, diet.in, diet.states
tree = tree_100_sub[[batch_tree_index]]
n.edge = length(tree$edge)/2
rm(tree_100_sub)

trees = read.table(paste0("../../tree100_batch", batch, ".txt"))
tree_id = trees[batch_tree_index, "V1"]

# This function adds the posterior samples of each chain that support each diet state per taxon to the current sample counts. To determine the diet state that each sample supports, the function checks where the liability lies relatively to the threshold values.
modify_diet_counts <- function(diet, estimates, ...)
{
  
  # Add the samples for the root node.
  
  ## state 0 is always included, so threshold 0 is always included
  # If the liability is below 0, this sample supports no such diet
  estimates[['root']]$Absent <- estimates[['root']]$Absent + length(
    which(fit$Sol[, paste0("trait",d)] < 0)
  )
  
  ## binary diets
  if(length(unique(sif[, d])) == 2){
    state = sort(unique(sif[, d]))[2] 
    estimates[['root']][names(diet.states[diet.states==state])][[1]] <- estimates[['root']][names(diet.states[diet.states==state])][[1]] + length(
      which(fit$Sol[, paste0("trait",d)] >= 0))
  }
  
  ## one non-zero cutpoint 
  if(length(unique(sif[, d])) == 3){
    state1 = sort(unique(sif[, d]))[2] 
    estimates[['root']][names(diet.states[diet.states==state1])][[1]] <- estimates[['root']][names(diet.states[diet.states==state1])][[1]] + length(
      which(fit$Sol[, paste0("trait",d)] >= 0 & 
              fit$Sol[, paste0("trait",d)] < fit$CP[,paste0("cutpoint.trait",d, ".1")])
    )
    ## note: Sol and CP have the same nrow corresponding to the sample size
    
    state2 = sort(unique(sif[, d]))[3] 
    estimates[['root']][names(diet.states[diet.states==state2])][[1]] <- estimates[['root']][names(diet.states[diet.states==state2])][[1]] + length(
      which(fit$Sol[, paste0("trait",d)] >= fit$CP[,paste0("cutpoint.trait",d, ".1")] ))
  }
  
  ## two non-zero cutpoints 
  if(length(unique(sif[, d])) == 4){
    state1 = sort(unique(sif[, d]))[2] 
    estimates[['root']][names(diet.states[diet.states==state1])][[1]] <- estimates[['root']][names(diet.states[diet.states==state1])][[1]] + length(
      which(fit$Sol[, paste0("trait",d)] >= 0 & 
              fit$Sol[, paste0("trait",d)] < fit$CP[,paste0("cutpoint.trait",d, ".1")])
    )
    
    state2 = sort(unique(sif[, d]))[3] 
    estimates[['root']][names(diet.states[diet.states==state2])][[1]] <- estimates[['root']][names(diet.states[diet.states==state2])][[1]] + length(
      which(fit$Sol[, paste0("trait",d)] >= fit$CP[,paste0("cutpoint.trait",d, ".1")] &
              fit$Sol[, paste0("trait",d)] < fit$CP[,paste0("cutpoint.trait",d, ".2")]))
    
    state3 = sort(unique(sif[, d]))[4] 
    estimates[['root']][names(diet.states[diet.states==state3])][[1]] <- estimates[['root']][names(diet.states[diet.states==state3])][[1]] + length(
      which(fit$Sol[, paste0("trait",d)] >= fit$CP[,paste0("cutpoint.trait",d, ".2")]))
  }
  
  # root is done. Do the same thing for all remaining taxa.
  for(x in 2:length(sp_names)){
    estimates[[sp_names[x]]]$Absent <- estimates[[sp_names[x]]]$Absent + length(
      which((fit$Sol[,paste0("trait",d)] + #this is the root value
               fit$Sol[, paste("trait",d, ".vertlife_name.", sp_names[x], sep = '')]) < 0))
    
    if(length(unique(sif[, d])) == 2){
      state = sort(unique(sif[, d]))[2] 
      estimates[[sp_names[x]]][names(diet.states[diet.states==state])][[1]] <- estimates[[sp_names[x]]][names(diet.states[diet.states==state])][[1]] + length(
        which((fit$Sol[,paste0("trait",d)] + #this is the root value
                 fit$Sol[, paste("trait",d, ".vertlife_name.", sp_names[x], sep = '')]) >= 0))
    }
    
    ## one non-zero cutpoint 
    if(length(unique(sif[, d])) == 3){
      state1 = sort(unique(sif[, d]))[2] 
      estimates[[sp_names[x]]][names(diet.states[diet.states==state1])][[1]] <- estimates[[sp_names[x]]][names(diet.states[diet.states==state1])][[1]] + length(
        which((fit$Sol[,paste0("trait",d)] + fit$Sol[, paste("trait",d, ".vertlife_name.", sp_names[x], sep = '')]) >= 0 & 
                (fit$Sol[,paste0("trait",d)] + fit$Sol[, paste("trait",d, ".vertlife_name.", sp_names[x], sep = '')]) < fit$CP[,paste0("cutpoint.trait",d, ".1")])
      )
      
      state2 = sort(unique(sif[, d]))[3] 
      estimates[[sp_names[x]]][names(diet.states[diet.states==state2])][[1]] <- estimates[[sp_names[x]]][names(diet.states[diet.states==state2])][[1]] + length(
        which((fit$Sol[,paste0("trait",d)] + fit$Sol[, paste("trait",d, ".vertlife_name.", sp_names[x], sep = '')]) >= fit$CP[,paste0("cutpoint.trait",d, ".1")] ))
    }
    
    ## two non-zero cutpoints 
    if(length(unique(sif[, d])) == 4){
      state1 = sort(unique(sif[, d]))[2] 
      estimates[[sp_names[x]]][names(diet.states[diet.states==state1])][[1]] <- estimates[[sp_names[x]]][names(diet.states[diet.states==state1])][[1]] + length(
        which((fit$Sol[,paste0("trait",d)] + fit$Sol[, paste("trait",d, ".vertlife_name.", sp_names[x], sep = '')]) >= 0 & 
                (fit$Sol[,paste0("trait",d)] + fit$Sol[, paste("trait",d, ".vertlife_name.", sp_names[x], sep = '')]) < fit$CP[,paste0("cutpoint.trait",d, ".1")])
      )
      
      state2 = sort(unique(sif[, d]))[3] 
      estimates[[sp_names[x]]][names(diet.states[diet.states==state2])][[1]] <- estimates[[sp_names[x]]][names(diet.states[diet.states==state2])][[1]] + length(
        which((fit$Sol[,paste0("trait",d)] + fit$Sol[, paste("trait",d, ".vertlife_name.", sp_names[x], sep = '')]) >= fit$CP[,paste0("cutpoint.trait",d, ".1")] &
                (fit$Sol[,paste0("trait",d)] + fit$Sol[, paste("trait",d, ".vertlife_name.", sp_names[x], sep = '')]) < fit$CP[,paste0("cutpoint.trait",d, ".2")]))
      
      state3 = sort(unique(sif[, d]))[4] 
      estimates[[sp_names[x]]][names(diet.states[diet.states==state3])][[1]] <- estimates[[sp_names[x]]][names(diet.states[diet.states==state3])][[1]] + length(
        which((fit$Sol[,paste0("trait",d)] + fit$Sol[, paste("trait",d, ".vertlife_name.", sp_names[x], sep = '')]) >= fit$CP[,paste0("cutpoint.trait",d, ".2")]))
    }
  }
  return(estimates)
}


# Look for model fit files in the working directory.
Rda_files_in_working_dir <- list.files(
  path = working_dir, pattern = '^\\d+\\.Rda$'
)

# Load the first .Rda file
cat('Now loading ', Rda_files_in_working_dir[1], ' ...\n', sep ='')
load(paste(working_dir, Rda_files_in_working_dir[1], sep = ''))
samples_per_chain <- nrow(fit$Sol)
### get the number of samples in the chain, and store the samples of thresholds in a vector.
thresholds = gsub("cutpoint.trait", "", colnames(fit$CP))
for(t in thresholds){
  data = rep(NA, samples_per_chain * length(Rda_files_in_working_dir))
  data[1:samples_per_chain] <- fit$CP[,which(thresholds == t)]
  assign(t, data)
}
# Get the names of extant species and internal nodes.
sp_names <- gsub(
  'traitArthropods.vertlife_name.', '', 
  grep('Arthropods', colnames(fit$Sol), value = TRUE)
)
sp_names[sp_names == 'traitArthropods'] <- 'root'
# Initialise a list to store the number of samples in support of each diet state per taxon per diet type.
for(d in diet.in){
  estimates <- list()
  for ( sp in sp_names ){
    estimates[[sp]] <- list(Absent = 0, Complementary = 0, Predominant = 0, Strict = 0)[(0:3) %in% unique(sif[,d])]
  }
  # Store the sample counts from the first chain. For each species & node, it is in state s in how many iterations. 
  assign(paste0(d, "_estimates"), modify_diet_counts(diet=d, estimates))
}

# Load all remaining .Rda files, extract their threshold samples, and append them to the vector.
for ( i in 2:length(Rda_files_in_working_dir) ){
  rm(fit)
  gc()
  
  cat('Now loading ', Rda_files_in_working_dir[i], ' ...\n', sep ='')
  load(paste(working_dir, Rda_files_in_working_dir[i], sep = ''))
  
  ## append the thresholds
  for(t in thresholds){
    data = get(t)
    data[(((i - 1) * samples_per_chain) + 1):(i * samples_per_chain)] <- fit$CP[, which(thresholds == t)]
    assign(t, data)
  }
  
  ## add counts of each diet state for probability calculations
  for(d in diet.in){
    assign(paste0(d, "_estimates"), 
           modify_diet_counts(diet=d, estimates=get(paste0(d, "_estimates"))))
  }
}

save(list=grep("_estimates", ls(), value=T), file="estimates.RData")

# Write the median and the 95% HPD interval of each threshold to an output file.
## perhaps not necessary
my.thresholds = NULL
for(t in thresholds){
  my.median = median(get(t))
  my.low = HPDinterval(mcmc(get(t)))[1]
  my.high = HPDinterval(mcmc(get(t)))[2]
  
  my.thresholds = rbind(my.thresholds, 
                        c(t, my.median, my.low, my.high))
}
my.thresholds = as.data.frame(my.thresholds)
names(my.thresholds) = c("threshold", "median", "95HPD_low", "95HPD_high")
write.csv(my.thresholds, paste0(working_dir, "threshold.csv"), row.names = F, quote = F)
#system(paste0("open ", working_dir, "threshold.csv"))

# write out the diet probabilities per diet type per taxon.
# if chains have the same lengths, you should combine those for the focal taxa using estimates 
results.all = data.frame(sp_names)
names(results.all) = "Taxon"
for(d in diet.in){
  states = sort(unique(sif[, d]))
  results = data.frame(matrix(nrow=length(sp_names), ncol=length(states)+1))
  names(results) = c("Taxon", paste0(d, "_Prob_", 
                                     names(diet.states[diet.states %in% states])))
  results$Taxon = sp_names
  
  # Calculate the probabilities
  estimates = get(paste0(d, "_estimates"))
  for ( i in 1:nrow(results) ){
    for(s in states){
      results[i, paste0(d, "_Prob_", names(diet.states[diet.states == s]))] = 
        estimates[[sp_names[i]]][names(diet.states[diet.states==s])][[1]] / sum(
          unlist(estimates[[sp_names[i]]]))
      # sum should be the same: the total number of iterations across chains
    }
  }
  results.all = merge(results.all, results, by="Taxon")
}
write.csv(results.all, paste0(working_dir, "node_prob.csv"), row.names = F, quote = F)


## get the raw estimates of each focal node and tips
## tree_100_sub, sif, diet.in, diet.states
for (d in diet.in){
  summary = NULL
  for(n in c(mono.nodes, sif$vertlife_name[1:179])){
    data = as.data.frame(get(paste0(d, "_estimates"))[[n]])
    data$Taxon = n
    data$tree = tree_id
    summary = rbind(summary, data)
  }
  assign(paste0(d, "_summary"), summary)
}
save(list=grep("_summary", ls(), value=T),file="summary.RData")

