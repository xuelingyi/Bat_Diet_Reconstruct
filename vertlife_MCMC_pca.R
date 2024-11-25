#!/usr/bin/env Rscript

library(ape)
library(MCMCglmm, lib="~/R_4.4")
library(phytools, lib="~/R_4.4")

########################################################################
#### In the dir of the input folder and the R scripts, create a folder for a specified diet_coding (e.g., raw; conserve if basal-insect coding)
#### in the diet_coding folder, create a folder for each tree to run using the name of the tree index in tree23_batch2.txt
#### in the folder of each tree, run the following script on a command line: 

# batch=1
# index=1  ## loop through 1-23, trees in the index file
# chain=1  ## loop through 1-10 chains
# nitt=2000000
# thin=50
# diet_coding=conserve  
### Rscript ../../vertlife_MCMC_pca.R ${batch} ${index} ${chain} ${nitt} ${thin} ${diet_coding}

########################################################################

# Read the user-provided command line arguments.
args <- commandArgs(TRUE)
batch=as.numeric(args[1])
batch_tree_index=as.numeric(args[2])
chain_id <- as.numeric(args[3])
my.nitt <- as.numeric(args[4])
# 2000000 default
my.thin <- as.numeric(args[5])
#50 default
diet.coding <- args[7]
## raw, conserve

out.dir = "./"
my.burnin = my.nitt*0.1
# 10% burin

load(paste0("../../inputs/sif176_ppca_topo11_index_", batch_tree_index, ".RData"))
# tree, scores, diet.in, diet.states, mono.nodes
n.edge = length(tree$edge)/2

## priors
my.prior = list(
  ## prior for the random effect
  G = list(
    G1 = list(
      V=diag(length(diet.in)), nu=length(diet.in),
      alpha.mu=rep(0, length(diet.in)), alpha.V=diag(rep(1000, length(diet.in)))
    )
  ),
  ## prior for the residual
  R = list(V = diag(length(diet.in)), nu=length(diet.in))
)

## use scale only when the tree is ultrametric 
if(is.ultrametric(tree)){
  my.scale=TRUE
  print("scale trees")
  cat("\n")
} else {
  my.scale=FALSE
  print("not ultrametric tree; do not scale")
  cat("\n")
}

# Set the seed and fit the model.
set.seed(chain_id)

## use ppca scores pc1 and pc2, continuous trait
# fit a model with MCMCglmm to estimate the correlation structure between pairs of diet items, accounting for phylogeny.
fit <- MCMCglmm(
  # Define the response variables, apply any needed transformations, and specify a distinct intercept per response.
  cbind(PC1, PC2) ~ trait - 1,
  
  # Specify a phylogenetic random effect on each intercept.
  # us: estimate both variance and covariance
  random =~ us(trait):vertlife_name,
  
  # Set the distribution for each response variable (gaussian or threshold)
  family = rep("gaussian", length(diet.in)),
  
  # only scale is ultrametric tree
  scale=my.scale,
  
  # Integrate the phylogenetic variance-covariance matrix into the model; scale=F as the tree is not ultrametric
  ginverse = list(vertlife_name=inverseA(tree, nodes = 'ALL', scale = my.scale)$Ainv),
  
  # Set relatively uninformative priors.
  prior = my.prior,
  
  # Set the data frame with all the needed data.
  data = scores,
  
  # formula for residual covariance structure. This has to be set up so that each data point is associated with a unique residual. For example a multi-response model might have the R-structure defined by ~us(trait):units
  # covariances are useless here as variances are fixed to 1 (see prior).
  rcov =~ us(trait):units,
  
  # Set the (total) number of iterations, the burn-in, and the sampling frequency.
  nitt = my.nitt,
  burnin = my.burnin,
  thin = my.thin,
  verbose = TRUE,
  
  DIC = TRUE,
  
  # Store the posterior distributions of random effects and latent variables.
  pr = TRUE,
  pl = TRUE,
  
  # Force the threshold liabilities to range from -7 to 7 to facilitate
  # their estimation.
  trunc = TRUE
)

# Wait for 10 seconds after fitting has finished (just in case), create an output directory if missing, and save the resulting fit to an .Rda file.
Sys.sleep(10)

save(fit, file = paste(out.dir, chain_id, '.Rda', sep = ''))
print(paste0("output file saved: ", out.dir, chain_id, ".Rda"))

