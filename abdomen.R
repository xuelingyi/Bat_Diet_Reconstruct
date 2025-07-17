# modified from https://github.com/BPerezLamarque/ABDOMEN

#data=elton143_diet5_tree23
#tree_id=$i
#iter=20000
#warmup=10000
#chains=4
#prior_Z0=empirical
#seed=1
#name=tree.${tree_id}.${prior_Z0}
#Rscript abdomen_sif_delta.R $data $tree_id $iter $warmup $chains $prior_Z0 $seed $name > log.${name}

##################################
library(ggplot2)
library(mvMORPH)
library(RPANDA)
library(rstan, lib="~/R_4.4")
rstan_options(auto_write = TRUE)
library(RColorBrewer)

source("ABDOMEN.R")
## available from https://github.com/BPerezLamarque/ABDOMEN
  
args <- commandArgs(TRUE)
data = args[1]
load(paste0(data, ".RData"))

tree_id = as.numeric(args[2])
tree = force.ultrametric(tree_100_sub[[tree_id]])

iter <-  as.numeric(args[3]) # total number of iterations in STAN
warmup <- as.numeric(args[4]) # number of warmup iterations in STAN
chains <-  as.numeric(args[5])  # number of chains for the inference
nb_cores = chains  # number of cores to run the analyses
prior_Z0 = args[6]
seed <- as.numeric(args[7])
name <- args[8]


code_path <- getwd() # where the ABDOMEN plots will be generated.
detection_threshold <- 0.01 # i.e., less than 1%. the detection threshold: below this threshold, we assume that we cannot detect the diet (very very rare or practially absent). Then, all relative abundances below this threshold are set to this threshold. 

sd_prior_logY <- 2
mean_prior_logY <- 0 

# hard coded thin=1
fit_summary <- ABDOMEN(tree, table, name, 
                       code_path = code_path, prior_Z0 = prior_Z0,
                       detection_threshold = detection_threshold, seed = seed, 
                       mean_prior_logY = mean_prior_logY, sd_prior_logY = sd_prior_logY,
                       nb_cores = nb_cores, chains = chains, warmup = warmup, iter = iter)


# The ancestral diet composition estimated for each node of the phylogenetic tree  (Z0_nodes):
Z0_nodes <- ABDOMEN_extract_Z0_nodes(tree, table, fit_summary, detection_threshold = detection_threshold) 

## change to mono.node names
Z0_nodes$tree_nodes = sapply(Z0_nodes$node, FUN=function(x){
  return(tree$node.label[as.numeric(x) - length(tree$tip.label)])
})
Z0_nodes$tree_id = tree_id
write.table(Z0_nodes[, c(grep("elton", names(Z0_nodes), value=T), "tree_nodes", "tree_id")], paste0(name, "_Z0_nodes.tsv"), row.names=F, quote=F, sep="\t")

save(list=ls(), file=paste0(name, ".RData"))


