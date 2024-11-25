# Reconstruct ancestral diets in bats
These scripts are used for bat ancestral diet reconstructions described in the manuscript [to add bioRxiv link]. 

# input data
The input diets and taxonomy updates are given in the supplementary Table S1. The input phylogenetic trees are the node-dated “DNA-only” mammalian trees downloaded from VertLife http://vertlife.org/phylosubsets (Upham et al., 2019). Indexes of the trees used in our analyses are provided in **data.zip** which includes: 
* **mammal.trees.all**: Names of all 10000 downloaded trees from Upham et al. 2019. Note that IDs in tree names range between 0-9999.
* **tree100_batch1.txt**: The index of the 100 randomly sampled trees from the 10000 total trees in Upham et al. 2019. Note that these indexes range between 1-10000.
* **tree23_batch2.txt**: The index of the 23 trees that have the same subfamily topology as the genome phylogeny. Indexes range between 1-10000.

The script **tree_summary_prep_input.R** is used to test clade monophyly, summarize tree topology, and prepare input files for MCMCglmm models. 

# scripts for ancestral reconstructions 
The scripts for ancestral reconstructions using MCMCglmm are adapted from https://github.com/dgkontopoulos/Kontopoulos_et_al_torpor_evolution_2024 (Kontopoulos et al. 2023). The ancestral diet reconstructions have three major steps, using the above input data and the following scripts:
1. **vertlife_MCMC_diet5.R**: fit the MCMCglmm model on five diets, used for the datasets sif176 and bat621. **vertlife_MCMC_diet6.R** includes blood feeding and vampire bats (sif179). **vertlife_MCMC_pca.R** fits MCMCglmm models on the first two PCs summarized using pPCA of the three main diets (sif176 with 23 trees). 
2. **check_ESS_PSRF_modified.R** and **check_ESS_PSRF_modified_ppca.R**: check chain convergence for each tree based on the effective sample size (ESS) and potential scale reduction factor (PSRF) of each model parameter.
3. **vertlife_extract_threshold_probability.R** and **vertlife_extract_threshold_probability_ppca**: for each tree, summarize results across chains.
4. **aggregate.R**: aggregate results of the labeled node of interest across trees. For a node not labeled in the input, **aggregate_one_node.R** can be used to get its estimated diet states across chains and trees. pPCA results are aggregated and plotted using **aggregate_ppca.R**

The script **output_summary.R** was used to visualize and compare results (e.g., tip input versus output states, results across topologies). 

# references
* Kontopoulos, D.-G., Levesque, D. L., & Hiller, M. (2023). Numerous independent gains of torpor and hibernation across endotherms, linked with adaptation to diverse environments (p. 2023.12.12.571278). bioRxiv. https://doi.org/10.1101/2023.12.12.571278
* Upham, N. S., Esselstyn, J. A., & Jetz, W. (2019). Inferring the mammal tree: Species-level sets of phylogenies for questions in ecology, evolution, and conservation. PLOS Biology, 17(12), e3000494. https://doi.org/10.1371/journal.pbio.3000494 

