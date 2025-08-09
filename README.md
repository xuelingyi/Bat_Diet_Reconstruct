# Reconstruct ancestral diets in bats
These scripts are used for bat ancestral diet reconstructions described in the following paper:

Yi, X., Kontopoulos, D. G., & Hiller, M. (2025). **Comprehensive phylogenetic trait estimations support ancestral omnivory in the ecologically diverse bat family Phyllostomidae.** Evolution, qpaf154. https://doi.org/10.1093/evolut/qpaf154

# input data
The input diets and taxonomy updates are given in the supplementary Table S1. The input phylogenetic trees are the node-dated “DNA-only” mammalian trees downloaded from VertLife http://vertlife.org/phylosubsets (Upham et al., 2019). Indexes of the trees used in our analyses are provided in **data.zip** which includes: 
* **mammal.trees.all**: Names of all 10000 downloaded trees from Upham et al. 2019. Note that IDs in tree names range between 0-9999.
* **tree100_batch1.txt**: The index of the 100 randomly sampled trees from the 10000 total trees in Upham et al. 2019. Note that these indexes range between 1-10000.
* **tree23_batch2.txt**: The index of the 23 trees that have the same subfamily topology as the genome phylogeny. Indexes range between 1-10000.

The script **input_tree_summary.R** is used to test clade monophyly, summarize tree topology, and prepare input files for MCMCglmm and ABDOMEN models. 
The summarized data are loaded for the following analyses.
* sif179_100tree_batch1_diet6.RData: the Phyllostomidae-focused analysis including vampire bats across 100 trees
* sif176_100tree_batch1_diet5.RData: the Phyllostomidae-focused analysis without vampire bats across 100 trees
* sif176_tree23_topo11_diet5.RData: the Phyllostomidae-focused analysis without vampire bats, using only the 23 trees having the optimal subfamily topology
* bat621_100tree_batch1.RData: the all-bat analysis without vampire bats across 100 trees
* elton143_diet5_tree23.RData: the Phyllostomidae-focused analysis using compositional diets for ABDOMEN       

# ancestral reconstructions 
The scripts for ancestral reconstructions using MCMCglmm are adapted from https://github.com/dgkontopoulos/Kontopoulos_et_al_torpor_evolution_2025 (Kontopoulos et al. 2025) using the above input data and the following scripts (1-4). Ancestral reconstructions using continuous compositional data follow the ABDOMEN (Perez-Lamarque et al. 2023) tutorial https://github.com/BPerezLamarque/ABDOMEN.
1. **vertlife_MCMC_diet5.R**: fit the MCMCglmm model on five diets, used for the datasets sif176 and bat621. **vertlife_MCMC_diet6.R** includes blood feeding and vampire bats (sif179). **vertlife_MCMC_pca.R** fits MCMCglmm models on the first two PCs summarized using pPCA of the three main diets (sif176 with 23 trees). 
2. **check_ESS_PSRF_modified.R** and **check_ESS_PSRF_modified_ppca.R**: check chain convergence for each tree based on the effective sample size (ESS) and potential scale reduction factor (PSRF) of each model parameter.
3. **vertlife_extract_threshold_probability.R** and **vertlife_extract_threshold_probability_ppca**: for each tree, summarize results across chains.
4. **aggregate.R**: aggregate results of the labeled node of interest across trees. For a node not labeled in the input, **aggregate_one_node.R** can be used to get its estimated diet states across chains and trees. pPCA results are aggregated and plotted using **aggregate_ppca.R**
5. **phyloEM.R**: automatically detect adaptive shifts using the ppca scores.
6. **abdomen.R**: run continuous compositional models using ABDOMEN.

The script **output_summary.R** was used to visualize and compare results (e.g., tip input versus output states, results across topologies). 

# references
* Kontopoulos, D. G., Levesque, D. L., & Hiller, M. (2025). Numerous independent gains of daily torpor and hibernation across endotherms, linked with adaptation to diverse environments. Functional Ecology. https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2435.14739
* Upham, N. S., Esselstyn, J. A., & Jetz, W. (2019). Inferring the mammal tree: Species-level sets of phylogenies for questions in ecology, evolution, and conservation. PLOS Biology, 17(12), e3000494. https://doi.org/10.1371/journal.pbio.3000494
* Perez-Lamarque, B., Sommeria-Klein, G., Duret, L., & Morlon, H. (2023). Phylogenetic comparative approach reveals evolutionary conservatism, ancestral composition, and integration of vertebrate gut microbiota. Molecular Biology and Evolution, 40(7), msad144. https://doi.org/10.1093/molbev/msad144

