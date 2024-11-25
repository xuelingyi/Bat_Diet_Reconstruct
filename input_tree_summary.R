## get libraries, functions
## load diets, diet.states, state.color
source("my.diet.functions.R")

## monophyletic nodes of interest (supported across all trees)
mono.nodes = c("Stenodermatinae", "Rhinophyllinae", "SR", 
               "Carolliinae", "Glyphonycterinae", "SRCG", 
               "Glossophaginae", "Phyllostominae", "Lonchophyllinae", "Lonchorhininae", "early_burst",
               "Micronycterinae", "EB_Mic", 
               "Desmodontinae", "Macrotinae",
               "Phyllostomidae", "Mormoops", "Pteronotus", "PM",
               "Noctilionidae", "Thyropteridae", "super_Noctilionoidea",
               "NM", 
               "Vespertilionidae", "Molossidae", "VM", 
               "Emballonuridae", "root")

## a subset of focus taxa
taxa = mono.nodes[!(mono.nodes %in% c("SR", "Mormoops", "Pteronotus", "super_Noctilionoidea", "VM"))]

## the input sif for each dataset is a subset from Table S1
sif = read.csv("vertlife_bat621.csv")

####################################################################################################################
############################### test monophyly and summarize topology #########################################
## download the full node-dated mammal tree from https://data.vertlife.org/
tree.files = read.table("mammal.trees.all")

## check monophyly of genus, subfamily, and family
mono_test = as.data.frame(matrix(nrow=length(c(unique(sif$Genus), unique(sif[sif$Family == "Phyllostomidae", "Subfamily"]), unique(sif$Family))), ncol=(10000 +2)))
names(mono_test) = c("level", "taxa", paste0("tree_", 1:10000))
mono_test$taxa = c(unique(sif$Genus), unique(sif[sif$Family == "Phyllostomidae", "Subfamily"]), unique(sif$Family))
mono_test[mono_test$taxa %in% sif$Genus, "level"] = "Genus"
mono_test[mono_test$taxa %in% sif$Subfamily, "level"] = "Subfamily"
mono_test[mono_test$taxa %in% sif$Family, "level"] = "Family"

## Prune trees using representatives of monophyletic subfamilies and families to summarize topologies
## Mormoopidae is not monophyletic and thus is represented by the two genera
topo_reps = c("Saccopteryx_leptura", "Eptesicus_fuscus", "Molossus_molossus",
               "Mystacina_tuberculata", 
               "Furipterus_horrens", "Mormoops_megalophylla", "Pteronotus_parnellii", "Pteronotus_davyi", "Noctilio_leporinus", "Thyroptera_tricolor",  
               "Carollia_perspicillata", "Desmodus_rotundus", "Glossophaga_soricina", "Glyphonycteris_daviesi", "Lionycteris_spurrelli", "Macrotus_waterhousii", "Micronycteris_megalotis", "Phyllostomus_discolor", "Rhinophylla_pumilio", "Sturnira_hondurensis", "Lonchorhina_aurita")

subf_trees = NULL
pg_trees = NULL
for(i in 1:10000){
  tree = read.tree(paste0("DNAonly_4098sp_topoFree_NDexp/", tree.files[i,]))
  tree = keep.tip(tree, sif$vertlife_name)
  
  # test monophyly
  for(t in mono_test$taxa) {
    my.tips = sif[sif$Family == t | sif$Subfamily == t | sif$Genus == t, "vertlife_name"]
    
    if(!(is.monophyletic(tree, tips=my.tips))){
      mono_test[mono_test$taxa == t, paste0("tree_",i)] = "not_monop"
   }
  }
  
  #summarize topology
  subf_trees = append(subf_trees, keep.tip(tree, topo_reps)) 
  pg_trees = append(pg_trees, keep.tip(tree, topo_reps[topo_reps %in% sif[sif$Family == "Phyllostomidae", "vertlife_name"]])) 
}

mono_test$count = apply(mono_test, 1, FUN=function(x){
  sum(!is.na(x[3:10002]))
})

### the subfamily-level topology (only considering divergence order of phyllostomid subfamilies)
pg.topology = unique.multiPhylo(pg_trees, use.edge.length = FALSE, use.tip.label = TRUE)

## count how many trees support each topology
pg.topology.summary = NULL
for(j in 1:length(pg.topology)){
  all.trees = NULL
  count.trees = 0
  
  for(i in 1:10000){
    if (all.equal.phylo(pg_trees[[i]], pg.topology[[j]], use.edge.length = FALSE)) {
      all.trees = paste0(all.trees, "; tree", i)
      count.trees = count.trees + 1
    }
  }
  
  pg.topology.summary = rbind(pg.topology.summary, c(all.trees, count.trees))
}
pg.topology.summary = as.data.frame(pg.topology.summary)
names(pg.topology.summary) = c("trees", "count")
pg.topology.summary$topology = row.names(pg.topology.summary)
pg.topology.summary$count = as.numeric(pg.topology.summary$count)
pg.topology.summary = pg.topology.summary[order(pg.topology.summary$count, decreasing = T),]

####################################################################################################################
############################### prep inputs for discrete trait reconstruction #########################################
tree.files = read.table("mammal.trees.all")

## RUN ONCE: a random sample of 100 trees
#trees = sample(1:10000, 100, replace = F)
#write.table(trees, "tree100_batch1.txt", row.names = F, quote = F, col.names = F)
trees = read.table("tree100_batch1.txt")

## sub-sample the 100 trees and label nodes of interest
taxa = c(unique(sif$Clade), unique(sif$Family), unique(sif$Subfamily))
taxa = taxa[taxa != "none"]
taxa = taxa[taxa != "Mormoopidae"]
taxa = c(taxa, "root")
taxa = c(taxa, "SRCG", "early_burst", "EB_Mic", "PM", "super_Noctilionoidea_M")
tree_100_sub = NULL
for(i in trees$V1){
  tree = read.tree(paste0("DNAonly_4098sp_topoFree_NDexp/", tree.files[i,]))
  tree = keep.tip(tree, sif$vertlife_name)
  tree$node.label <- (length(tree$tip.label) + 1):((length(tree$tip.label)) + tree$Nnode)
  
  # label the mrca of taxonomic groups  
  for(t in taxa) {
    if(t == "root"){
      my.tips = sif$vertlife_name
    } else {
      if(t == "SRCG"){ my.tips = sif[sif$Subfamily %in% c("Stenodermatinae", "Rhinophyllinae", "Carolliinae", "Glyphonycterinae"), "vertlife_name"] }
      if(t == "early_burst"){ my.tips = sif[sif$Family == "Phyllostomidae" & !(sif$Subfamily %in% c("Desmodontinae", "Micronycterinae", "Macrotinae")), "vertlife_name"] }
      if(t == "EB_Mic"){ my.tips = sif[sif$Family == "Phyllostomidae" & !(sif$Subfamily %in% c("Desmodontinae", "Macrotinae")), "vertlife_name"]  }
      if(t == "PM"){ my.tips = sif[sif$Family %in% c("Phyllostomidae", "Mormoopidae"), "vertlife_name"] }
      if(t == "super_Noctilionoidea_M"){ my.tips = sif[sif$Family %in% c("Phyllostomidae", "Mormoopidae", "Noctilionidae", "Furipteridae", "Thyropteridae", "Mystacinidae"), "vertlife_name"] }
      if(!(t %in% c("SRCG", "early_burst", "EB_Mic", "PM", "super_Noctilionoidea_M"))){
        my.tips = sif[sif$Clade == t | sif$Family == t | sif$Subfamily == t, "vertlife_name"] 
      }
    }

    # only label monophyletic nodes
    if(is.monophyletic(tree, my.tips)){
      tree$node.label[tree$node.label == getMRCA(tree, my.tips)] = t 
    }  else {
      print(t)
    }
  }
  
  ## save to the final tree set
  tree_100_sub = append(tree_100_sub, tree) 
}

diet.states = c("Absent" = "0", "Complementary" = "1", "Predominant" = "2", "Strict" = "3")
diet.in = diets[!(diets %in% c("Leaves.and.flower.pieces", "Seed"))]

## save as the inputs for MCMCglmm analyses
save(tree_100_sub, taxa, diet.in, diet.states, sif, file="bat621_100tree_batch1.RData")
## also for the other datasets sif176 (Phyllostomidae-focus) and sif179 (Phyllostomidae-focus including vampire bats)

####################################################################################################################
############################### prep inputs of the 23 trees having the genome subfamily topology #########################################
## get the 23 trees sharing the topology in the genome phylogeny (topology #11 in the summary)
tree_list = gsub("tree", "", unlist(strsplit(pg.topology.summary[pg.topology.summary$topology == 11, "trees"], split = "; ")))
tree_list = tree_list[tree_list != ""]
write.table(tree_list, "tree23_batch2.txt", row.names = F, quote = F, col.names = F)

tree.files = read.table("mammal.trees.all")
trees_topo11 = NULL
for(i in tree_list){
  tree = read.tree(paste0("DNAonly_4098sp_topoFree_NDexp/", tree.files[i,]))
  tree = keep.tip(tree, sif$vertlife_name)
  tree$node.label <- (length(tree$tip.label) + 1):((length(tree$tip.label)) + tree$Nnode)
  
  for(t in unique(sif$Family)){
    my.tips = sif[sif$Family == t, "vertlife_name"]
    
    if(is.monophyletic(tree, my.tips)) {
      node = getMRCA(tree, my.tips)
      tree$node.label[tree$node.label == node] = t
    }
  }
  
  ## label mrca nodes of subfamily
  for(t in unique(sif[sif$Family == "Phyllostomidae", "Subfamily"])){
    my.tips = sif[sif$Family == "Phyllostomidae" & sif$Subfamily == t, "vertlife_name"]
    
    if(is.monophyletic(tree, my.tips)) {
      node = getMRCA(tree, my.tips)
      tree$node.label[tree$node.label == node] = t
    }
  }
  
  ## label all monophyletic genus
  for(g in sif$Genus){
    if(is.monophyletic(tree, sif[sif$Genus == g, "vertlife_name"])){
      node = getMRCA(tree, sif[sif$Genus == g, "vertlife_name"])
      tree$node.label[tree$node.label == node] = g
    }
  }
  
  ## label the other internal nodes of interest
  for(t in mono.nodes[!(mono.nodes %in% c(sif$Family, sif$Subfamily, "Mormoops", "Pteronotus"))]){
    if(t == "SR"){ my.tips = sif[sif$Subfamily %in% c("Stenodermatinae", "Rhinophyllinae"), "vertlife_name"] }
    if(t == "SRCG"){ my.tips = sif[sif$Subfamily %in% c("Stenodermatinae", "Rhinophyllinae", "Carolliinae", "Glyphonycterinae"), "vertlife_name"] }
    if(t == "early_burst"){ my.tips = sif[sif$Family == "Phyllostomidae" & !(sif$Subfamily %in% c("Desmodontinae", "Micronycterinae", "Macrotinae")), "vertlife_name"] }
    if(t == "EB_Mic"){ my.tips = sif[sif$Family == "Phyllostomidae" & !(sif$Subfamily %in% c("Desmodontinae", "Macrotinae")), "vertlife_name"]  }
    if(t == "PM"){ my.tips = sif[sif$Family %in% c("Phyllostomidae", "Mormoopidae"), "vertlife_name"] }
    if(t == "super_Noctilionoidea"){ my.tips = sif[sif$Family %in% c("Phyllostomidae", "Mormoopidae", "Noctilionidae", "Furipteridae", "Thyropteridae"), "vertlife_name"] }
    if(t == "NM"){ my.tips = sif[sif$Family %in% c("Phyllostomidae", "Mormoopidae", "Noctilionidae", "Furipteridae", "Thyropteridae", "Mystacinidae"), "vertlife_name"] }
    if(t == "VM"){ my.tips = sif[sif$Family %in% c("Vespertilionidae", "Molossidae" ), "vertlife_name"] }
    if(t == "root"){ my.tips = sif$vertlife_name }
    
    my.tips = my.tips[!is.na(my.tips)]
    if(is.monophyletic(tree, my.tips)) {
      node = getMRCA(tree, my.tips)
      tree$node.label[tree$node.label == node] = t
    }
  }
  
  ## save to the final tree set
  trees_topo11 = append(trees_topo11, tree) 
}

### to make it easier with the downstream scripts, use the same file name although only 23 trees in this dataset
tree_100_sub = trees_topo11
save(tree_100_sub, sif, diet.in, diet.states, mono.nodes,
     file=paste0("sif176_tree23_topo11_diet", length(diet.in), ".RData"))

####################################################################################################################
############################### pPCA and the inputs for continuous trait reconstruction #########################################
# use the basal-inset coding
for(sp in sif[sif$conservative_coding != "NO", "vertlife_name"]){
  sif[sif$vertlife_name == sp, c("Arthropods", "Fruit")] = c(3, 0) 
  if(sp == "Micronycteris_microtis"){sif[sif$vertlife_name == sp, "Arthropods"] = 2 }
}

# use the three main diets
diet3 = c("Arthropods", "Fruit", "Pollen.and.nectar")

## ppca  
my.data = sif[, c("vertlife_name", diet3)]
row.names(my.data) = my.data$vertlife_name
my.data = my.data[diet3]
diet.in = c("PC1", "PC2")
for(i in 1:23){
  tree = tree_100_sub[[i]]
  ppca <- phyl.pca(tree, my.data, method = 'lambda', mode = 'cor', opt = 'REML')

  scores = as.data.frame(ppca$S)
  scores$vertlife_name = row.names(scores)
  
  loadings = as.data.frame(ppca$L)
  loadings$diet = row.names(loadings)

  ## save the inputs for mcmcglmm models
  save(tree, scores, ppca, loadings, diet.in, diet.states, mono.nodes,
       file=paste0("sif176_ppca_topo11_index_", i, ".RData"))
}

####################################################################################################################

