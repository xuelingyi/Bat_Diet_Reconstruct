## get libraries, functions
## load diets, diet.states, state.color
source("../../my.diet.functions.R")

load("sif179_100tree_batch1_diet6.RData")
# tree_100_sub, sif, diet.in, diet.states, mono.nodes

##########################################################
###################### more to load ###############################
### maybe load below if needed
load("sif179_100tree_batch1_topology.RData")
# topology100.summary, topology100

sif = read.csv("sif179.csv")
## only tips with the compiled diet data (original coding)
subfamilies = unique(sif[sif$Family == "Phyllostomidae", "Subfamily"])
## target taxa
ingroup.f = c("Phyllostomidae", "Mormoopidae", "Noctilionidae", "Furipteridae", "Thyropteridae", "Mystacinidae") 
outgroup.f = c("Vespertilionidae", "Molossidae", "Emballonuridae")

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

##########################################################
###### COPIED: subsample 100 trees - batch1 - topology - node.label ##################
tree_100_ids = sample(1:10000, 100, replace = F)
#8741:8743 %in% tree_100_ids
#[1] FALSE
## all subfamilies monophyletic
write.table(tree_100_ids, "tree_100_ids.txt", row.names = F, quote = F, col.names = F)
# renamed to tree100_batch1

## sub-sample 100 trees: label nodes, prep inputs for MCMC
taxa = c(subfamilies, ingroup.f[which(ingroup.f != "Mormoopidae")],
         "Mormoops", "Pteronotus",
         outgroup.f, "root")
## note: "Furipteridae" and "Mystacinidae" have only 1 sp, so no internal nodes labeled

tree_100_sub = NULL
for(i in tree_100_ids){
  tree = read.tree(paste0("~/Desktop/DNAonly_4098sp_topoFree_NDexp/", tree.files[i,]))
  tree = keep.tip(tree, sif$vertlife_name)
  
  tree$node.label <- (length(tree$tip.label) + 1):((length(tree$tip.label)) + tree$Nnode)
  
  # label the mrca of taxonomic groups  
  for(t in taxa) {
    if(t == "root"){
      my.tips = sif$vertlife_name
    } else {
      my.tips = sif[sif$Subfamily == t | sif$Family == t | sif$Genus == t, "vertlife_name"]
    }
    my.tips = my.tips[!is.na(my.tips)]
    
    if(!is.monophyletic(tree, my.tips)){print(paste0(t, "is not monophyly in tree ", i))} else {
      ## name the mrca node
      tree$node.label[tree$node.label == getMRCA(tree, my.tips)] = t 
    }
  }
  
  ## save to the final tree set
  tree_100_sub = append(tree_100_sub, tree) 
}
#ggtree(tree_100_sub[[1]]) + geom_tiplab(size=2) + geom_nodelab(size=2) + hexpand(0.3, 1)


## represents in the 100 trees
topo_reps = c("Carollia_perspicillata", "Desmodus_rotundus", "Glossophaga_soricina", "Glyphonycteris_daviesi", "Lionycteris_spurrelli", "Macrotus_waterhousii", "Micronycteris_megalotis", "Phyllostomus_discolor", "Rhinophylla_pumilio", "Sturnira_hondurensis", "Lonchorhina_aurita")
topo_rep_tree = NULL
for(i in 1:100){
  topo_rep_tree = append(topo_rep_tree, keep.tip(tree_100_sub[[i]], topo_reps))
}
topology100 = unique.multiPhylo(topo_rep_tree, use.edge.length = FALSE, use.tip.label = TRUE)
# 11 topologies
topology100.summary = NULL
for(j in 1:length(topology100)){
  all.trees = NULL
  count.trees = 0
  
  for(i in 1:100){
    if (all.equal.phylo(topo_rep_tree[[i]], topology100[[j]], use.edge.length = FALSE)) {
      all.trees = paste0(all.trees, "; tree", i)
      count.trees = count.trees + 1
    }
  }
  
  topology100.summary = rbind(topology100.summary, c(all.trees, count.trees))
}
topology100.summary = as.data.frame(topology100.summary)
names(topology100.summary) = c("trees", "count")
topology100.summary$topology = row.names(topology100.summary)
sum(as.numeric(topology100.summary$count))

pdf("sif179/topology100.pdf")
for(i in 1:length(topology100)){
  print(ggtree(groupOTU(topology100[[i]], 
                        list(Glossophaginae = "Glossophaga_soricina", 
                             Phyllostominae="Phyllostomus_discolor", 
                             Lonchorhininae="Lonchorhina_aurita")), 
               aes(color=group),
               branch.length="none") + 
          scale_color_manual(values=c("black", "orange", "royalblue3","purple4"),
                             name=paste0(as.numeric(topology100.summary[i, "count"]), "% trees")) +
          geom_tiplab(aes(label=gsub("_.*", "", label))) +
          hexpand(0.2,1))
}
dev.off()

sum(as.numeric(topology100.summary$count))

#load("sif179_100tree_batch1.RData")

diet.states = c("Absent" = "0", "Complementary" = "1", "Predominant" = "2", "Strict" = "3")
my.taxa = c("Stenodermatinae", "Glossophaginae",   "Carolliinae",      "Phyllostominae",   "Desmodontinae",   
            "Glyphonycterinae", "Lonchophyllinae",  "Micronycterinae",  "Lonchorhininae",   "Macrotinae" ,     
            "Rhinophyllinae" ,  "Phyllostomidae",   "Noctilionidae" ,   "Furipteridae" ,    "Thyropteridae" ,  
            "Mystacinidae",     "Mormoops",         "Pteronotus"  ,     "Vespertilionidae", "Molossidae" ,     
            "Emballonuridae" ,  "root")

save(tree_100_sub, topology100.summary, topology100,
     diets, diet.states, sif, file="sif179_100tree_batch1.RData")



### relabel the nodes
load("sif179_100tree_batch1.RData")

for(i in 1:100){
  
  ## label all monophyletic genus
  for(g in sif$Genus){
    if(is.monophyletic(tree_100_sub[[i]], sif[sif$Genus == g, "vertlife_name"])){
      node = getMRCA(tree_100_sub[[i]], sif[sif$Genus == g, "vertlife_name"])
      tree_100_sub[[i]]$node.label[tree_100_sub[[i]]$node.label == node] = g
    }
  }
  
  ## label taxa (same above)
  for(t in my.nodes[1:21]){
    my.tips = sif[sif$Subfamily == t | sif$Family == t | sif$Genus == t, "vertlife_name"]
    my.tips = my.tips[!is.na(my.tips)]
    
    if(!is.monophyletic(tree_100_sub[[i]], my.tips)) {
      print(paste0("not mono of ", t, " in tree ", i))
    } else {
      node = getMRCA(tree_100_sub[[i]], my.tips)
      tree_100_sub[[i]]$node.label[tree_100_sub[[i]]$node.label == node] = t
    }
  }
  
  ## label the other internal nodes 
  for(t in my.nodes[22:30]){
    if(t == "SR"){ my.tips = sif[sif$Subfamily %in% c("Stenodermatinae", "Rhinophyllinae"), "vertlife_name"] }
    if(t == "SRCG"){ my.tips = sif[sif$Subfamily %in% c("Stenodermatinae", "Rhinophyllinae", "Carolliinae", "Glyphonycterinae"), "vertlife_name"] }
    if(t == "early_burst"){ my.tips = sif[sif$Family == "Phyllostomidae" & !(sif$Subfamily %in% c("Desmodontinae", "Micronycterinae", "Macrotinae")), "vertlife_name"] }
    if(t == "EB_Mic"){ my.tips = sif[sif$Family == "Phyllostomidae" & !(sif$Subfamily %in% c("Desmodontinae", "Macrotinae")), "vertlife_name"]  }
    if(t == "PM"){ my.tips = sif[sif$Family %in% c("Phyllostomidae", "Mormoopidae"), "vertlife_name"] }
    if(t == "super_Noctilionoidea"){ my.tips = sif[sif$Family %in% c("Phyllostomidae", "Mormoopidae", "Noctilionidae", "Furipteridae", "Thyropteridae"), "vertlife_name"] }
    if(t == "NM"){ my.tips = sif[sif$Family %in% ingroup.f, "vertlife_name"] }
    if(t == "VM"){ my.tips = sif[sif$Family %in% c("Vespertilionidae", "Molossidae" ), "vertlife_name"] }
    if(t == "root"){ my.tips = sif$vertlife_name }
    
    my.tips = my.tips[!is.na(my.tips)]
    
    if(!is.monophyletic(tree_100_sub[[i]], my.tips)) {
      print(paste0("not mono of ", t, " in tree ", i))
    } else {
      node = getMRCA(tree_100_sub[[i]], my.tips)
      tree_100_sub[[i]]$node.label[tree_100_sub[[i]]$node.label == node] = t
    }
  }
  
}

save(tree_100_sub, topology100.summary, topology100,
     diets, diet.states, sif, file="sif179_100tree_batch1.RData")


## double-check that the labeled nodes are the same nodes across trees
# "Furipteridae" "Mystacinidae" are single-taxa --> no internal node labels
my.nodes.label = my.nodes[!(my.nodes %in% c("Furipteridae", "Mystacinidae"))]
for(i in 2:100){
  cmp = comparePhylo(tree_100_sub[[1]], tree_100_sub[[i]])
  n = cmp$NODES
  n[,1] = gsub(" \\(.*", "", n[,1])
  n[,2] = gsub(" \\(.*", "", n[,2])
  if(!all(my.nodes.label %in% n[n[,1] == n[,2], 1])){
    print(i)
  }
}
# checked!

##################################################################
####################### RUNONCE: prep inputs ################################
diet.data = read.csv("../../diet_data.csv")
## updated 20241004

## ONCE
#load("sif179_100tree_batch1.RData")
## update sif and diets in the load data; this was the saved RData used in versions/diet_run1 
#save(topology100.summary, topology100, file="sif179_100tree_batch1_topology.RData")
#load("versions/diet7/sif179_100tree_batch1_diet7.RData")

length(unique(diet.data$Taxon))
#276
diet.data$Taxon = gsub("_", " ", diet.data$Taxon)
all(sif$ScientificName %in% diet.data$Taxon)
# TRUE

for(i in 1:nrow(sif)){
  if(any(sif[i, diets] != diet.data[diet.data$Taxon == sif[i, "ScientificName"], diets])){
    ## these species have different diets now
    print("changed diets")
    print(sif[i, "ScientificName"])
  }
  if(sif[i, "Family"] != diet.data[diet.data$Taxon == sif[i, "ScientificName"], "Family"]){print(sp)}
  if(sif[i, "Genus"] != diet.data[diet.data$Taxon == sif[i, "ScientificName"], "Genus"]){print(sp)}
  if(sif[i, "Family"] == "Phyllostomidae"){
    if(sif[i, "Subfamily"] != diet.data[diet.data$Taxon == sif[i, "ScientificName"], "Subfamily"]){print(sp)}
  }
}
#[1] "changed diets" compared to versions diet_run1
#[1] "Chiroderma villosum"
#[1] "Macrotus waterhousii"
#[1] "Micronycteris megalotis"
#[1] "Micronycteris microtis"

## taxonomy is identical
sif = sif[, c("ScientificName", "vertlife_name","genes", "gene_list")]
sif = merge(sif, diet.data[, c("Family", "Subfamily", "Genus", "Taxon", diets, "conservative_coding", "source")], by.x="ScientificName", by.y="Taxon", all.x=T)

for(d in diets){
  print(paste0(d, ": ", paste0(unique(sif[,d]), collapse = " ")))
}
#[1] "Arthropods: 0 2 3 1"
#[1] "Blood: 0 3"
#[1] "Terrestrial.vertebrates: 0 2 1"
#[1] "Fish: 0 1 2"
#[1] "Leaves.and.flower.pieces: 0 1"
#[1] "Pollen.and.nectar: 0 2 1"
#[1] "Fruit: 3 1 2 0"
#[1] "Seed: 0 2"

unique(sif[sif$Leaves.and.flower.pieces !=0, c("Pollen.and.nectar", "Fruit")])
#Pollen.and.nectar Fruit
#7                  1     2
#54                 2     2
#57                 2     1
#74                 1     1
#80                 0     1
diet.in = diets[diets != "Leaves.and.flower.pieces"]

unique(sif[sif$Seed !=0, c("vertlife_name", "Pollen.and.nectar", "Fruit")])
#         vertlife_name Pollen.and.nectar Fruit
#25   Chiroderma_doriae                 1     1
#28 Chiroderma_villosum                 1     2
sif[sif$vertlife_name == "Chiroderma_doriae", "source"] = paste0(sif[sif$vertlife_name == "Chiroderma_doriae", "source"], "_fruit2_in_mcmc")
sif[sif$Seed != 0, "Fruit"] = 2
unique(sif[sif$Seed !=0, c("vertlife_name", "Pollen.and.nectar", "Fruit")])
#vertlife_name Pollen.and.nectar Fruit
#25   Chiroderma_doriae                 1     2
#28 Chiroderma_villosum                 1     2
write.csv(sif, "sif179.csv", row.names = F, quote = F)
diet.in = diet.in[diet.in != "Seed"]

save(tree_100_sub, sif, diet.in, diet.states, mono.nodes,
     file=paste0("sif179_100tree_batch1_diet", length(diet.in), ".RData"))


## test monophyly
for(i in 1:100){
  tree = tree_100_sub[[i]]
  if(!is.monophyletic(tree, c("Desmodus_rotundus", "Diaemus_youngi"))){
    print(i)
  }
  if(!is.monophyletic(tree, c("Desmodus_rotundus", "Diaemus_youngi"))){
    print(i)
  }
}
## all show "Desmodus_rotundus" and "Diaemus_youngi" as sister
##########################################################
########################### densi tree conserve ####################################
## only keep the mono node lab
for(i in 1:100){
  ## relabel the species with changed diets
  for(j in 1:179){
    if(tree_100_sub[[i]]$tip.label[j] %in% sif[sif$conservative_coding != "NO", "vertlife_name"]) {
      tree_100_sub[[i]]$tip.label[j] = paste0(tree_100_sub[[i]]$tip.label[j], " *")
    }
  }
  
  for(j in 1:178){
    if(!(tree_100_sub[[i]]$node.label[j] %in% mono.nodes)){
      tree_100_sub[[i]]$node.label[j] = j+179
    }
  }
}

#conservative coding on these species 
sif[sif$conservative_coding != "NO", "vertlife_name"]
#[1] "Lonchorhina_aurita"        "Macrotus_californicus"     "Macrotus_waterhousii"     
#[4] "Micronycteris_brosseti"    "Micronycteris_matses"      "Micronycteris_megalotis"  
#[7] "Micronycteris_microtis"    "Micronycteris_schmidtorum"
sif[sif$conservative_coding != "NO", "vertlife_name"] = paste0(sif[sif$conservative_coding != "NO", "vertlife_name"], " *")

##if(diet.coding == "conserve") {
for(sp in sif[sif$conservative_coding != "NO", "vertlife_name"]){
  sif[sif$vertlife_name == sp, c("Arthropods", "Fruit")] = c(3, 0) 
  if(sp == "Micronycteris_microtis"){sif[sif$vertlife_name == sp, "Arthropods"] = 2 }
}
## 

diet.heatmap = sif[, diet.in]
rownames(diet.heatmap) = sif$vertlife_name
for(i in 1:ncol(diet.heatmap)){
  diet.heatmap[,i] = as.character(diet.heatmap[,i])
}

p0 = ggdensitree(tree_100_sub, alpha=.3, colour='steelblue') + 
  geom_tiplab(size=1) + geom_nodelab(aes(label=gsub('([0-9])', '', label)), size=3) +
  hexpand(.3) + vexpand(.3)

pdf("sif179_100tree_densi_diet6_conserve.pdf", width = 7, height = 10)
gheatmap(p0, diet.heatmap, offset=7, width=0.2, font.size=2.5, color="black",
         colnames_angle=90, colnames_position = "top", hjust = 0) + 
  scale_fill_manual(values=setNames(state.color, diet.states), name="Diet code") +
  theme(legend.key.size = unit(0.8, "line"),
        legend.text = element_text(size=6), legend.title = element_text(size=8),
        legend.position = "inside", legend.position.inside = c(0.1, 0.6),
        legend.background = element_blank())
dev.off()


png("sif179_100tree_densi.png", width = 7, height = 10, res = 600, units = "in")
ggdensitree(tree_100_sub, alpha=.3, colour='steelblue') + 
  hexpand(.3) + vexpand(.3)
dev.off()

pdf("tree1_diet6_conserve.pdf", width = 7, height = 10)
gheatmap(ggdensitree(tree_100_sub[1]) + geom_tiplab(size=1) + geom_nodelab(aes(label=gsub('([0-9])', '', label)), size=3) + hexpand(.3) + vexpand(.3), 
         diet.heatmap, offset=7, width=0.2, font.size=2.5, color="black",
         colnames_angle=90, colnames_position = "top", hjust = 0) + 
  scale_fill_manual(values=setNames(state.color, diet.states), name="Diet code") +
  theme(legend.key.size = unit(0.8, "line"),
        legend.text = element_text(size=6), legend.title = element_text(size=8),
        legend.position = "inside", legend.position.inside = c(0.1, 0.6),
        legend.background = element_blank())
dev.off()

##################################################################
####################### tree27 (Rojas2018 topology) #################################
## pie plots only using trees of certain topology
#load("sif179_100tree_batch1_topology.RData")
#rojas.topo.trees = unlist(strsplit(topology100.summary[topology100.summary$topology==8, "trees"], split="; "))[-1]
## [1] "tree27" only one tree in the batch 1 set shows this topology 
trees = read.table("tree100_batch1.txt")
id = trees[27, "V1"]
rm(trees)
tree = tree_100_sub[[27]]
my.node.label = tree$node.label
## redo node labels
tree$node.label <- (length(tree$tip.label) + 1):((length(tree$tip.label)) + tree$Nnode)

## get results
for(diet.code in c("raw", "conserve")){
  ## RUNONCE
  #system(paste0("scp $delta:/projects/project-xyi/diet_reconstruct/vertlife/sif179/trees_", diet.code,
  #              "/tree_", id, "/node_prob.csv tree_", id, "_", diet.code, ".csv"))
  
  node_prob = read.csv(paste0("tree_", id, "_", diet.code, ".csv"))
  tips = sif
  tips$label = tips$vertlife_name
  
  pdf(paste0("tree_", id, "_", diet.code, "_MCMC_results.pdf"), width = 12, height = 25)
  for(d in diet.in){
    prob = node_prob[, c("Taxon", grep(d, names(node_prob), value=T))]
    names(prob) = gsub(paste0(d, "_Prob_"), "", names(prob))
    prob$node = prob$Taxon
    
    ## change label back to number 
    prob$label = sapply(prob$node, FUN = function(x){
      if(x %in% tree$tip.label){
        return(which(tree$tip.label == x))
      } else{
        return(which(my.node.label == x) + 179)
      }
    })
    prob$node = prob$label
    
    pies <- nodepie(prob, cols = 2:(ncol(prob)-2))
    pies <- lapply(pies, function(g) g+scale_fill_manual(values = setNames(state.color, names(diet.states))))
    
    tips$diet = as.vector(sapply(tips[,d], FUN=function(x){
      names(diet.states)[diet.states == x]
    }))
    
    tree2 = left_join(tree, tips[,c("label", "diet")], by="label")
    p0 = ggtree(tree2) + 
      geom_tiplab(size=1.5, offset=1) + 
      ## tips will be covered by the pie charts; this is only to get the legends
      geom_tippoint(aes(color=diet), size=0.0001) +
      scale_color_manual(values = setNames(state.color, names(diet.states)), name = d, breaks=names(diet.states)) +
      guides(color = guide_legend(override.aes = list(size = 1.5))) +
      hexpand(0.2, 1)
    
    p1 = p0 + geom_inset(pies, width = 0.03, height = 0.03) +
      theme(legend.key.size = unit(0.8, "line"),
            legend.text = element_text(size=6), legend.title = element_text(size=8),
            legend.position = "inside", legend.position.inside = c(0.2, 0.8),
            legend.background = element_blank())
    
    print(p1)
  }
  
  dev.off()
}

## moved
setwd("tree_459")
for(diet.code in c("raw", "conserve")){
  node_prob = read.csv(paste0("tree_459_", diet.code, ".csv"))
  for(col in 2:20){
    node_prob[,col] = round(node_prob[,col], 3)
  }
  names(node_prob) = gsub("_Prob", "", names(node_prob))
  write.csv(node_prob[node_prob$Taxon %in% c("early_burst", "Phyllostomidae", "PM"),],paste0(diet.code, ".csv"), row.names = F, quote = F)
}
setwd("../")
##########################################################
####################### check tip species ################################

for(diet.code in c("raw", "conserve")){
  load(paste0("mcmc_aggregate_diet6_", diet.code, "_1_100.RData"))
  
  ## save the estimated tip states
  sif_estimate = as.data.frame(matrix(nrow=179, ncol=(length(diet.in) + 1)))
  sif_estimate[,1] = sif$vertlife_name
  names(sif_estimate) = c("vertlife_name", paste0(diet.in, "_", diet.code))

  for(d in diet.in){
    for(i in 1:179){
      sp = sif_estimate[i, "vertlife_name"]
      sif_estimate[i, paste0(d, "_", diet.code)] = get_tip_posterior(sp, d)
    }
  }
  
  if(!all(sif$vertlife_name == sif_estimate$vertlife_name)){
    print("species not in the same order!!")
  } else {
    row_check = NULL
    diet.mismatch = NULL
    
    for(i in 1:179){
      
      # strict should not have other diets
      if(any(sif_estimate[i, ] == 3)){
        ### this species has a strict diet
        
        if(sum(sif_estimate[i, ] == 3) > 1){
          print("more than one strict diet! check!")
          print(sif_estimate[i, ])
        } else {
          ## only one strict diet; all the other diets should be absent
          if(!all(sif_estimate[i, -1] %in% c(0, 3))){
            print(sif_estimate[i, ])
          }
        }
        
      }
      
      ## compare posterior with the input
      if(!all(sif[i, diet.in] == sif_estimate[i, paste0(diet.in, "_", diet.code)])){
        row_check = c(row_check, i)
        
        mismatch = which(sif[i, diet.in] != sif_estimate[i, paste0(diet.in, "_", diet.code)])
        diet.mismatch = c(diet.mismatch, diet.in[mismatch])
        for(d in diet.in[mismatch]){
          if(!(sif_estimate[i, paste0(d, "_", diet.code)] %in% (sif[i, d]-1):(sif[i, d]+1))){
            print("the esimate and input states differ more than 1")
            print(sif$vertlife_name[i])
            print(sif_estimate[i, paste0(diet.in, "_", diet.code)])
            print(sif[i, diet.in])
          }
        }
      }
    }
    
    ### DONOTUSE this plot is too large and having lots of uninformative values (ie identical)
    ## plot the sp that have different input and posterior
    #posterior = sif_estimate[row_check, c("vertlife_name", paste0(diet.in, "_", diet.code))]
    #names(posterior) = c("vertlife_name", diet.in)
    #posterior$vertlife_name = paste0(posterior$vertlife_name, "_post")
    #all = rbind(sif[row_check, c("vertlife_name", diet.in)], posterior)
    #all = melt(all, id="vertlife_name")
    #all$value = as.character(all$value)
    #all = all[order(all$vertlife_name),]
    #pdf(paste0("tip_post_mismatch_", diet.code, ".pdf"), height = 12, width = 6)
    #ggplot(data = all, aes(x = variable, y = vertlife_name)) +
    #  geom_tile(aes(fill = value), color="black") +
    #  scale_fill_manual(values = setNames(state.color, diet.states), name="Diet code") +
    #  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle=90))
    #dev.off()
  }
  
  write.table(row_check, paste0("sif_mismatch_rows_", diet.code, ".txt"), row.names = F, quote = F)
  write.csv(sif_estimate, paste0("sif_estimate_", diet.code, ".csv"), row.names = F, quote = F)
  assign(paste0("sif_estimate_", diet.code), sif_estimate)
}

# 1. species with most likely states having prob < 0.5
## raw
#[1] "Chrotopterus_auritus, Terrestrial.vertebrates: Complementary 0.463; Predominant 0.283"
#[1] "Trachops_cirrhosus, Terrestrial.vertebrates: Complementary 0.496; Absent 0.311"
#[1] "Vampyrum_spectrum, Terrestrial.vertebrates: Complementary 0.465; Predominant 0.273"
# conserve
#[1] "Macrophyllum_macrophyllum, Arthropods: Strict 0.492; Predominant 0.488"
#[1] "Trachops_cirrhosus, Arthropods: Predominant 0.496; Complementary 0.484"
#[1] "Chrotopterus_auritus, Terrestrial.vertebrates: Complementary 0.462; Predominant 0.283"
#[1] "Mimon_cozumelae, Terrestrial.vertebrates: Absent 0.499; Complementary 0.426"
#[1] "Trachops_cirrhosus, Terrestrial.vertebrates: Complementary 0.496; Absent 0.308"
#[1] "Vampyrum_spectrum, Terrestrial.vertebrates: Complementary 0.464; Predominant 0.275"
#[1] "Artibeus_fraterculus, Pollen.and.nectar: Complementary 0.5; Absent 0.499"
#[1] "Phyllostomus_discolor, Fruit: Complementary 0.498; Predominant 0.491"


# 2. species having strict diet and another diet (same in raw and conserve)
#vertlife_name      Arthropods_raw Blood_raw Terrestrial.vertebrates_raw Fish_raw Pollen.and.nectar_raw Fruit_raw
#Diphylla_ecaudata         1        3 (>0.5 as this is not reported above)   0       0            0         0


## 3. the esimate and input states differ more than 1" (same in raw and conserve)
#[1] "Glyphonycteris_daviesi"
#Arthropods_raw Blood_raw Terrestrial.vertebrates_raw Fish_raw Pollen.and.nectar_raw Fruit_raw
# 2             0                  0         0                     0         1
#Arthropods Blood Terrestrial.vertebrates Fish Pollen.and.nectar Fruit
# 2             0                  2         0                 0     1


## 4. table(diet.mismatch)
#diet.mismatch (raw/conserve)
#Arthropods                   Fish                   Fruit        Pollen.and.nectar 
#22 / 27                      5 / 5                  16 / 26              12 / 12
#Terrestrial.vertebrates 
#14 / 14


## focus on Arthropods, Fruit, and nectar
## plot the tips that differ
sif = read.csv("sif179.csv")
for(diet.code in c("raw", "conserve")){
  load(paste0("mcmc_aggregate_diet6_", diet.code, "_1_100.RData"))
  estimates = read.csv(paste0("sif_estimate_", diet.code, ".csv"))
  #all(estimates$vertlife_name == sif$vertlife_name) 
  #row = read.table(paste0("sif_mismatch_rows_", diet.code, ".txt"), header=T)

  if(diet.code == "conserve") {
    for(sp in sif[sif$conservative_coding != "NO", "vertlife_name"]){
      sif[sif$vertlife_name == sp, c("Arthropods", "Fruit")] = c(3, 0) 
      if(sp == "Micronycteris_microtis"){sif[sif$vertlife_name == sp, "Arthropods"] = 2 }
    }
  }
  
  pdf(paste0("tip_mismatch_", diet.code, ".pdf"))
  sp_25 = NULL
  for(d in c("Arthropods", "Pollen.and.nectar", "Fruit", "Terrestrial.vertebrates", "Fish", "Blood")){
    species = sif$vertlife_name[sif[,d] != estimates[,paste0(d, "_", diet.code)]]
    if(length(species) >0){
      data.all = NULL
      for(sp in species){
        data = get_tip_posterior(sp, d, prob = T)
        data$species =sp
        if(data[data$state == names(diet.states)[diet.states == sif[sif$vertlife_name ==sp, d]], "prob"] < 0.25){
          sp_25 = c(sp_25, sp)
        }
        data.all = rbind(data.all, data)
      }
      print(plot_tip_bar(d, data.all))
    } else {
      #print(paste0("all species have the same input and esitmate states in ", d))
      print(ggplot()+theme_void())
    }
  }
  dev.off()
  assign(paste0("sp_25_", diet.code), sp_25)
}
# species with posterior probability of the input state < 0.25 
## raw has 23 species
## conserve has 22 species
#sp_25_raw[!(sp_25_raw %in% sp_25_conserve)]
#"Vampyressa_nymphaea"

## total species with different input and posterior states (determined by the state having the highest prob)
## raw has 57 species
## conserve has 50 species


load("mcmc_aggregate_diet6_raw_1_100.RData")
load("mcmc_aggregate_diet6_conserve_1_100.RData")
for(sp in sif[sif$Blood == 3, "vertlife_name"]){
  print(sp)
  print(get_tip_posterior(sp, "Blood", prob = T))
}


## only plot the tips that posterior of the input state is <0.25
sif = read.csv("sif179.csv")
for(diet.code in c("raw", "conserve")){
  load(paste0("mcmc_aggregate_diet6_", diet.code, "_1_100.RData"))
  estimates = read.csv(paste0("sif_estimate_", diet.code, ".csv"))
  
  if(diet.code == "conserve") {
    for(sp in sif[sif$conservative_coding != "NO", "vertlife_name"]){
      sif[sif$vertlife_name == sp, c("Arthropods", "Fruit")] = c(3, 0) 
      if(sp == "Micronycteris_microtis"){sif[sif$vertlife_name == sp, "Arthropods"] = 2 }
    }
  }

  for(d in c("Arthropods", "Pollen.and.nectar", "Fruit", "Terrestrial.vertebrates", "Fish")){
    species = sif$vertlife_name[sif[,d] != estimates[,paste0(d, "_", diet.code)]]
    if(length(species) >0){
      data.all = NULL
      for(sp in species){
        data = get_tip_posterior(sp, d, prob = T)
        data$species =sp
        if(data[data$state == names(diet.states)[diet.states == sif[sif$vertlife_name ==sp, d]], "prob"] < 0.25){
          data.all = rbind(data.all, data)
        }
      }
      assign(paste0(d, "_plot"), plot_tip_bar(d, data.all) + 
              scale_y_discrete(expand = c(0, 0)) + theme_bw() +
              theme(panel.grid = element_line(color="white")))
    } 
  }
  
  pdf(paste0("tip_mismatch_25_", diet.code, ".pdf"), width = 12, height = 6)
  print(ggarrange(nrow=2, ncol=2,
            Arthropods_plot, Fruit_plot, 
            Terrestrial.vertebrates_plot, 
            ggarrange(nrow=3, ncol=1, heights = c(1,1,0.5),
                      Pollen.and.nectar_plot, Fish_plot, 
                      ggplot()+geom_text(aes(x=1, y=1, label=diet.code), size=10) + theme_void())))
  dev.off()
}


##########################################################
####################### pie plot of nodes (including all subf) ################################

## pie plots across all trees & node states
for(diet.code in c("raw", "conserve")){
  node_estimate = as.data.frame(matrix(nrow=length(mono.nodes), ncol=(length(diet.in) + 1)))
  node_estimate[,1] = mono.nodes
  names(node_estimate) = c("mono.nodes", paste0(diet.in, "_", diet.code))
  
  load(paste0("mcmc_aggregate_diet6_", diet.code, "_1_100.RData"))
  
  pdf(paste0("mono.node_pie_diet6_", diet.code, ".pdf"), height = 3, width=5)
  for(n in mono.nodes){
    
    for(d in diet.in){
      assign(paste0(d, ".plot"), get_node_pie_plot(n, d))
      node_estimate[node_estimate$mono.nodes == n, paste0(d, "_", diet.code)] = get_tip_posterior(n, d)
    }

    p=ggarrange(nrow=2, ncol=3,
                Arthropods.plot, Pollen.and.nectar.plot, Fruit.plot, 
                Terrestrial.vertebrates.plot, Fish.plot, Blood.plot)
    
    print(annotate_figure(p, top=paste0(n, "_", diet.code)))
  }
  dev.off()
  
  write.csv(node_estimate, paste0("mono_node_estimate_", diet.code, ".csv"), row.names = F, quote = F)
  assign(paste0("node_estimate_", diet.code), node_estimate)
  
  for(n in c("root", "PM", "Phyllostomidae", "early_burst")){
    for(d in c("Arthropods", "Pollen.and.nectar", "Fruit")){
      print(paste0(n, "_", d))
      data = get_tip_posterior(n, d, prob = T)
      data[,3] = round(data[,3], 3)
      print(data[order(data[,1]),c(1,3)])
    }
  } 

}
#conserve
#"Noctilionidae, Arthropods: Predominant 0.492; Complementary 0.488"

## see which mono.nodes are different in raw and conserve
raw = read.csv("mono_node_estimate_raw.csv")
conserve = read.csv("mono_node_estimate_conserve.csv")
all(raw$mono.nodes == conserve$mono.nodes)
for(i in 1:length(mono.nodes)){
  if(any(raw[i, ] != conserve[i,])){
    print(raw[i, c(1, which(raw[i, ] != conserve[i,]))])
    print(conserve[i, c(1, which(raw[i, ] != conserve[i,]))])
  }
}
#mono.nodes Arthropods_raw
#13     EB_Mic              1
#mono.nodes Arthropods_conserve
#13     EB_Mic                   2
#mono.nodes Arthropods_raw Fruit_raw
#15 Macrotinae              2         1
#mono.nodes Arthropods_conserve Fruit_conserve
#15 Macrotinae                   3              0
#mono.nodes Fruit_raw
#16 Phyllostomidae         1
#mono.nodes Fruit_conserve
#16 Phyllostomidae              0
#mono.nodes Arthropods_raw
#22 super_Noctilionoidea              2
#mono.nodes Arthropods_conserve
#22 super_Noctilionoidea                   3
#mono.nodes Arthropods_raw
#23         NM              2
#mono.nodes Arthropods_conserve
#23         NM                   3



load("mcmc_aggregate_diet6_raw_1_100.RData")
load("mcmc_aggregate_diet6_conserve_1_100.RData")
get_node_pie_plot("Lonchorhininae", "Arthropods", prob = T)

load("mcmc_aggregate_diet6_raw_GENUS_1_100.RData")
load("mcmc_aggregate_diet6_conserve_GENUS_1_100.RData")
get_node_pie_plot("Micronycteris", "Arthropods", prob = T)

##########################################################
######################## genus mrca ###########################
## pie plots across all trees & node states

for(diet.code in c("raw", "conserve")){
  
  load(paste0("mcmc_aggregate_diet6_", diet.code, "_GENUS_1_100.RData"))
  # summary, not_mono
  
  genus_estimate = as.data.frame(matrix(nrow=length(unique(Arthropods_summary_aggregate$Taxon)), ncol=(length(diet.in) + 1)))
  genus_estimate[,1] = unique(Arthropods_summary_aggregate$Taxon)
  names(genus_estimate) = c("Genus", paste0(diet.in, "_", diet.code))
  
  pdf(paste0("Genus_pie_diet6_", diet.code, ".pdf"), height = 3, width=5)
  for(n in genus_estimate$Genus){
    
    for(d in diet.in){
      assign(paste0(d, ".plot"), get_node_pie_plot(n, d))
      genus_estimate[genus_estimate$Genus == n, paste0(d, "_", diet.code)] = get_tip_posterior(n, d)
    }
    if(n %in% not_mono$Genus){
      if(nrow(not_mono[not_mono$Genus == n,]) + nrow(Arthropods_summary_aggregate[Arthropods_summary_aggregate$Taxon == n, ]) != 100){
        print(paste0(n, " rows do not sum up!!"))
      } else {
        print(paste0(n, "is esitmated using ", nrow(Arthropods_summary_aggregate[Arthropods_summary_aggregate$Taxon ==n, ]), " trees"))
      }
    }
    
    p=ggarrange(nrow=2, ncol=3,
                Arthropods.plot, Pollen.and.nectar.plot, Fruit.plot, 
                Terrestrial.vertebrates.plot, Fish.plot, Blood.plot)
    
    print(annotate_figure(p, top=paste0(n, "_", diet.code)))
  }
  dev.off()
  
  write.csv(genus_estimate, paste0("Genus_estimate_", diet.code, ".csv"), row.names = F, quote = F)
  assign(paste0("Genus_estimate_", diet.code), genus_estimate)
}
#[1] "Choeroniscusis esitmated using 99 trees"
#[1] "Pipistrellusis esitmated using 7 trees"
#conserve
#"Noctilionidae, Arthropods: Predominant 0.492; Complementary 0.488"

## these genera are not monophyletic in any trees, so they were not reconstructed
unique(not_mono$Genus)[!(unique(not_mono$Genus) %in% unique(Arthropods_summary_aggregate$Taxon))]
#[1] "Lonchophylla" "Mimon"        "Phyllostomus"



## see which mono.nodes are different in raw and conserve
all(Genus_estimate_conserve$Genus == Genus_estimate_raw$Genus)
for(i in 1:nrow(Genus_estimate_raw)){
  if(any(Genus_estimate_raw[i, ] != Genus_estimate_conserve[i,])){
    print(Genus_estimate_raw[i, ])
    print(Genus_estimate_conserve[i,])
  }
}
#       Genus Arthropods_raw Blood_raw Terrestrial.vertebrates_raw Fish_raw Pollen.and.nectar_raw Fruit_raw
# Macrotinae              2           0                           0             0                     0         1
#Genus Arthropods_conserve Blood_conserve Terrestrial.vertebrates_conserve Fish_conserve Pollen.and.nectar_conserve Fruit_conserve
# Macrotinae              3           0                           0             0                     0              0
#Genus Arthropods_raw Blood_raw Terrestrial.vertebrates_raw Fish_raw Pollen.and.nectar_raw Fruit_raw
# Micronycteris              2             0                           0        0                     0         1
#Genus Arthropods_conserve Blood_conserve Terrestrial.vertebrates_conserve Fish_conserve Pollen.and.nectar_conserve Fruit_conserve
# Micronycteris              2             0                           0        0                      0        0

##########################################################


########################### major topology plots #################################
load("sif179_100tree_batch1_topology.RData")
load("sif179_all_10k/sif179_10k_trees_monophyly_internal.RData")
load("sif179_all_10k/sif179_10k_trees_monophyly.RData")

trees = read.table("tree100_batch1.txt")
### these are the 41 topologies of the 10k trees
pg.topology.summary$trees100 = sapply(pg.topology.summary$trees, FUN=function(x){
   data = unlist(strsplit(x, split = "; "))
   return(paste0(data[data %in% paste0("tree", trees$V1)], collapse = "; "))
})

nrow(pg.topology.summary[pg.topology.summary$trees100 != "", ])
## 11
pg.topology.summary = pg.topology.summary[order(pg.topology.summary$count, decreasing = T),]
pg.topology.summary = pg.topology.summary[1:6,]
sum(pg.topology.summary$count)/10000
#  0.9471


sif.t = sif[sif$vertlife_name %in% pg_trees[[1]]$tip.label,]
sif.t$label = sif.t$vertlife_name
sif.t$color_group = sapply(sif.t$Subfamily, FUN=function(x){
  if(x %in% c("Carolliinae","Rhinophyllinae", "Stenodermatinae", "Glyphonycterinae")){
    return("fruit")
  }
  if(x %in% c("Glossophaginae", "Lonchophyllinae")){return("nectar")}
  if(x == "Desmodontinae"){return("blood")}
  if(x == "Phyllostominae"){return("omnivore")}
  if(x %in% c( "Lonchorhininae",   "Macrotinae", "Micronycterinae" )){return("insect")}
})

for(i in 1:6){
  topo = pg.topology[[as.numeric(pg.topology.summary[i, "topology"])]]
  topo = left_join(topo, sif.t, by="label")
  
  p10k = paste0(round((pg.topology.summary[i, "count"] / 10000)*100, 2), "%")
  p100 = paste0(length(unlist(strsplit(pg.topology.summary[i, "trees100"], split = "; "))), "%")
  
  assign(paste0("topo", i), 
         ggtree(topo, aes(color=color_group)) + geom_tiplab(aes(label=Subfamily)) + hexpand(0.2,1) + xlim(0,50)+
    scale_color_manual(values=c("fruit"="#028484", "nectar"="orange", "blood"="red3", "insect"="black", "omnivore"="purple4"),guide="none") + 
      labs(title=paste0("topology in ", p10k, " & ", p100, " trees"))
    )
}


pdf("topology_final.pdf", width = 10, height = 5)
ggarrange(nrow=2, ncol=3,
          topo1, topo2, topo3, 
          topo4, topo5, topo6)
dev.off()

##################################################################
########################## some more test ####################################
load("mcmc_aggregate_diet6_conserve_1_100.RData")

for (i in 1:100){
  tree = tree_100_sub[[i]]
  if(!is.monophyletic(tree, sif[sif$Genus == "Mormoops" | sif$Family == "Phyllostomidae", "vertlife_name"])){
    print(i)
  }
}
# 51

#tree$node.label getMRCA(tree, sif[sif$Genus == "Mormoops" | sif$Family == #"Phyllostomidae", "vertlife_name"])
get_node_pie_plot("PM", "Arthropods", prob = T)

##################################################################
##################################################################

So for example why we have Desmodus and Diaemus (from Desmodontinae) in the middle of (Phyllostominae), i.e. Macrotus and Macrophyllum in Xueling ordering.

sif[sif$Genus %in% c("Desmodus", "Diaemus", "Macrotus", "Macrophyllum"), c("Family", "Subfamily", "Genus")]


##################################################################
##################################################################