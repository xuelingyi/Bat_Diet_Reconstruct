## get libraries, functions
## load diets, diet.states, state.color
source("my.diet.functions.R")

load("sif176_100tree_batch1_diet5.RData")
tree_list = read.table("tree100_batch1.txt")

## if using the basal-insect coding
diet.code = "conserve"
for(sp in sif[sif$conservative_coding != "NO", "vertlife_name"]){
  sif[sif$vertlife_name == sp, c("Arthropods", "Fruit")] = c(3, 0) 
  if(sp == "Micronycteris_microtis"){sif[sif$vertlife_name == sp, "Arthropods"] = 2 }
}

##################################################################
########################### plot trees and diets ####################################
## only keep the mono node lab & label the species with changed diets
for(i in 1:23){
  for(j in 1:(nrow(sif)-1)){
    if(!(tree_100_sub[[i]]$node.label[j] %in% mono.nodes)){
      tree_100_sub[[i]]$node.label[j] = j+nrow(sif)
    }
  }
}
sif[sif$conservative_coding != "NO", "ScientificName"] = paste0(sif[sif$conservative_coding != "NO", "ScientificName"], " *")
sif$label = sif$vertlife_name
tree = left_join(tree_100_sub[[1]], sif[, c("label", "ScientificName")], by="label")

diet.heatmap = sif[, diet.in]
rownames(diet.heatmap) = sif$vertlife_name
for(i in 1:ncol(diet.heatmap)){
  diet.heatmap[,i] = as.character(diet.heatmap[,i])
}

## diets with tree1
## duplicate tree1 to use the densitree function 
tree = append(tree, tree)
gheatmap(ggdensitree(tree) + geom_tiplab(aes(label=ScientificName), size=1.8) + 
           geom_nodelab(aes(label=gsub('([0-9])', '', label)), size=3) + 
           hexpand(.3) + vexpand(0.05), 
         diet.heatmap, offset=12, width=0.15, font.size=2.5, color="black",
         colnames_angle=90, colnames_position = "top", hjust = 0) + 
  scale_fill_manual(values=setNames(state.color, diet.states), name="Diet code", 
                    labels = names(diet.states)) +
  theme(legend.key.size = unit(0.8, "line"),
        legend.text = element_text(size=6), legend.title = element_text(size=8),
        legend.position = "inside", legend.position.inside = c(0.1, 0.6),
        legend.background = element_blank())

## densitree of all trees
ggdensitree(tree_100_sub, alpha=.3, colour='steelblue') + 
  hexpand(.3) + vexpand(0.05)


##################################################################
####################### node pie charts ################################
## only the nodes of interest
taxa = mono.nodes[!(mono.nodes %in% c("SR", "Mormoops", "Pteronotus", "super_Noctilionoidea", "VM"))]

for(diet.code in c("raw", "conserve")){
  node_estimate = as.data.frame(matrix(nrow=length(mono.nodes), ncol=(length(diet.in) + 1)))
  node_estimate[,1] = mono.nodes
  names(node_estimate) = c("mono.nodes", paste0(diet.in, "_", diet.code))
  
  load(paste0("mcmc_aggregate_diet", length(diet.in), "_", diet.code, "_1_100.RData"))
  for(n in mono.nodes){
    
    for(d in diet.in){
      assign(paste0(d, ".plot"), get_node_pie_plot(n, d, title="no_food", lab.size=5))
      node_estimate[node_estimate$mono.nodes == n, paste0(d, "_", diet.code)] = get_tip_posterior(n, d)
    }
    
    assign(paste0(n, ".plot"), 
           ggarrange(nrow=1, ncol=6,
                     ggplot()+ geom_text(aes(x=1, y=1, label=n), size=2.5) + theme_void(),
                     Arthropods.plot, Terrestrial.vertebrates.plot, Fish.plot, Pollen.and.nectar.plot, Fruit.plot)
           )
  }
  assign(paste0("node_estimate_", diet.code), node_estimate)
  names.plot = ggarrange(nrow=1, ncol=6,
                         ggplot()+ geom_text(aes(x=1, y=1, label=paste0("Node", "\n", diet.code)), size=2.5) + theme_void(), 
                         ggplot()+ geom_text(aes(x=1, y=1, label="Arthropods"), size=2.5) + theme_void(),
                         ggplot()+ geom_text(aes(x=1, y=1, label="Terrestrial\nvertebrates"), size=2.5) + theme_void(),
                         ggplot()+ geom_text(aes(x=1, y=1, label="Fish"), size=2.5) + theme_void(),
                         ggplot()+ geom_text(aes(x=1, y=1, label="Pollen.nectar"), size=2.5) + theme_void(),
                         ggplot()+ geom_text(aes(x=1, y=1, label="Fruit"), size=2.5) + theme_void())
  
  print(ggarrange(ncol = 1, nrow = (length(taxa)+1),
            names.plot,
            Stenodermatinae.plot, Rhinophyllinae.plot, Carolliinae.plot, Glyphonycterinae.plot,
            SRCG.plot,
            Glossophaginae.plot, Phyllostominae.plot, Lonchophyllinae.plot, Lonchorhininae.plot,
            early_burst.plot, Micronycterinae.plot, EB_Mic.plot, Macrotinae.plot,
            Phyllostomidae.plot, PM.plot, 
            Noctilionidae.plot, Thyropteridae.plot, NM.plot, 
            Vespertilionidae.plot, Molossidae.plot, Emballonuridae.plot, root.plot))
}
dev.off()

### tip of families Mystacinidae and Furipteridae
for(f in c("Mystacinidae", "Furipteridae")){
  sp = sif[sif$Family == f, "vertlife_name"]
  for(d in diet.in[c(1,4,5)]){
    assign(paste0(d, "_plot"), get_node_pie_plot(sp, d))
  }
  assign(paste0(sp, "_plot"), ggarrange(nrow=1, ncol = 4, 
                                        ggplot()+geom_text(aes(x=1, y=1), label=sp, size=3) + theme_void(),
                                        Arthropods_plot, Pollen.and.nectar_plot, Fruit_plot
  ))
}

##################################################################
####################### compare input and output for tips ################################
diet.code = "conserve"
for(sp in sif[sif$conservative_coding != "NO", "vertlife_name"]){
    sif[sif$vertlife_name == sp, c("Arthropods", "Fruit")] = c(3, 0) 
    if(sp == "Micronycteris_microtis"){sif[sif$vertlife_name == sp, "Arthropods"] = 2 }
  }

load(paste0("mcmc_aggregate_diet5_", diet.code, "_1_100.RData"))
  
sif_estimate = as.data.frame(matrix(nrow=176, ncol=(length(diet.in) + 1)))
sif_estimate[,1] = sif$vertlife_name
names(sif_estimate) = c("vertlife_name", paste0(diet.in, "_", diet.code))
for(d in diet.in){
    for(i in 1:176){
      sp = sif_estimate[i, "vertlife_name"]
      sif_estimate[i, paste0(d, "_", diet.code)] = get_tip_posterior(sp, d)
    }
}

## focus on the main diets
for(d in c("Arthropods", "Pollen.and.nectar", "Fruit")){
  ## species with different input and output
  species = sif$vertlife_name[sif[,d] != sif_estimate[,paste0(d, "_", diet.code)]]
  
    data.all = NULL
    for(sp in species){
      data = get_tip_posterior(sp, d, prob = T)
      data$species = sif[sif$vertlife_name ==sp, "ScientificName"]
      data.all = rbind(data.all, data)
    }
    assign(paste0(d, "_data"), data.all)
    assign(paste0(d, "_plot"), plot_tip_bar_HPD(d, data.all, tip.lab = "ScientificName") + labs(title=d))
}

ggarrange(nrow=1, ncol=2, widths = c(1.2, 1),
          Arthropods_plot + theme(legend.position = "bottom", legend.title = element_blank()), 
          ggarrange(nrow=2, ncol=1, heights = c(1.5, 1),
                    Fruit_plot + theme(legend.position = "none"), 
                    Pollen.and.nectar_plot + theme(legend.position = "none"), 
                    labels = c("B", "C")),
          labels = c("A", ""))


##########################################################
############################# extract node probabilities per topology ##################################
### topology in the 100 trees with 179 bats (including vampires)
## all subfamilies are monophyletic in all trees; prune to representative species
load("sif179_100tree_batch1_diet6.RData")
topo_reps = c("Carollia_perspicillata", "Desmodus_rotundus", "Glossophaga_soricina", "Glyphonycteris_daviesi", "Lionycteris_spurrelli", "Macrotus_waterhousii", "Micronycteris_megalotis", "Phyllostomus_discolor", "Rhinophylla_pumilio", "Sturnira_hondurensis", "Lonchorhina_aurita")
topo_rep_tree = NULL
for(i in 1:100){
  topo_rep_tree = append(topo_rep_tree, keep.tip(tree_100_sub[[i]], topo_reps))
}
topology100 = unique.multiPhylo(topo_rep_tree, use.edge.length = FALSE, use.tip.label = TRUE)
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
topology100.summary = topology100.summary[order(topology100.summary$count, decreasing = T),]


diet.code = "conserve"
load(paste0("mcmc_aggregate_diet5_", diet.code, "_1_100.RData"))

topo_diet_summary = NULL
for(i in 1:11){
  trees = gsub("tree", "", unlist(strsplit(topology100.summary[i, "trees"], split="; ")))
  trees = trees[trees != ""]
  ## these are index, so correspond to row numbers
  trees = tree_list[trees,]
  
  for (n in mono.nodes){
    for(d in diet.in){
      data = get_node_pie_plot(n, d, prob = T, trees=trees)
      topo_diet_summary = rbind(topo_diet_summary, 
                                c(n, d, as.character(data[1, c(1,3)]), paste0("topo100_", topology100.summary[i, "topology"])))
      
      data = get_node_pie_plot(n, d, prob = T)
      topo_diet_summary = rbind(topo_diet_summary, 
                                c(n, d, as.character(data[1, c(1,3)]), "topo100_all"))
    }
  }
}
topo_diet_summary = as.data.frame(topo_diet_summary)
names(topo_diet_summary) = c("node", "diet", "state", "prob", "topology")

state.color2[4] = "grey20"
taxa = mono.nodes[!(mono.nodes %in% c("SR", "Mormoops", "Pteronotus", "super_Noctilionoidea", "VM"))]
data = topo_diet_summary[topo_diet_summary$node %in% taxa,]
data$node = factor(data$node, levels=taxa)
data = data[order(data$node),]
for(d in diet.in){
  assign(paste0(d, "_plot"), ggplot() + 
           geom_point(data=data[data$diet == d & data$topology != "topo100_all",],
                      aes(x=node, y=as.numeric(prob), color=state), alpha=0.8, shape=4) + 
           scale_color_manual(values=setNames(state.color2, names(diet.states)), name=d) +
           geom_point(data=data[data$diet == d & data$topology == "topo100_all",],
                      aes(x=node, y=as.numeric(prob), color=state), alpha=0.8, shape=21, size=3) + 
           labs(y="Prob") +
           theme_bw() + 
           theme(axis.text.x = element_text(angle=90)) 
  )
}

print(ggarrange(nrow = 3, ncol=1, 
                labels = c("A", "B", "C"), heights = c(1,1,1.5),
                Arthropods_plot + labs(x="") + theme(axis.text.x = element_blank(), 
                                                     axis.ticks.x = element_blank()), 
                Pollen.and.nectar_plot + labs(x="") + theme(axis.text.x = element_blank(), 
                                                           axis.ticks.x = element_blank()), 
                Fruit_plot))

