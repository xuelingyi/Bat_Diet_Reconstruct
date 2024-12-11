library(PhylogeneticEM)
library(phytools)
library(ggtree)
library(ggplot2)
library(tidyr)
library(ggpubr)

##### basal-diet coding 
# ppca data have been generated in the script input_tree_summary.R

shifts=NA
for(i in 1:23) {
  load(paste0("inputs/sif176_ppca_topo11_index_", i, ".RData"))

  ## use PC1 and PC2
  data = as.data.frame(t(scores[, 1:2]))
  tree.nnls<-force.ultrametric(tree)
  
  res = PhyloEM(Y_data = data, ## data
                phylo = tree.nnls, ## phylogeny
                process = "scOU", ## scalar OU
                random.root = TRUE, ## root is stationary
                stationary.root = TRUE,
                K_max = 20, ## maximal number of shifts
                parallel_alpha = TRUE, ## parallelize on alpha values
                Ncores = 6)
  
  ## Plot k
  print(plot_criterion(res, main=paste0("tree_", i)))
  ## Plot selected solution (LINselect)
  print(plot(res, show.tip.label=T, label_cex = 0.5, margin_plot = c(0,0,0.5,0),
       main=paste0("tree_", i)))

  tips.nodes = c(tree.nnls$tip.label, tree.nnls$node.label)
  ## shift detected on the branches leading to these nodes/tips
  shifts = c(shifts, list(tips.nodes[tree.nnls$edge[params_process(res)$shifts$edges,2]]))
  assign(paste0("res_", i), res)
}
shifts = shifts[2:24]

data = as.data.frame(table(unlist(shifts)))
names(data) = c("node", "count")
data$freq = data$count / 23
data = data[order(data$freq, decreasing = T),]
data$node = factor(data$node, levels = data$node)
## only extract labeled nodes (non-numeric) that are shared across trees
data = data[is.na(extract_numeric(data$node)),]

p1 = ggplot(data) + geom_bar(aes(x=node, y=freq), stat="identity") + 
  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle=45)) 

######## raw coding
### run ppca on raw coding of the three main diets
load("sif176_tree23_topo11_diet5.RData")
tree_list = read.table("tree23_batch2.txt")
diet3 = c("Arthropods", "Fruit", "Pollen.and.nectar")

my.data = sif[, c("vertlife_name", diet3)]
row.names(my.data) = my.data$vertlife_name
my.data = my.data[diet3]
diet.in = c("PC1", "PC2")

shifts=NA
for(i in 1:23){
tree = tree_100_sub[[i]]
ppca <- phyl.pca(
  tree, 
  my.data, 
  method = 'lambda', 
  mode = 'cor', opt = 'REML')

scores = as.data.frame(ppca$S)
data = as.data.frame(t(scores[, 1:2]))
tree.nnls<-force.ultrametric(tree)

res = PhyloEM(Y_data = data, ## data
              phylo = tree.nnls, ## phylogeny
              process = "scOU", ## scalar OU
              random.root = TRUE, ## root is stationary
              stationary.root = TRUE,
              K_max = 20, ## maximal number of shifts
              parallel_alpha = TRUE, ## parallelize on alpha values
              Ncores = 6)

print(plot_criterion(res, main=paste0("tree_", i)))
print(plot(res, show.tip.label=T, label_cex = 0.5, margin_plot = c(0,0,0.5,0),
           main=paste0("tree_", i)))

tips.nodes = c(tree.nnls$tip.label, tree.nnls$node.label)
shifts = c(shifts, list(tips.nodes[tree.nnls$edge[params_process(res)$shifts$edges,2]]))
assign(paste0("res_", i), res)
}

shifts = shifts[2:24]

data = as.data.frame(table(unlist(shifts)))
names(data) = c("node", "count")
data$freq = data$count / 23
data = data[order(data$freq, decreasing = T),]
data$node = factor(data$node, levels = data$node)
data = data[is.na(extract_numeric(data$node)),]

p2 = ggplot(data) + geom_bar(aes(x=node, y=freq), stat="identity") + 
  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle=45)) 


ggarrange(ncol=2, nrow=1, 
          p1 + theme(axis.text.x = element_text(angle=90)), 
          p2 + theme(axis.text.x = element_text(angle=90)), 
          labels = c("A)", "B)"))


