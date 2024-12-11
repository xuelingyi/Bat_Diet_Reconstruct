#https://cran.r-project.org/web/packages/PhylogeneticEM/vignettes/tutorial.html
#https://pbastide.github.io/PhylogeneticEM/articles/monkeys.html
library(PhylogeneticEM)
library(phytools)
library(ggtree)
library(ggplot2)
library(tidyr)
library(ggpubr)

#################################################################
########################## phyloEM ###################################
## PC1,2,3 results are not good
shifts=NA
pdf("phyloEM_conserve_PC12_tree11.pdf", width = 8, height = 8)

for(i in 1:23) {
  load(paste0("inputs/sif176_ppca_topo11_index_", i, ".RData"))
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
  
  ## Plot selected solution (LINselect)
  print(plot_criterion(res, main=paste0("tree_", i)))
  #plot(res, method.selection = "DDSE")
  print(plot(res, show.tip.label=T, label_cex = 0.5, margin_plot = c(0,0,0.5,0),
       main=paste0("tree_", i)))
  #plot(equivalent_shifts(tree.nnls, params_process(res)))
  #plot(res, show.tip.label=T, label_cex = 0.5, margin_plot = c(0,0,0,0),params = params_process(res, K = 5))
  
  tips.nodes = c(tree.nnls$tip.label, tree.nnls$node.label)
  ## shift detected on the branches leading to these nodes/tips
  shifts = c(shifts, list(tips.nodes[tree.nnls$edge[params_process(res)$shifts$edges,2]]))
  assign(paste0("res_", i), res)
}
dev.off()
shifts = shifts[2:24]

data = as.data.frame(table(unlist(shifts)))
names(data) = c("node", "count")
data$freq = data$count / 23
data = data[order(data$freq, decreasing = T),]
data$node = factor(data$node, levels = data$node)

data = data[is.na(extract_numeric(data$node)),]

p1 = ggplot(data) + geom_bar(aes(x=node, y=freq), stat="identity") + 
  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle=45)) 

save(list = c(ls(pattern = "res_"), "shifts"), file="phyloEM_conserve_PC12_tree11.RData")

#################################################################
########################### phyloEM - raw #################################
### run ppca on raw diets
load("../topo11_mcmc/sif176_tree23_topo11_diet5.RData")
tree_list = read.table("../topo11_mcmc/tree100_batch2.txt")
diet3 = c("Arthropods", "Fruit", "Pollen.and.nectar")

my.data = sif[, c("vertlife_name", diet3)]
row.names(my.data) = my.data$vertlife_name
my.data = my.data[diet3]
diet.in = c("PC1", "PC2")


shifts=NA
pdf("phyloEM_raw_PC12_tree11.pdf", width = 8, height = 8)

for(i in 1:23){
tree = tree_100_sub[[i]]
ppca <- phyl.pca(
  tree, 
  my.data, 
  method = 'lambda', 
  mode = 'cor', opt = 'REML')

scores = as.data.frame(ppca$S)
#scores$vertlife_name = row.names(scores)

#loadings = as.data.frame(ppca$L)
#loadings$diet = row.names(loadings)
#save(tree, scores, ppca, loadings, diet.in, diet.states, mono.nodes,
#     file=paste0("topo11_ppca/inputs/sif176_ppca_topo11_index_", i, ".RData"))
#scores = merge(scores, sif, sort=F)
#ggplot(scores) + 
#  geom_point(aes(x=PC1-1, y=PC2, color=as.character(Arthropods))) +
#  geom_point(aes(x=PC1, y=PC2, color=as.character(Fruit))) +
#  geom_point(aes(x=PC1+1, y=PC2, color=as.character(Pollen.and.nectar))) +
  #scale_color_manual(values=setNames(state.color, diet.states), name="Arthropods_Fruit_Nectar") +
#  labs(x=paste0("PC1: ", round(ppca$Eval[1,1]/sum(ppca$Eval)*100, 2), "%"),
#       y=paste0("PC2: ", round(ppca$Eval[2,2]/sum(ppca$Eval)*100, 2), "%")) +
#  theme_bw()
#}

## use pca for shift modeling
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

## Plot selected solution (LINselect)
print(plot_criterion(res, main=paste0("tree_", i)))
#plot(res, method.selection = "DDSE")
print(plot(res, show.tip.label=T, label_cex = 0.5, margin_plot = c(0,0,0.5,0),
           main=paste0("tree_", i)))
#plot(equivalent_shifts(tree.nnls, params_process(res)))
#plot(res, show.tip.label=T, label_cex = 0.5, margin_plot = c(0,0,0,0),params = params_process(res, K = 5))

tips.nodes = c(tree.nnls$tip.label, tree.nnls$node.label)
## shift detected on the branches leading to these nodes/tips
shifts = list(shifts, tips.nodes[tree.nnls$edge[params_process(res)$shifts$edges,2]])
assign(paste0("res_", i), res)
}
dev.off()

## shifts were formatted wrongly
shifts=NA
for(i in 1:23){
  tree = tree_100_sub[[i]]
  tree.nnls<-force.ultrametric(tree)
  res = get(paste0("res_", i))
  tips.nodes = c(tree.nnls$tip.label, tree.nnls$node.label)
  
  shifts = c(shifts, list(tips.nodes[tree.nnls$edge[params_process(res)$shifts$edges,2]]))
}
shifts = shifts[2:24]

data = as.data.frame(table(unlist(shifts)))
names(data) = c("node", "count")
data$freq = data$count / 23
data = data[order(data$freq, decreasing = T),]
#data[data$freq > 0.5, ]
data$node = factor(data$node, levels = data$node)

data = data[is.na(extract_numeric(data$node)),]

p2 = ggplot(data) + geom_bar(aes(x=node, y=freq), stat="identity") + 
  ylim(0,1) + theme_bw() + theme(axis.text.x = element_text(angle=45)) 

save(list = c(ls(pattern = "res_"), "shifts"), file="phyloEM_raw_PC12_tree11.RData")


pdf("phyloEM_plot.pdf", width = 8, height = 5)
png("phyloEM_plot.png", width = 8, height = 5, res=600, units="in")
ggarrange(ncol=2, nrow=1, 
          p1 + theme(axis.text.x = element_text(angle=90)), 
          p2 + theme(axis.text.x = element_text(angle=90)), 
          labels = c("A)", "B)"))
dev.off()
#################################################################
#################################################################
########################## phyloEM diets ###################################
load("../topo11_mcmc/sif176_tree23_topo11_diet5.RData")
## conservative coding
for(sp in sif[sif$conservative_coding != "NO", "vertlife_name"]){
  sif[sif$vertlife_name == sp, c("Arthropods", "Fruit")] = c(3, 0) 
  if(sp == "Micronycteris_microtis"){sif[sif$vertlife_name == sp, "Arthropods"] = 2 }
}
rownames(sif) = sif$vertlife_name
data = as.data.frame(t(sif[, diet.in]))
res2 <- PhyloEM(Y_data = data, ## data
                phylo = tree.nnls, ## phylogeny
                process = "scOU", ## scalar OU
                random.root = TRUE, ## root is stationary
                stationary.root = TRUE,
                K_max = 10, ## maximal number of shifts
                parallel_alpha = TRUE, ## parallelize on alpha values
                Ncores = 4)
# More than half of the observations lie on a hyperplane

data = as.data.frame(t(sif[, c("Fruit", "Pollen.and.nectar")]))
res3 <- PhyloEM(Y_data = data, ## data
                phylo = tree.nnls, ## phylogeny
                process = "scOU", ## scalar OU
                random.root = TRUE, ## root is stationary
                stationary.root = TRUE,
                K_max = 10, ## maximal number of shifts
                parallel_alpha = TRUE, ## parallelize on alpha values
                Ncores = 4)

pdf("phyloEM_test3.pdf")
## Plot selected solution (LINselect)
plot(res2, show.tip.label=T, label_cex = 0.2) 
plot(res3, show.tip.label=T, label_cex = 0.2) 
dev.off()


#################################################################
########################## mvMORPH ###################################
#https://cran.r-project.org/web/packages/mvMORPH/vignettes/How_to_use_mvMORPH.pdf
library(mvMORPH)

i=1
set.seed(i)
load(paste0("inputs/sif176_ppca_topo11_index_", i, ".RData"))

## require continuous traits; use PC1 and 2
data = scores[, c("PC1", "PC2")]
#fitting of multivariate multiple rates of evolution under a Brownian Motion model.
## "BM1" for a unique rate of evolution per trait.
fit_BM <- mvBM(tree, data, model="BM1")
## "BMM" for multi-rate and multi-selective regimes
## No selective regimes mapped on the tree, only a BM1 model could be estimated
fit_BMM <- mvBM(tree, data, model="BMM")
## fits to a multivariate dataset of continuous traits a multivariate Early Burst (EB) 
fit_EB <- mvEB(tree, data)
## fitting of a multivariate Ornstein-Uhlenbeck (OU1) model with possibly multiple optima (OUM) for different "selective regimes". 
fit_OU1 <- mvOU(tree, data, model="OU1")


### compare likelihoods 
LRT(fit_BM, fit_BMM)
#p=1
LRT(fit_EB, fit_OU1)

## compare models using Akaike weights:
results <- list(fit_BM, fit_BMM, fit_EB, fit_OU1)
weights <- aicw(results)
weights
#-- Akaike weights -- 
#.            Rank  AIC   diff     wi     AICw
#OU1 4            1 2066  0.0 1.00e+00 1.00e+00
#BM1 default 1    2 2096 29.9 3.26e-07 3.26e-07
#BM1 default 2    3 2096 29.9 3.26e-07 3.26e-07
#ACDC 3           4 2098 31.9 1.20e-07 1.20e-07
### note: this is when using all 176 bats with insect outgroups!!!


## test raw codes
load("../topo11_mcmc/sif176_tree23_topo11_diet5.RData")
## conservative coding
for(sp in sif[sif$conservative_coding != "NO", "vertlife_name"]){
  sif[sif$vertlife_name == sp, c("Arthropods", "Fruit")] = c(3, 0) 
  if(sp == "Micronycteris_microtis"){sif[sif$vertlife_name == sp, "Arthropods"] = 2 }
}
rownames(sif) = sif$vertlife_name
## note: using raw codes of all 5 diets is super slow!!!
data = sif[, c("Arthropods", "Fruit", "Pollen.and.nectar")]


## all sif176
BM176 <- mvBM(tree, data, model="BM1")
EB176 <- mvEB(tree, data)
OU176 <- mvOU(tree, data, model="OU1")
results176 <- list(BM176, EB176, OU176)
aicw(results176)
#-- Akaike weights -- 
#                Rank AIC diff wi AICw
#BM1 default 1    1  915    0  1    1
#OU1 3            2 1120  204  0    0
#ACDC 2           3 1457  542  0    0

## only phyllostomids
tree = keep.tip(tree, tip=sif[sif$Family == "Phyllostomidae", "vertlife_name"])
data = data[sif[sif$Family == "Phyllostomidae", "vertlife_name"],]
BM138 <- mvBM(tree, data, model="BM1")
EB138 <- mvEB(tree, data)
OU138 <- mvOU(tree, data, model="OU1")
results138 <- list(BM138, EB138, OU138)
aicw(results138)
#-- Akaike weights -- 
#               Rank AIC diff wi AICw
#OU1 3            1 609    0  1    1
#BM1 default 1    2 714  105  0    0
#ACDC 2           3 716  107  0    0


## only after early-burst
tree = keep.tip(tree, tip=sif[sif$Family == "Phyllostomidae" & !(sif$Subfamily %in% c("Macrotinae", "Micronycterinae")), "vertlife_name"])
data = data[tree$tip.label,]
BM128 <- mvBM(tree, data, model="BM1")
EB128 <- mvEB(tree, data)
OU128 <- mvOU(tree, data, model="OU1")
results128 <- list(BM128, EB128, OU128)
aicw(results128)
#-- Akaike weights -- 
#Rank AIC diff wi AICw
#OU1 3            1 565  0.0  1    1
#BM1 default 1    2 650 85.7  0    0
#ACDC 2           3 652 87.7  0    0

data = scores[tree$tip.label, c("PC1", "PC2")]
BM <- mvBM(tree, data, model="BM1")
EB <- mvEB(tree, data)
OU <- mvOU(tree, data, model="OU1")
results <- list(BM, EB, OU)
aicw(results)
#-- Akaike weights -- 
#Rank  AIC  diff    wi  AICw
#OU1 3            1 1508 0.000 1.000 0.538
#BM1 default 1    2 1509 0.934 0.627 0.338
#ACDC 2           3 1511 2.934 0.231 0.124


data = sif[tree$tip.label, c("Fruit", "Pollen.and.nectar")]
BM128 <- mvBM(tree, data, model="BM1")
EB128 <- mvEB(tree, data)
OU128 <- mvOU(tree, data, model="OU1")
results128 <- list(BM128, EB128, OU128)
aicw(results128)
#-- Akaike weights -- 
#Rank AIC diff     wi   AICw
#OU1 3            1 423 0.00 1.0000 0.8132
#BM1 default 1    2 426 3.57 0.1679 0.1366
#ACDC 2           3 428 5.57 0.0618 0.0502

#################################################################
#################################################################
