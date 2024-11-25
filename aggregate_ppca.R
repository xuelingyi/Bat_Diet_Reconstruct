## get libraries, functions
## load diets, diet.states, state.color
source("my.diet.functions.R")

load("sif176_tree23_topo11_diet5.RData")
tree_list = read.table("tree23_batch2.txt")

## conservative coding
for(sp in sif[sif$conservative_coding != "NO", "vertlife_name"]){
  sif[sif$vertlife_name == sp, c("Arthropods", "Fruit")] = c(3, 0) 
  if(sp == "Micronycteris_microtis"){sif[sif$vertlife_name == sp, "Arthropods"] = 2 }
}

## main diets
diet3 = c("Arthropods", "Fruit", "Pollen.and.nectar")


######### aggregate mcmcglmm results 
i=1
id = tree_list[i, "V1"]
estimates_all = read.csv(paste0("trees_conserve/tree_", id, "/estimates.csv"))

for(i in 2:23){
  id = tree_list[i, "V1"]
  estimates = read.csv(paste0("trees_conserve/tree_", id, "/estimates.csv"))
  estimates_all = rbind(estimates_all, estimates)
}
# write.csv(estimates_all, "estimates.all23.csv", quote = F, row.names = F)



######### plot original and reconstructed pPCA
taxa = mono.nodes[!(mono.nodes %in% c("EB_Mic", "VM"))]

## mean scores and loadings of pCPA from 23 trees
load(paste0("inputs/sif176_ppca_topo11_index_", 1, ".RData"))
loadings_mean = as.data.frame(unique(loadings$diet)) 
names(loadings_mean) = "diet"
loadings_mean$PC1 = 0
loadings_mean$PC2 = 0
Eval_PC1 = 0
Eval_PC2 = 0

for(i in 1:23){
  load(paste0("inputs/sif176_ppca_topo11_index_", i, ".RData"))
  for(d in loadings_mean$diet){
    loadings_mean[loadings_mean$diet == d, "PC1"] = sum(loadings_mean[loadings_mean$diet == d, "PC1"], loadings[loadings$diet == d, "PC1"])
    loadings_mean[loadings_mean$diet == d, "PC2"] = sum(loadings_mean[loadings_mean$diet == d, "PC2"], loadings[loadings$diet == d, "PC2"])
  }
  
  Eval_PC1 = sum(Eval_PC1, ppca$Eval[1,1]/sum(ppca$Eval)*100)
  Eval_PC2 = sum(Eval_PC2, ppca$Eval[2,2]/sum(ppca$Eval)*100)
}
loadings_mean$PC1 = loadings_mean$PC1 / 23
loadings_mean$PC2 = loadings_mean$PC2 / 23
Eval_PC1 = Eval_PC1 / 23
Eval_PC2 = Eval_PC2 / 23


## mean posterior PC scores across 23 trees
# estimates_all = read.csv("estimates.all23.csv")
estimates_all = merge(estimates_all, sif[, c(diet3, "vertlife_name")], all.x = T)

estimates_mean = as.data.frame(unique(estimates_all[, c("vertlife_name", "type")]))
estimates_mean = estimates_mean[estimates_mean$vertlife_name %in% c(taxa, sif$vertlife_name),]
estimates_mean$PC1 = 0
estimates_mean$PC2 = 0
for(i in 1:nrow(estimates_mean)){
  sp = estimates_mean[i, "vertlife_name"]
  type = estimates_mean[i, "type"]
  estimates_mean[i, "PC1"] = mean(estimates_all[estimates_all$vertlife_name == sp & estimates_all$type == type, "PC1"])
  estimates_mean[i, "PC2"] = mean(estimates_all[estimates_all$vertlife_name == sp & estimates_all$type == type, "PC2"])
  if(nrow(estimates_all[estimates_all$vertlife_name == sp & estimates_all$type == type, ]) != 23){print(paste0(i, ":nrows not 23!!"))}
}

## node input states
sif$AFN = paste0(sif$Arthropods, sif$Fruit, sif$Pollen.and.nectar)
estimates_mean = merge(estimates_mean, sif[, c("vertlife_name", "AFN")], by="vertlife_name", all.x=T)
p0 = ggplot() + 
  geom_segment(data=loadings_mean, 
               aes(x=0, xend=PC1*5, y=0, yend=PC2*5, color=diet), arrow = arrow(length = unit(0.03, "npc"))) +
  geom_text(data=loadings_mean, aes(x=PC1*6, y=PC2*6, label = diet, color=diet), size=4) +
  scale_color_manual(values=c("Arthropods"="black", "Fruit"="#2B823A", "Pollen.and.nectar"="orange"), guide="none") +
  scale_fill_manual(values=c("input"="black", "output"="grey75"), name="mean_23trees") +
  geom_line(data=estimates_mean[estimates_mean$vertlife_name %in% sif$vertlife_name,], aes(x=PC1, y=PC2, group=vertlife_name), color="grey50", alpha=0.5) +
  labs(x=paste0("PC1: ", round(Eval_PC1, 2), "%"),
       y=paste0("PC2: ", round(Eval_PC2, 2), "%")) +
  theme_bw()


ggarrange(nrow = 1, ncol = 2, 
          p0 + geom_point(data=estimates_mean[estimates_mean$vertlife_name %in% sif$vertlife_name & estimates_mean$type == "output",], aes(x=PC1, y=PC2, fill=type), size=3, alpha=0.7, shape=21, color="grey75") + 
            geom_text(data=estimates_mean[estimates_mean$vertlife_name %in% sif$vertlife_name & estimates_mean$type == "input",], aes(x=PC1, y=PC2, label=AFN), size=4.5), 
          p0 + geom_point(data=estimates_mean[estimates_mean$vertlife_name %in% sif$vertlife_name,], aes(x=PC1, y=PC2, fill=type), size=3, alpha=0.7, shape=21, color="grey75") +
            geom_point(data=estimates_mean[estimates_mean$vertlife_name %in% taxa,], aes(x=PC1, y=PC2), size=3, color="royalblue3") +
            geom_text(data=estimates_mean[estimates_mean$vertlife_name %in% taxa,], aes(x=PC1, y=PC2, label=vertlife_name), size=4.5, color="royalblue3") )

