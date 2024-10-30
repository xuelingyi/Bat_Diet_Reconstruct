## get libraries, functions
## load diets, diet.states, state.color
source("../my.diet.functions.R")

## download the full node-dated mammal tree from https://data.vertlife.org/
#tree.files = read.table("mammal.trees.all")

#vertlife = read.csv("Upham2019_NCBI.csv")
# all 5911 species
sif = read.csv("sif179/sif179.csv")
## TAXA!!!
#sif[sif$ScientificName %in% c("Mimon crenulatum", "Tonatia saurophila", "Lophostoma evotis"), c("ScientificName", diets, "source")]
#bat624 = read.csv("bat624/vertlife_bat624.csv")
#bat624[bat624$ScientificName %in% c("Gardnerycteris crenulata", "Tonatia bakeri", "Lophostoma evote"), c("ScientificName", diets, "diet_source")]

### below was not included in sif179 but SHOULD!!! (slight diet difference on Elton; should be fine)
sif[sif$ScientificName == "Mimon crenulatum", c("ScientificName", "Genus")] = c("Gardnerycteris crenulata", "Gardnerycteris")
sif[sif$ScientificName == "Tonatia saurophila", c("ScientificName")] = c("Tonatia bakeri")
sif[sif$ScientificName == "Lophostoma evotis", c("ScientificName")] = c("Lophostoma evote")

unique(sif[sif$Family == "Phyllostomidae", "Genus"])
# 55

write.csv(sif, "sif179/sif179.csv", row.names = F, quote = F)

subfamilies = unique(sif[sif$Family == "Phyllostomidae", "Subfamily"])

## target taxa
ingroup.f = c("Phyllostomidae", "Mormoopidae", "Noctilionidae", "Furipteridae", "Thyropteridae", "Mystacinidae") 
outgroup.f = c("Vespertilionidae", "Molossidae", "Emballonuridae")

##################################################################
######################### RUNONCE: get outgroups 26 #############################
## get outgroups from sif63
sif63 = read.csv("../../genomes/sif63/sif63.csv",row.names = 1)
outgroup = unique(sif63[!(sif63$Family %in% ingroup.f), "ScientificName"])

outgroup = gsub(" ", "_", outgroup)
# add the myotis that might be fish-eating; Eptesicus fuscus included
outgroup = c(outgroup, "Myotis_capaccinii", "Myotis_macropus")
# remove the two species not in vertlife
outgroup = outgroup[!(outgroup %in% c("Cnephaeus_nilssonii", "Molossus_nigricans"))]
# 26 species 
rm(sif63)

write.table(outgroup, "outgroup26", quote = F, row.names = F, col.names = F)

##################################################################
####################### RUNONCE taxa selection sif179 ##################################
## check vertlife inputs
vertlife = read.csv("Upham2019_NCBI.csv")
vertlife$gene_list = apply(vertlife, 1, FUN=function(x){
  paste0((names(vertlife)[12:42])[!is.na(x[12:42])], collapse = "; ")
})

## outgroups
outgroup = read.table("outgroup26")
vertlife_outgroup = vertlife[vertlife$Species_Name %in% outgroup$V1,]
# 26
vertlife_outgroup[vertlife_outgroup$genes <=1, c(1, 43)]
#Species_Name gene_list
#3119 Molossus_alvarezi          
#3488     Myotis_vivesi      CYTB
vertlife_outgroup = vertlife_outgroup[vertlife_outgroup$genes >1,]
# 24 species


vertlife_ingroup = vertlife[vertlife$fam %in% toupper(ingroup.f),]
#225
length(unique(vertlife_ingroup$gen))
#66
nrow(vertlife_ingroup[vertlife_ingroup$genes <=1, c(1,9,43)])
# 70
length(unique(vertlife_ingroup[vertlife_ingroup$genes >1, "gen"]))
#59
## total 7 genera only had one species and had only one gene, thus were removed
#all = unique(vertlife_ingroup$gen)
#vertlife_ingroup[vertlife_ingroup$gen %in% all[!(all %in% unique(vertlife_ingroup[vertlife_ingroup$genes >1, "gen"]))], c("Species_Name", "fam")]
#"Amorphochilus_schnablii" "Dryadonycteris_capixaba" "Lichonycteris_obscura"   "Neonycteris_pusilla"    
#[5] "Platalina_genovensium"   "Scleronycteris_ega"      "Xeronycteris_vieirai" 
#rm(all)
vertlife_ingroup = vertlife_ingroup[vertlife_ingroup$genes >1, ]
# 155
table(vertlife_ingroup$fam)
#FURIPTERIDAE    MORMOOPIDAE   MYSTACINIDAE  NOCTILIONIDAE PHYLLOSTOMIDAE  THYROPTERIDAE
#1              8              1              2            141              2 


vertlife_sif179 = rbind(vertlife_ingroup, vertlife_outgroup)
vertlife_sif179 = vertlife_sif179[, c("Species_Name", "fam", "gen", "genes", "gene_list")]
#rm(vertlife, vertlife_ingroup, vertlife_outgroup)
names(vertlife_sif179)[1] = "vertlife_name"
vertlife_sif179$ScientificName = gsub("_", " ", vertlife_sif179$vertlife_name)

##correct species names for merging with diet data
vertlife_sif179[vertlife_sif179$ScientificName == "Diaemus youngi", "ScientificName"] = "Diaemus youngii"
vertlife_sif179[vertlife_sif179$ScientificName == "Lophostoma silvicolum", "ScientificName"] = "Lophostoma silvicola"
vertlife_sif179[vertlife_sif179$ScientificName == "Dermanura cinereus", "ScientificName"] = "Dermanura cinerea"
vertlife_sif179[vertlife_sif179$ScientificName == "Dermanura glaucus", "ScientificName"] = "Dermanura glauca"
vertlife_sif179[vertlife_sif179$ScientificName == "Dermanura toltecus", "ScientificName"] = "Dermanura tolteca"
vertlife_sif179[vertlife_sif179$ScientificName == "Vampyressa bidens", c("ScientificName", "gen")] = c("Vampyriscus bidens", "Vampyriscus")
vertlife_sif179[vertlife_sif179$ScientificName == "Vampyressa brocki", c("ScientificName", "gen")] = c("Vampyriscus brocki", "Vampyriscus")
vertlife_sif179[vertlife_sif179$ScientificName == "Vampyressa nymphaea", c("ScientificName", "gen")] = c("Vampyriscus nymphaeus", "Vampyriscus")


## add diet
rojas = read.csv("../../Rojas2018_updated.csv")
vertlife_sif179[!(vertlife_sif179$ScientificName %in% rojas$Taxon), "ScientificName"]
# none 
vertlife_sif179 = merge(vertlife_sif179, rojas, by.x = "ScientificName", by.y="Taxon", all.x = T)
vertlife_sif179$fam = str_to_sentence(vertlife_sif179$fam)
vertlife_sif179$gen = str_to_sentence(vertlife_sif179$gen)
## double-check taxonomy names
all(vertlife_sif179$fam == vertlife_sif179$Family, na.rm = T)
#TRUE
all(vertlife_sif179$gen == vertlife_sif179$Genus, na.rm = T)
#TRUE
vertlife_sif179 = vertlife_sif179[, c("ScientificName", "vertlife_name", "fam", "Subfamily", "gen", "genes", "gene_list", "Arthropods", "Blood", "Terrestrial.vertebrates", "Fish", "Leaves.and.flower.pieces", "Pollen.and.nectar", "Fruit", "Seed", "Update_note")]
names(vertlife_sif179)[3:5] = c("Family", "Subfamily", "Genus")


vertlife_sif179[vertlife_sif179$Family == "Phyllostomidae" & is.na(vertlife_sif179$Subfamily), "vertlife_name"]
# all have subfamily names

table(vertlife_sif179$Family)
#Emballonuridae     Furipteridae       Molossidae      Mormoopidae     Mystacinidae    Noctilionidae 
#4                1                3                8                1                2 
#Phyllostomidae    Thyropteridae Vespertilionidae 
#141                2               17


system("mkdir sif179")
write.csv(vertlife_sif179, "sif179/vertlife_sif179.csv", quote = F, row.names = F)

##################################################################
####### RUNONCE test monophyly & summarize topology #################
## download the full node-dated mammal tree from https://data.vertlife.org/
tree.files = read.table("mammal.trees.all")

## check monophyly of genus, subfamily, and family
mono_test = as.data.frame(matrix(nrow=length(c(unique(sif$Genus), subfamilies, ingroup.f, outgroup.f)), ncol=(10000 +2)))
names(mono_test) = c("level", "taxa", paste0("tree_", 1:10000))
mono_test$taxa = c(unique(sif$Genus), subfamilies, ingroup.f, outgroup.f)
mono_test[mono_test$taxa %in% sif$Genus, "level"] = "Genus"
mono_test[mono_test$taxa %in% sif$Subfamily, "level"] = "Subfamily"
mono_test[mono_test$taxa %in% sif$Family, "level"] = "Family"

table(mono_test$level)
#Family     Genus Subfamily 
#9        74        11 

## representatives of subfamilies and families (except for Mormoopidae which is not monophyletic and thus had reps for each genus)
topo_reps = c("Saccopteryx_leptura", "Eptesicus_fuscus", "Molossus_molossus",
               "Mystacina_tuberculata", 
               "Furipterus_horrens", "Mormoops_megalophylla", "Pteronotus_parnellii", "Pteronotus_davyi", "Noctilio_leporinus", "Thyroptera_tricolor",  
               "Carollia_perspicillata", "Desmodus_rotundus", "Glossophaga_soricina", "Glyphonycteris_daviesi", "Lionycteris_spurrelli", "Macrotus_waterhousii", "Micronycteris_megalotis", "Phyllostomus_discolor", "Rhinophylla_pumilio", "Sturnira_hondurensis", "Lonchorhina_aurita")


subf_trees = NULL
pg_trees = NULL
#pdf("sif179_monophyly.pdf")
for(i in 1:10000){
  tree = read.tree(paste0("~/Desktop/DNAonly_4098sp_topoFree_NDexp/", tree.files[i,]))
  tree = keep.tip(tree, sif$vertlife_name)
  
  # test monophyly
  for(t in mono_test$taxa) {
    my.tips = sif[sif$Family == t | sif$Subfamily == t | sif$Genus == t, "vertlife_name"]
    
    if(!(is.monophyletic(tree, tips=my.tips))){
      mono_test[mono_test$taxa == t, paste0("tree_",i)] = "not_monop"
      
      #only print if not monophyly on sub- or family levels
  #    if(mono_test[mono_test$taxa == t, "level"] != "Genus"){
  #      print(ggtree(groupOTU(tree, my.tips), aes(color=group)) + 
  #              geom_tiplab(size=1) + 
  #              scale_color_manual(values=c("black", "red"), name=paste0(t,"_tree_", i)) + 
  #              hexpand(0.1,1))
  #    } 
      
   }
  }
  
  #summarize topology
  #subf_trees = append(subf_trees, keep.tip(tree, topo_reps)) 
  #pg_trees = append(pg_trees, keep.tip(tree, topo_reps[topo_reps %in% sif[sif$Family == "Phyllostomidae", "vertlife_name"]])) 
}
#dev.off()
#dev.off()


mono_test$count = apply(mono_test, 1, FUN=function(x){
  sum(!is.na(x[3:10002]))
})
mono_test[mono_test$count != 0, c(1,2,10003)]
#level           taxa count
#11     Genus   Choeroniscus    22
#25     Genus Glyphonycteris     1
#31     Genus   Lonchophylla 10000
#38     Genus          Mimon 10000
#51     Genus   Phyllostomus 10000
#52     Genus   Pipistrellus  9171
#66     Genus        Tonatia     2
#70     Genus     Vampyressa     2
#79 Subfamily  Desmodontinae     1
#87    Family    Mormoopidae  9895


unique(sif[sif$Genus %in% mono_test[mono_test$count != 0, "taxa"], c("Family", "Subfamily", "Genus")])
#.      Family        Subfamily          Genus
#29    Phyllostomidae   Glossophaginae   Choeroniscus
#55    Phyllostomidae Glyphonycterinae Glyphonycteris
#64    Phyllostomidae  Lonchophyllinae   Lonchophylla
#86    Phyllostomidae   Phyllostominae          Mimon
#109   Phyllostomidae   Phyllostominae   Phyllostomus
#112 Vespertilionidae                    Pipistrellus
#165   Phyllostomidae   Phyllostominae        Tonatia
#171   Phyllostomidae  Stenodermatinae     Vampyressa


## representative topologies
#ggdensitree(subf_trees, alpha=.3, colour='steelblue') + 
#  geom_tiplab(size=1) + 
#  hexpand(.25) 
#for(i in 1:20){
#  print(ggtree(subf_trees[[i]]) + geom_tiplab(aes(label=gsub("_.*", "", label)), size=2) + hexpand(.25))
#}
all.topology = unique.multiPhylo(subf_trees, use.edge.length = FALSE, use.tip.label = TRUE)
#162 phylogenetic trees
all.topology.summary = data.frame(matrix(nrow=length(all.topology), ncol=3))
names(all.topology.summary) = c("topology", "trees", "count")
all.topology.summary$topology = 1:length(all.topology)
for(j in 1:nrow(all.topology.summary)){
  all.trees = NULL
  count.trees = 0
  
  for(i in 1:10000){
    if (all.equal.phylo(subf_trees[[i]], all.topology[[j]], use.edge.length = FALSE)) {
      all.trees = paste0(all.trees, "; tree", i)
      count.trees = count.trees + 1
    }
  }
  
  all.topology.summary[j, "trees"] = all.trees
  all.topology.summary[j, "count"] = count.trees
}
#ggplot(all.topology.summary) + geom_histogram(aes(x=count))
range(all.topology.summary$count)
#1 1883
sum(all.topology.summary$count)
#10000
all.topology.summary$percent = all.topology.summary$count / 10000
range(all.topology.summary$percent)
#[1] 0.0001 0.1883
all.topology.summary[all.topology.summary$percent > 0.05, c("topology", "count")]
#   topology count
#2         2   558
#3         3   872
#4         4  1172
#5         5  1883
#23       23   746
pdf("topology162.pdf")
for(i in 1:length(all.topology)){
  print(ggtree(groupOTU(all.topology[[i]], 
                        list(Glossophaginae = "Glossophaga_soricina", 
                             Phyllostominae="Phyllostomus_discolor", 
                             Lonchorhininae="Lonchorhina_aurita")), 
               aes(color=group),
               branch.length="none") + 
          scale_color_manual(values=c("black", "orange", "royalblue3","purple4"),
                             name=paste0(all.topology.summary[i, "percent"]*100, "% trees")) +
          geom_tiplab(aes(label=gsub("_.*", "", label))) +
          hexpand(0.2,1))
}
dev.off()



pg.topology = unique.multiPhylo(pg_trees, use.edge.length = FALSE, use.tip.label = TRUE)
#41
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
sum(as.numeric(pg.topology.summary$count))
pdf("topology_pg41.pdf")
for(i in 1:length(pg.topology)){
  print(ggtree(groupOTU(pg.topology[[i]], 
                        list(Glossophaginae = "Glossophaga_soricina", 
                             Phyllostominae="Phyllostomus_discolor", 
                             Lonchorhininae="Lonchorhina_aurita")), 
               aes(color=group),
               branch.length="none") + 
          scale_color_manual(values=c("black", "orange", "royalblue3","purple4"),
                             name=paste0(as.numeric(pg.topology.summary[i, "count"])/100, "% trees")) +
          geom_tiplab(aes(label=gsub("_.*", "", label))) +
          hexpand(0.2,1))
}
dev.off()

## summarize who diverged first
L = c(3,4,10,19,23,27,31,36,37)
P = c(9,11,26,29,33,40,41)
LP = c(6,14,32)
LG = c(8,13,20,22,28)
GP = c(12)
G = (1:41)[!((1:41) %in% L)]
G = G[!(G %in% c(P, LP, LG, GP))]
Lionycteris = c(13, 14, 15,16,18,19,20,21,22,23,27,28,29,30,31,32,34,35,40,41)

pg.topology.summary$count = as.numeric(pg.topology.summary$count)
## two sumf diverged first?
sum(pg.topology.summary[LG, "count"])/100
#[1] 1.93
sum(pg.topology.summary[LP, "count"])/100
#[1] 7.7
sum(pg.topology.summary[GP, "count"])/100
#[1] 0.03

## most had Glossophaginae diverged first
sum(pg.topology.summary[G, "count"])/100
#[1] 75.1
# Lonchorhininae diverged first
sum(pg.topology.summary[L, "count"])/100
#[1] 14.59
# Phyllostominae first
sum(pg.topology.summary[P, "count"])/100
#[1] 0.65
sum(pg.topology.summary[Lionycteris, "count"])/100
#[1] 2.8

pg.topology.summary$topology = row.names(pg.topology.summary)

save(subf_trees, all.topology, all.topology.summary, 
     pg_trees, pg.topology, pg.topology.summary,
     mono_test, 
     file="sif179/sif179_10k_trees_monophyly.RData")
## moved to all_10k



load("sif179/all_10k/sif179_10k_trees_monophyly.RData")

pg.topology.summary$count = as.numeric(pg.topology.summary$count)
pg.topology.summary = pg.topology.summary[order(pg.topology.summary$count, decreasing = T),]
sum(pg.topology.summary[1:5, "count"])
#9156
pg.topology.top5 = pg.topology.summary[1:5, "topology"]
#[1] "2"  "1"  "10" "5"  "6" 

##################################################################
############################ didNOTrun compare with Rojas 2016 ##############################
#compare topology with Rojas 2016
i=10000
tree = read.tree(paste0("~/Desktop/DNAonly_4098sp_topoFree_NDexp/", tree.files[i,]))
tree = keep.tip(tree, sif211$vertlife_name)
rojas2016 = read.nexus("../../paper_supplementary/Rojas2016_SystBiol_doi_10.5061_dryad.s533p__v2/MLtree_1_1430845747_Noctilionoidea_ML.tree")

rojas_tips = rojas2016$tip.label
sum(gsub(" ", "_", sif211$Species) %in% rojas_tips)
# 149
tips.in = rojas_tips[rojas_tips %in% sif211$vertlife_name]
# 150 shared tips
pdf("comp_tree9999_rojas2016.pdf")
plot(cophylo(tree, rojas2016, assoc = cbind(tips.in, tips.in)),
     link.lwd=1, link.type="curved",
     link.lty="solid", fsize=c(0.3,0.3))
dev.off()

##################################################################



######################## RUNONCE monophyly of internal nodes ####################################
## download the full node-dated mammal tree from https://data.vertlife.org/
tree.files = read.table("mammal.trees.all")

## check monophyly of five key internal nodes + 2 subf nodes
my.nodes = c("P_M", "super_Noctilionoidea", "N_M", "V_M", "V_M_E", "SR", "CG", "SRCG_fruit", "early_burst", "eb_mi", "eb_mi_de")

mono_internal_test = as.data.frame(matrix(nrow=length(my.nodes), ncol=(10000 +1)))
names(mono_internal_test) = c("node", paste0("tree_", 1:10000))
mono_internal_test$node = my.nodes

for(i in 1:10000){
  tree = read.tree(paste0("~/Desktop/DNAonly_4098sp_topoFree_NDexp/", tree.files[i,]))
  tree = keep.tip(tree, sif$vertlife_name)
  
  # test monophyly of internal nodes (this can be added above)
  for(t in mono_internal_test$node) {
    if(t == "P_M"){ my.tips = sif[sif$Family %in% c("Phyllostomidae", "Mormoopidae"), "vertlife_name"] }
    if(t == "super_Noctilionoidea"){ 
      my.tips = sif[sif$Family %in% c("Phyllostomidae", "Mormoopidae", "Noctilionidae", "Furipteridae", "Thyropteridae"), "vertlife_name"] }
    if(t == "N_M"){ my.tips = sif[sif$Family %in% ingroup.f, "vertlife_name"] }
    if(t == "V_M"){ my.tips = sif[sif$Family %in% c("Vespertilionidae", "Molossidae"), "vertlife_name"] }
    if(t == "V_M_E"){ my.tips = sif[sif$Family %in% outgroup.f, "vertlife_name"] }
    if(t == "SR"){ my.tips = sif[sif$Subfamily %in% c("Stenodermatinae", "Rhinophyllinae"), "vertlife_name"] }
    if(t == "CG"){ my.tips = sif[sif$Subfamily %in% c("Carolliinae", "Glyphonycterinae"), "vertlife_name"] }
    if(t == "SRCG_fruit"){ my.tips = sif[sif$Subfamily %in% c("Stenodermatinae", "Rhinophyllinae", "Carolliinae", "Glyphonycterinae"), "vertlife_name"] }
    if(t == "early_burst"){ my.tips = sif[sif$Family == "Phyllostomidae" & !(sif$Subfamily %in% c("Desmodontinae", "Micronycterinae", "Macrotinae")), "vertlife_name"] }
    if(t == "eb_mi"){ my.tips = sif[sif$Family == "Phyllostomidae" & !(sif$Subfamily %in% c("Desmodontinae", "Macrotinae")), "vertlife_name"] }
    if(t == "eb_mi_de"){ my.tips = sif[sif$Family == "Phyllostomidae" & sif$Subfamily != "Macrotinae", "vertlife_name"] }
    
    
    if(!(is.monophyletic(tree, tips=my.tips))){
      mono_internal_test[mono_internal_test$node == t, paste0("tree_",i)] = "not_monop"
    }
  }
}

mono_internal_test$count = apply(mono_internal_test, 1, FUN=function(x){
  sum(!is.na(x[2:10001]))
})
mono_internal_test[mono_internal_test$count != 0, c(1,10002)]
# node count
#5     V_M_E  1168
#7        CG     5
#11 eb_mi_de    15

#1168/10000*100
#[1] 11.68

## check topologies regarding the non-mono outgroups
vme = mono_internal_test[mono_internal_test$node == "V_M_E", 2:10001]
treeID = gsub("tree_", "", names(vme[,!is.na(vme)]))
## if Emballonuridae is sister to NM
E_NM = NULL
## if Emballonuridae is the basal
E_base = NULL
## if else
E_other = NULL
for(i in treeID){
  tree = read.tree(paste0("~/Desktop/DNAonly_4098sp_topoFree_NDexp/", tree.files[i,]))
  tree = keep.tip(tree, sif$vertlife_name)
  
  # test monophyly of internal nodes (this can be added above)
  if(!(is.monophyletic(tree, 
                       tips=sif[sif$Family %in% c("Emballonuridae", ingroup.f), "vertlife_name"]))){
    E_NM = c(E_NM, paste0("tree_", i))
  } else {
    if(!(is.monophyletic(tree, 
                         tips=sif[sif$Family != "Emballonuridae", "vertlife_name"]))){
      E_base = c(E_base, paste0("tree_", i))
    } else {
      E_other = c(E_other, paste0("tree_", i))
    }
  } 
}
length(treeID)
length(E_NM)
#[1] 12
length(E_base)
#[1] 1156
length(E_other)
#[1] 0


save(mono_internal_test, E_base, E_NM,
     file="sif179/all_10k/sif179_10k_trees_monophyly_internal.RData")


##################################################################
######################## RUNONCE monophyly of CORRECTED genus ####################################
tree.files = read.table("mammal.trees.all")

## check monophyly of genus
my.nodes = unique(sif$Genus)

mono_internal_test = as.data.frame(matrix(nrow=length(my.nodes), ncol=(10000 +1)))
names(mono_internal_test) = c("node", paste0("tree_", 1:10000))
mono_internal_test$node = my.nodes

for(i in 1:10000){
  tree = read.tree(paste0("~/Desktop/DNAonly_4098sp_topoFree_NDexp/", tree.files[i,]))
  tree = keep.tip(tree, sif$vertlife_name)
  
  # test monophyly of internal nodes (this can be added above)
  for(t in mono_internal_test$node) {
    my.tips = sif[sif$Genus == t, "vertlife_name"] 

    if(!(is.monophyletic(tree, tips=my.tips))){
      mono_internal_test[mono_internal_test$node == t, paste0("tree_",i)] = "not_monop"
    }
  }
}

mono_internal_test$count = apply(mono_internal_test, 1, FUN=function(x){
  sum(!is.na(x[2:10001]))
})
mono_internal_test[mono_internal_test$count != 0, c(1,10002)]
#           node count
#11   Choeroniscus    22
#25 Glyphonycteris     1
#31   Lonchophylla 10000
#52   Phyllostomus 10000
#53   Pipistrellus  9171 --> the only one not phyllostomids
#67        Tonatia     2
#71     Vampyressa     2

unique(sif[sif$Genus %in% mono_internal_test[mono_internal_test$count != 0, 1], c("Family", "Genus")])

##################################################################