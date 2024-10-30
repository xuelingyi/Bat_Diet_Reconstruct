## get libraries, functions
## load diets, diet.states, state.color
source("../my.diet.functions.R")

sif = read.csv("bat624/vertlife_bat624.csv")
# 21 families
table(sif$Family)
#Cistugidae Craseonycteridae   Emballonuridae     Furipteridae   Hipposideridae   Megadermatidae 
#1                1               24                1               40                3 
#Miniopteridae       Molossidae      Mormoopidae     Mystacinidae      Myzopodidae        Natalidae 
#11               41                8                1                1                8 
#Noctilionidae       Nycteridae   Phyllostomidae     Pteropodidae    Rhinolophidae  Rhinonycteridae 
#2                4              141               85               32                6 
#Rhinopomatidae    Thyropteridae Vespertilionidae 
#2                2              210

unique(sif[sif$Clade == "Yinpterochiroptera", "Family"])
# 7
unique(sif[sif$Clade == "Yangochiroptera", "Family"])
#14
##################################################################
########################## RUNONCE: get all bats having >1 gene #####################################
vertlife = read.csv("Upham2019_NCBI.csv")
# all 5911 species
vertlife = vertlife[vertlife$ord == "CHIROPTERA",]
# all 1287 bats; 18 families
range(vertlife$genes)
# 0, 30
vertlife$ord = str_to_sentence(vertlife$ord)
vertlife$fam = str_to_sentence(vertlife$fam)
## add subfamilies
diet = read.csv("../diet_data.csv")
vertlife$subf = sapply(vertlife$gen, FUN = function(x){
  if(x %in% diet[diet$Family == "Phyllostomidae", "Genus"]){
    return(unique(diet[diet$Genus == x, "Subfamily"]))
  } else {
    return("none")
  }
})

sif = vertlife[vertlife$genes > 1,]
# 624 bats, 18 families 
sif= sif[, c("Species_Name", "ord", "clade", "fam", "subf", "gen", "genes")]
names(sif) = c("vertlife_name", "Order", "Clade", "Family", "Subfamily", "Genus", "genes")

#write.csv(sif, "bat624/vertlife_bat624.csv", row.names = F, quote = F)
table(sif$Family)
#Craseonycteridae   Emballonuridae     Furipteridae   Hipposideridae   Megadermatidae       Molossidae 
#1               24                1               46                3               41 
#Mormoopidae     Mystacinidae      Myzopodidae        Natalidae    Noctilionidae       Nycteridae 
#8                1                1                8                2                4 
#Phyllostomidae     Pteropodidae    Rhinolophidae   Rhinopomatidae    Thyropteridae Vespertilionidae 
#141               85               32                2                2              222 


## check what was left out (gene==1)
sum(vertlife$genes == 0)
# 434
data.exc = vertlife[vertlife$genes == 1,]
# 229
## these genera were not included
unique(data.exc[!(data.exc$gen %in% sif$Genus), c("fam", "gen")])
#                 fam            gen
#640    Megadermatidae    Cardioderma
#674    Emballonuridae Centronycteris
#793        Molossidae    Cheiromeles
#925    Emballonuridae        Cormura
#2047 Vespertilionidae  Hesperoptenus
#2247 Vespertilionidae        Hypsugo
#2270 Vespertilionidae   Idionycteris
#2519   Phyllostomidae  Lichonycteris ok
#3100     Pteropodidae       Mirimiri ok
#3765 Vespertilionidae    Nyctophilus
#4022     Pteropodidae  Paranyctimene ok
#4376   Phyllostomidae      Platalina ok
#4617     Pteropodidae     Pteralopex ok
#4997   Emballonuridae    Saccolaimus
#5372 Vespertilionidae    Submyotodon
#5881   Phyllostomidae   Xeronycteris ok

## how many species excluded per family
table(data.exc$fam)
#Emballonuridae   Hipposideridae   Megadermatidae       Molossidae      Myzopodidae       Nycteridae 
#16               15                1               19                1                2 
#Phyllostomidae     Pteropodidae    Rhinolophidae   Rhinopomatidae    Thyropteridae Vespertilionidae 
#32               37               26                2                1               77

##################################################################
########################### RUNONCE: check with batnames ###########################
## need to update species taxonomy!!
sif$ScientificName = gsub("_", " ", sif$vertlife_name)
##correct phyllostomid names 
sif[sif$ScientificName == "Diaemus youngi", "ScientificName"] = "Diaemus youngii"
sif[sif$ScientificName == "Lophostoma silvicolum", "ScientificName"] = "Lophostoma silvicola"
sif[sif$ScientificName == "Dermanura cinereus", "ScientificName"] = "Dermanura cinerea"
sif[sif$ScientificName == "Dermanura glaucus", "ScientificName"] = "Dermanura glauca"
sif[sif$ScientificName == "Dermanura toltecus", "ScientificName"] = "Dermanura tolteca"
sif[sif$ScientificName == "Vampyressa bidens", c("ScientificName", "Genus")] = c("Vampyriscus bidens", "Vampyriscus")
sif[sif$ScientificName == "Vampyressa brocki", c("ScientificName", "Genus")] = c("Vampyriscus brocki", "Vampyriscus")
sif[sif$ScientificName == "Vampyressa nymphaea", c("ScientificName", "Genus")] = c("Vampyriscus nymphaeus", "Vampyriscus")


### below was not included in sif179 but SHOULD!!!
sif[sif$ScientificName == "Mimon crenulatum", c("ScientificName", "Genus")] = c("Gardnerycteris crenulata", "Gardnerycteris")
sif[sif$ScientificName == "Tonatia saurophila", c("ScientificName")] = c("Tonatia bakeri")
sif[sif$ScientificName == "Lophostoma evotis", c("ScientificName")] = c("Lophostoma evote")



### below is additional for other bats
sif[sif$ScientificName == "Desmalopex leucopterus", c("ScientificName")] = c("Desmalopex leucoptera")
sif[sif$ScientificName == "Megaerops wetmorei", c("ScientificName", "Genus")] = c("Ptenochirus wetmorei", "Ptenochirus")
sif[sif$ScientificName == "Melonycteris fardoulisi", c("ScientificName", "Genus")] = c("Nesonycteris fardoulisi", "Nesonycteris")
sif[sif$ScientificName == "Melonycteris woodfordi", c("ScientificName", "Genus")] = c("Nesonycteris woodfordi", "Nesonycteris")
sif[sif$ScientificName == "Micropteropus pusillus", c("ScientificName", "Genus")] = c("Epomophorus pusillus", "Epomophorus")
sif[sif$ScientificName == "Pteropus giganteus", c("ScientificName")] = c("Pteropus vampyrus") # synonymy


sif[sif$ScientificName == "Penthetor lucasi", c("ScientificName")] = c("Penthetor lucasii")
sif[sif$ScientificName == "Ptenochirus jagori", c("ScientificName")] = c("Ptenochirus jagorii")
sif[sif$ScientificName == "Rousettus bidens", c("ScientificName", "Genus")] = c("Boneia bidens", "Boneia")
sif[sif$ScientificName == "Rousettus lanosus", c("ScientificName", "Genus")] = c("Stenonycteris lanosa", "Stenonycteris")
sif[sif$ScientificName == "Hipposideros commersoni", c("ScientificName", "Genus")] = c("Macronycteris commersonii", "Macronycteris")
sif[sif$ScientificName == "Hipposideros vittatus", c("ScientificName", "Genus")] = c("Macronycteris vittata", "Macronycteris")
sif[sif$ScientificName == "Hipposideros cyclops", c("ScientificName", "Genus")] = c("Doryrhina cyclops", "Doryrhina")
sif[sif$ScientificName == "Hipposideros pendelburyi", c("ScientificName")] = c("Hipposideros pendleburyi")
sif[sif$ScientificName == "Hipposideros pratti", c("ScientificName")] = c("Hipposideros swinhoii")
## synonymy on batnames
sif[sif$ScientificName == "Diclidurus isabellus", c("ScientificName")] = c("Diclidurus isabella")
sif[sif$ScientificName == "Emballonura atrata", c("ScientificName", "Genus")] = c("Paremballonura atrata", "Paremballonura")
sif[sif$ScientificName == "Emballonura tiavato", c("ScientificName", "Genus")] = c("Paremballonura tiavato", "Paremballonura") 
sif[sif$ScientificName == "Chaerephon atsinanana", c("ScientificName", "Genus")] = c("Mops atsinanana", "Mops")
sif[sif$ScientificName == "Chaerephon nigeriae", c("ScientificName", "Genus")] = c("Mops nigeriae", "Mops")
sif[sif$ScientificName == "Chaerephon plicatus", c("ScientificName", "Genus")] = c("Mops plicatus", "Mops")
sif[sif$ScientificName == "Chaerephon pumilus", c("ScientificName", "Genus")] = c("Mops pumilus", "Mops")
sif[sif$ScientificName == "Tadarida jobimena", c("ScientificName", "Genus")] = c("Mops jobimena", "Mops")
sif[sif$ScientificName == "Cynomops paranus", c("ScientificName")] = c("Cynomops planirostris") ## synonymy on batnames
sif[sif$ScientificName == "Molossops mattogrossensis", c("ScientificName", "Genus")] = c("Neoplatymops mattogrossensis", "Neoplatymops")
sif[sif$ScientificName == "Natalus saturatus", c("ScientificName")] = c("Natalus stramineus") # subspecies based on batnames
sif[sif$ScientificName == "Megaderma lyra", c("ScientificName", "Genus")] = c("Lyroderma lyra", "Lyroderma")

sif[sif$ScientificName == "Eptesicus bottae", c("ScientificName", "Genus")] = c("Cnephaeus bottae", "Cnephaeus")
sif[sif$ScientificName == "Eptesicus brasiliensis", c("ScientificName", "Genus")] = c("Neoeptesicus brasiliensis", "Neoeptesicus")
sif[sif$ScientificName == "Eptesicus diminutus", c("ScientificName", "Genus")] = c("Neoeptesicus diminutus", "Neoeptesicus")
sif[sif$ScientificName == "Eptesicus dimissus", c("ScientificName", "Genus")] = c("Cassistrellus dimissus", "Cassistrellus")
sif[sif$ScientificName == "Eptesicus furinalis", c("ScientificName", "Genus")] = c("Neoeptesicus furinalis", "Neoeptesicus")
sif[sif$ScientificName == "Eptesicus gobiensis", c("ScientificName", "Genus")] = c("Cnephaeus gobiensis", "Cnephaeus")
sif[sif$ScientificName == "Eptesicus hottentotus", c("ScientificName", "Genus")] = c("Cnephaeus hottentotus", "Cnephaeus")
sif[sif$ScientificName == "Eptesicus isabellinus", c("ScientificName", "Genus")] = c("Cnephaeus isabellinus", "Cnephaeus")
sif[sif$ScientificName == "Eptesicus nilssonii", c("ScientificName", "Genus")] = c("Cnephaeus nilssonii", "Cnephaeus")
sif[sif$ScientificName == "Eptesicus serotinus", c("ScientificName", "Genus")] = c("Cnephaeus serotinus", "Cnephaeus")
sif[sif$ScientificName == "Eptesicus pachyomus", c("ScientificName", "Genus")] = c("Cnephaeus pachyomus", "Cnephaeus")

sif[sif$ScientificName == "Murina cineracea", c("ScientificName")] = c("Murina feae") # batname synonymy
sif[sif$ScientificName == "Murina huttoni", c("ScientificName")] = c("Murina huttonii")
sif[sif$ScientificName == "Murina tiensa", c("ScientificName")] = c("Murina harrisoni") # batname synonymy
sif[sif$ScientificName == "Myotis melanorhinus", c("ScientificName")] = c("Myotis ciliolabrum") # batname synonymy
sif[sif$ScientificName == "Myotis aurascens", c("ScientificName")] = c("Myotis davidii") # batname synonymy
sif[sif$ScientificName == "Myotis annamiticus", c("ScientificName")] = c("Myotis laniger") # batname synonymy
sif[sif$ScientificName == "Neoromicia brunnea", c("ScientificName")] = c("Pseudoromicia brunnea") 
sif[sif$ScientificName == "Neoromicia rendalli", c("ScientificName")] = c("Pseudoromicia rendalli") 
sif[sif$ScientificName == "Neoromicia tenuipinnis", c("ScientificName", "Genus")] = c("Pseudoromicia tenuipinnis", "Pseudoromicia")
sif[sif$ScientificName == "Neoromicia capensis", c("ScientificName", "Genus")] = c("Laephotis capensis", "Laephotis")
sif[sif$ScientificName == "Neoromicia nana", c("ScientificName", "Genus")] = c("Afronycteris nanus", "Afronycteris")

sif[sif$ScientificName == "Pipistrellus alaschanicus", c("ScientificName", "Genus")] = c("Hypsugo alaschanicus", "Hypsugo")
sif[sif$ScientificName == "Pipistrellus cadornae", c("ScientificName", "Genus")] = c("Hypsugo cadornae", "Hypsugo")
sif[sif$ScientificName == "Pipistrellus savii", c("ScientificName", "Genus")] = c("Hypsugo savii", "Hypsugo")
sif[sif$ScientificName == "Pipistrellus eisentrauti", c("ScientificName", "Genus")] = c("Nycticeinops eisentrauti", "Nycticeinops")
sif[sif$ScientificName == "Pipistrellus deserti", c("ScientificName")] = c("Pipistrellus kuhlii")  # synonymy
sif[sif$ScientificName == "Pipistrellus rueppellii", c("ScientificName", "Genus")] = c("Vansonia rueppellii", "Vansonia")
sif[sif$ScientificName == "Pipistrellus subflavus", c("ScientificName", "Genus")] = c("Perimyotis subflavus", "Perimyotis")
sif[sif$ScientificName == "Rhogeessa aeneus", c("ScientificName")] = c("Rhogeessa aenea")  
sif[sif$ScientificName == "Rhogeessa alleni", c("ScientificName", "Genus")] = c("Baeodon alleni", "Baeodon")
sif[sif$ScientificName == "Rhogeessa gracilis", c("ScientificName", "Genus")] = c("Baeodon gracilis", "Baeodon")
sif[sif$ScientificName == "Nycticeinops schlieffeni", c("ScientificName")] = c("Nycticeinops schlieffenii")
sif[sif$ScientificName == "Harpiocephalus mordax", c("ScientificName", "Genus")] = c("Hypsugo mordax", "Hypsugo")
sif[sif$ScientificName == "Falsistrellus petersi", c("ScientificName", "Genus")] = c("Hypsugo petersi", "Hypsugo")
sif[sif$ScientificName == "Arielulus aureocollaris", c("ScientificName", "Genus")] = c("Thainycteris aureocollaris", "Thainycteris")


#add this family Cistugidae
sif[sif$ScientificName == "Cistugo seabrae", c("ScientificName", "Family")] = c("Cistugo seabrae", "Cistugidae")

# add this family
#batnames = read.csv("~/OneDrive/project/Ling/batnames/Miniopteridae2024-10-16.csv")
#batnames$ScientificName = paste0(batnames$Genus, " ", batnames$Species)
#sif[sif$ScientificName %in% batnames$ScientificName, c("Family")] = c("Miniopteridae")
sif[sif$Genus == "Miniopterus", c("Family")] = c("Miniopteridae")

# add this family
batnames = read.csv("~/OneDrive/project/Ling/batnames/Rhinonycteridae2024-10-16.csv")
batnames$ScientificName = paste0(batnames$Genus, " ", batnames$Species)
sif[sif$ScientificName %in% batnames$ScientificName, c("Family")] = c("Rhinonycteridae")
#batnames[batnames$Genus == "Paratriaenops", "ScientificName"]
sif[sif$ScientificName == "Paratriaenops furculus", c("ScientificName", "Family")] = c("Paratriaenops furcula", "Rhinonycteridae")
rm(batnames)

sif[sif$ScientificName == "Nyctimene albiventer", c("ScientificName")] = c("Nyctimene albiventris")
## subspecies 


for(f in unique(sif$Family)){
  data = sif[sif$Family ==f, ]
  #assign(paste0(f, "_data"), check_bat_names(data, f, data=T))
  check_bat_names(data, f)
  ## this function checks if "ScientificName" is consistent with batnames
}
##clean!!

## check and add subfamilies (phyllostomids were based on the diet data rojas2018; double check)
subfamilies = list.files("~/OneDrive/project/Ling/batnames/subf/")
for(f in subfamilies){
  batnames = read.csv(paste0("~/OneDrive/project/Ling/batnames/subf/", f))
  batnames$Species = gsub(" ", "", batnames$Species)
  batnames$ScientificName = paste0(batnames$Genus, " ", batnames$Species)
  
  subf = gsub("2024.*", "", f)
  
  if(subf %in% sif$Subfamily){
    if(!all(sif[sif$Subfamily == subf, "ScientificName"] %in% batnames$ScientificName)){
      print(f)
      data = sif[sif$Subfamily == subf,]
      data[!(data$ScientificName %in% batnames$ScientificName),]
      # Lophostoma evotis. evote
    }
  } else {
    ## pteropodids
    #print("not in!")
    sif[sif$ScientificName %in% batnames$ScientificName, "Subfamily"] = subf
  }
}
## clean!

write.csv(sif[, c("ScientificName", "vertlife_name", "Order", "Clade", "Family", "Subfamily", "Genus", "genes")],
          "bat624/vertlife_bat624.csv", row.names = F, quote = F)

sif[sif$ScientificName %in% sif[which(duplicated(sif$ScientificName)), "ScientificName"], 
    c("vertlife_name", "ScientificName")]
#vertlife_name      ScientificName
#3389       Myotis_aurascens      Myotis davidii --> synonymy
#3409         Myotis_davidii      Myotis davidii
#3404     Myotis_ciliolabrum  Myotis ciliolabrum
#3447    Myotis_melanorhinus  Myotis ciliolabrum --> synonymy
#3437         Myotis_laniger      Myotis laniger
#3385     Myotis_annamiticus      Myotis laniger --> synonymy
#3527      Natalus_saturatus  Natalus stramineus --> subspecies of Natalus stramineus
#3528     Natalus_stramineus  Natalus stramineus 
#4313   Pipistrellus_deserti Pipistrellus kuhlii --> synonymy
#4323    Pipistrellus_kuhlii Pipistrellus kuhlii
#4644     Pteropus_giganteus   Pteropus vampyrus --> synonymy
#4688      Pteropus_vampyrus   Pteropus vampyrus

##################################################################
############################### check monophyly (no genus) ##################################
## download the full node-dated mammal tree from https://data.vertlife.org/
tree.files = read.table("mammal.trees.all")

families = unique(sif$Family) 
subf = unique(sif$Subfamily)
subf = subf[subf != "none"]
genus = unique(sif$Genus)

## check monophyly of family, phyllostomid subfamily, and genus
mono_test = as.data.frame(matrix(nrow=length(c(families, subf, genus)), ncol=(10000 +2)))
names(mono_test) = c("level", "taxa", paste0("tree_", 1:10000))
mono_test$taxa = c(families, subf, genus)
mono_test[mono_test$taxa %in% families, "level"] = "Family"
mono_test[mono_test$taxa %in% subf, "level"] = "Subfamily"
mono_test[mono_test$taxa %in% genus, "level"] = "Genus"


for(i in 1:10000){
  tree = read.tree(paste0("~/Desktop/DNAonly_4098sp_topoFree_NDexp/", tree.files[i,1]))
  tree = keep.tip(tree, sif$vertlife_name)
  
  # test monophyly
  for(t in mono_test$taxa) {
    my.tips = sif[sif$Family == t | sif$Subfamily == t | sif$Genus == t, "vertlife_name"]
    
    if(!(is.monophyletic(tree, tips=my.tips))){
      mono_test[mono_test$taxa == t, paste0("tree_",i)] = "not_monop"
    }
  }
}

mono_test$count = apply(mono_test, 1, FUN=function(x){
  sum(!is.na(x[3:10002]))
})
mono_test[mono_test$count != 0, c(1,2,10003)]
#       level           taxa count
#14     Family    Mormoopidae  9895
#26  Subfamily  Desmodontinae     1
write.csv(sif[, c("ScientificName", "vertlife_name", "Order", "Clade", "Family", "Subfamily", "Genus", "genes")],
          "bat624/vertlife_bat624.csv", row.names = F, quote = F)



mono_subf = as.data.frame(matrix(nrow=length(subfamilies), ncol=(10000 +2)))
names(mono_subf) = c("level", "taxa", paste0("tree_", 1:10000))
mono_subf$taxa = gsub("2024.*", "", subfamilies)
mono_subf$level = "Subfamily"

for(i in 1:10000){
  tree = read.tree(paste0("~/Desktop/DNAonly_4098sp_topoFree_NDexp/", tree.files[i,1]))
  tree = keep.tip(tree, sif$vertlife_name)
  
  # test monophyly
  for(t in mono_subf$taxa) {
    my.tips = sif[sif$Subfamily == t, "vertlife_name"]
    
    if(!(is.monophyletic(tree, tips=my.tips))){
      mono_subf[mono_subf$taxa == t, paste0("tree_",i)] = "not_monop"
    }
  }
}

mono_subf$count = apply(mono_subf, 1, FUN=function(x){
  sum(!is.na(x[3:10002]))
})
mono_subf[mono_subf$count != 0, c(1,2,10003)]
#level          taxa count
#3    NA Desmodontinae     1

mono_test = read.csv("bat624/mono_test.csv")
all(mono_test[mono_test$taxa == "Desmodontinae", ] == mono_subf[mono_subf$taxa == "Desmodontinae", ], na.rm = T)
#T
mono_test = rbind(mono_test, mono_subf)
mono_test = mono_test[mono_test$level != "Genus",]
write.csv(mono_test, "bat624/mono_test.csv", row.names = F, quote = F)

## monophyly of all families except    Mormoopidae  9895
# monophyly of all fruit eating subfamilies except  Desmodontinae     1


##############################################################
############################### topology ############################
taxa = c(unique(sif$Family), unique(sif$Subfamily))
taxa = taxa[!(taxa %in% c("Mormoopidae", "Pteropodidae", "Phyllostomidae", "none", "Desmodontinae"))]
taxa = c(taxa, unique(sif[sif$Family == "Mormoopidae", "Genus"]))
taxa = c(taxa, "Desmodus")
## 37 taxa
sif$topology = "NO"
for(x in taxa){
  data = sif[sif$Family == x | sif$Subfamily == x | sif$Genus == x, ]
  sif[sif$vertlife_name == data[1, "vertlife_name"], "topology"] = x
}
nrow(sif[sif$topology != "NO",])
write.csv(sif, "bat624/vertlife_bat624.csv", row.names = F, quote = F)


## keep represents
taxa.trees = NULL
for(i in 1:10000){
  tree = read.tree(paste0("~/Desktop/DNAonly_4098sp_topoFree_NDexp/", tree.files[i,1]))
  taxa.trees = append(taxa.trees, keep.tip(tree, sif[sif$topology != "NO", "vertlife_name"])) 
}
all.topology = unique.multiPhylo(taxa.trees, use.edge.length = FALSE, use.tip.label = TRUE)
#274 unique trees
all.topology.summary = data.frame(matrix(nrow=length(all.topology), ncol=3))
names(all.topology.summary) = c("topology", "trees", "count")
all.topology.summary$topology = 1:length(all.topology)
for(j in 1:nrow(all.topology.summary)){
  all.trees = NULL
  count.trees = 0
  
  for(i in 1:10000){
    if (all.equal.phylo(taxa.trees[[i]], all.topology[[j]], use.edge.length = FALSE)) {
      all.trees = paste0(all.trees, "; tree", i)
      count.trees = count.trees + 1
    }
  }
  
  all.topology.summary[j, "trees"] = all.trees
  all.topology.summary[j, "count"] = count.trees
}
sum(all.topology.summary$count)
# 10000
all.topology.summary$percent = all.topology.summary$count / 10000
all.topology.summary = all.topology.summary[order(all.topology.summary$percent, decreasing = T),]
sum(all.topology.summary[1:5, "percent"])

sif$label = sif$vertlife_name
  
pdf("bat624/topology_taxa.pdf")
for(i in 1:length(all.topology)){
  tree = left_join(all.topology[[i]], sif[, c("label", "topology")], by="label")
  print(ggtree(tree) + 
          labs(title=paste0(all.topology.summary[i, "percent"]*100, "% trees")) +
          geom_tiplab(aes(label=topology)) +
          hexpand(0.2,1))
}
dev.off()

## some higher-level divergence is different from 21F trees because they were contraint in the patches
## but this should not impact the reconstructions of families & subfamilies (almost all monophyletic)
## genus will be done too but this seems more messy

save(taxa.trees, all.topology, all.topology.summary, 
     file="bat624/bat264_10k_topology.RData")



load("bat624/bat264_10k_topology.RData")
trees = read.table("sif179/tree100_batch1.txt")
trees = paste0("tree", trees$V1)

all.topology.summary$sampled = sapply(all.topology.summary$trees, FUN=function(x){
  if(any(trees %in% unlist(strsplit(x, "; ")))){
    return(sum(trees %in% unlist(strsplit(x, "; "))))
  } else {
    return(0)
  }
})
nrow(all.topology.summary[all.topology.summary$sampled != 0 , ])
## 39

all.topology.summary[1:4, ]
sum(as.numeric(all.topology.summary[1:4, "count"]))/10000 * 100


##############################################################
######################### add diets in sif ###############################
### exported excel sheet
diet_data = read.csv("~/Desktop/bat624_diet.csv")
all(sif$vertlife_name %in% diet_data$vertlife_name)
bat624 = merge(sif, diet_data[, 1:12], by="vertlife_name")

names(bat624)
all(bat624$ScientificName.x == bat624$ScientificName.y)
#TRUE
bat624$ScientificName = bat624$ScientificName.x

write.csv(bat624[, c(1,3:10, 12:22)], "bat624/vertlife_bat624.csv", row.names = F, quote = F)

nrow(bat624[bat624$conservative_coding != "NO", ])
#8

##############################################################
###### 100 trees batch1: node.label, prep input ##################
trees = read.table("sif179/tree100_batch1.txt")
tree.files = read.table("mammal.trees.all")

## sub-sample 100 trees: label nodes, prep inputs for MCMC
taxa = c(unique(sif$Clade), unique(sif$Family), unique(sif$Subfamily))
taxa = taxa[taxa != "none"]
taxa = taxa[taxa != "Mormoopidae"]
taxa = c(taxa, "root")
taxa = c(taxa, "SRCG", "early_burst", "EB_Mic", "PM", "super_Noctilionoidea_M")
  
  
tree_100_sub = NULL
# do not include genus
for(i in trees$V1){
  tree = read.tree(paste0("~/Desktop/DNAonly_4098sp_topoFree_NDexp/", tree.files[i,]))
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
    
    if(is.monophyletic(tree, my.tips)){
      ## name the mrca node
      tree$node.label[tree$node.label == getMRCA(tree, my.tips)] = t 
    }  else {
      print(t)
    }
  }
  
  ## save to the final tree set
  tree_100_sub = append(tree_100_sub, tree) 
}
## all the above monophyly


diet.states = c("Absent" = "0", "Complementary" = "1", "Predominant" = "2", "Strict" = "3")
diet.in = diets[!(diets %in% c("Leaves.and.flower.pieces", "Seed"))]

save(tree_100_sub, taxa, diet.in, diet.states, sif, file="bat624/bat624_100tree_batch1.RData")


unique(sif[sif$Leaves.and.flower.pieces != 0 | sif$Seed != 0, c("Pollen.and.nectar", "Fruit")])
#   Pollen.and.nectar Fruit
#13                  1     2
#51                  1     1
#144                 2     2
#190                 2     1
#248                 0     1


load("bat624/bat624_100tree_batch1.RData")

for(i in 1:nrow(sif)){
  if(any(sif[i, diet.in] ==3)){
    if(sum(sif[i, diet.in] ==0) != 5){
      print(sif[i])
    }
  }
  if(all(sif[i, diet.in] < 2)){
    print(sif[i,])
  }
}

sif[sif$ScientificName == "Chiroderma doriae", c("Fruit", "diet_source")] = c(2, "Rojas2018_fruit_update_to_2")

write.csv(sif, "bat624/vertlife_bat624.csv", row.names = F, quote = F)

save(tree_100_sub, taxa, diet.in, diet.states, sif, file="bat624/bat624_100tree_batch1.RData")


##############################################################
##################################################################