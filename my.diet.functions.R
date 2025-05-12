library(phytools)
library(ggtree)
library(treeio)
library(ggplot2)
library(ggpubr)
library(khroma)
library(tidytree)
library(ggimage)
library(ape)
library(stringr)
library(reshape2)
library(coda)


diets = c("Arthropods", "Blood", "Terrestrial.vertebrates", 
          "Fish", "Leaves.and.flower.pieces", "Pollen.and.nectar", "Fruit", "Seed")
diet.states = c("Absent" = "0", "Complementary" = "1", "Predominant" = "2", "Strict" = "3")
state.color = c("#d2dae1", "#6baed6", "#08519c", "black")
state.color2 = c("grey60", "#6baed6", "#08519c", "white")


## functions to plot MCMC results 

# the input raw counts have been aggregated across trees (diet_summary_aggregate)
# returns the pie plot of a certain node (n) and a certain diet (d)
get_node_pie_plot = function(n, d, prob=F, trees = NULL, title = NULL, lab.size=7, ...){
  summary = get(paste0(d, "_summary_aggregate"))
  summary = summary[summary$Taxon == n, ]
  
  if(!is.null(trees)){
    summary = summary[summary$tree %in% as.numeric(trees), ]
  }
  
  data = as.data.frame(names(diet.states)[diet.states %in% unique(sif[,d])])
  names(data) = "state"
  data$estimate = sapply(data$state, FUN = function(x){
    sum(as.numeric(summary[, x]))
  })
  
  data$prob = data$estimate/sum(data$estimate)
  data = data[order(data$estimate, decreasing = T),]
  
  if(is.null(title)){
    my.title = paste0(d, "_", data[1,"state"], "\n", round(data[1,"prob"],3))
  } else {
    if(title == "no_food"){
      my.title = paste0(data[1,"state"], "\n", round(data[1,"prob"],3))
    }
  }
  
  if(prob){
    return(data)
  } else {
    return(ggplot(data, aes(x="", y=estimate, fill=state)) +
             scale_fill_manual(values=setNames(state.color, names(diet.states)), guide="none") +
             geom_bar(stat="identity", width=1) +
             coord_polar("y", start=0) + 
             labs(title=my.title) +
             theme_void() + theme(title = element_text(size=lab.size)))
  }
}


# input is diet_summary_aggregate
# returns the posterior state (counts or probabilities) of a tip species (sp; could also be a node) and a diet (d)
# print a note if the most likely state has the probability < 0.5
get_tip_posterior_HPD = function(sp, d, prob = F, ...){
  summary = get(paste0(d, "_summary_aggregate"))
  summary = summary[summary$Taxon == sp, ]
  
  data = as.data.frame(names(diet.states)[diet.states %in% unique(sif[,d])])
  names(data) = "state"
  data$estimate = sapply(data$state, FUN = function(x){
    sum(as.numeric(summary[, x]))
  })
  #sum(data$estimate)
  # 36000000 = 100 trees * 10 chains * (2000000*0.9)/50 iterations
  data$prob = data$estimate/sum(data$estimate)
  data = data[order(data$estimate, decreasing = T),]
  #print(data)
  
  ## print a note if this state has the probability < 0.5
  if(data[1,"prob"] < 0.5){
    print(paste0(sp, ", ", d, ": ", 
                 data[1,"state"], " ", round(data[1,"prob"], 3), "; ",
                 data[2,"state"], " ", round(data[2,"prob"], 3)))
  }
  
  ## the above equals to the mean value
  #data$mean = sapply(data$state, FUN = function(x){ mean(summary[, x])/360000 })

  ## estimate median and 95% HPD: aggregate across trees, each tree representing counts from 360000 iterations 
  data$median = sapply(data$state, FUN = function(x){
    median(summary[, x])/360000
  })
  data$HPD95_lower = sapply(data$state, FUN = function(x){
    HPDinterval(mcmc(summary[, x]))[1]/360000
  })
  data$HPD95_upper = sapply(data$state, FUN = function(x){
    HPDinterval(mcmc(summary[, x]))[2]/360000
  })
  
  if(prob == T){
    ## return state probabilities
    return(data)
  }
  
  ## return the state having the highest likelihood as the posterior state
  return(diet.states[data[1,"state"]])
}


# this function plots the tips that have different states of a diet between input and posterior
# update tip names to ScientificName 
plot_tip_bar = function(d, tip.posterior.data, tip.lab = "vertlife_name", ...) {
  data = NULL
  for(s in 0:3){
    if(s %in% unique(sif[,d])){
      species = sif[sif[,d] == s, tip.lab]
      
      data.plot = tip.posterior.data[tip.posterior.data$species %in% species,]
      if(nrow(data.plot) > 0){
        data.order = data.plot[data.plot$state == names(diet.states)[diet.states == s],]
        my.order = unique(data.order[order(data.order$prob, decreasing = T), "species"])
        data.plot$species = factor(data.plot$species, levels = my.order)
        data.plot = data.plot[order(data.plot$species),]
        data.plot$input = s
        data = rbind(data, data.plot)
      }
    }
  }
  data$input = as.character(data$input)
  return(ggplot(data, aes(y=species, x=prob, fill=state, label=input)) + 
           geom_bar(stat="identity") +
           scale_fill_manual(values=setNames(state.color, names(diet.states)), name=d, 
                             labels = paste0(diet.states, "_", names(diet.states))) +
           geom_text(aes(x = -0.05, colour = input)) + 
           scale_color_manual(values=c("0"="grey60", "1"="#6baed6", "2"="#08519c", "3"="black"), guide="none")+ theme_classic() + theme(panel.grid = element_line(color="white")) 
  )
}

