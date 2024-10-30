#!/usr/bin/env Rscript

## NB! run this in the dir where the results are saved

# This script checks if a model that was fitted with MCMCglmm was 
# executed for a sufficient number of MCMC generations. It does so
# by examining the effective sample size and potential scale reduction 
# factor of each model parameter.

# The script needs to be run from the command line, with the user providing the directory where the model outputs (*.Rda files) are found, and the number of traits included, e.g.:

# Rscript check_ESS_PSRF.R ../Results/MCMCglmm_fits/ 6

library(MCMCglmm, lib="~/R_4.4")

#############################
# F  U  N  C  T  I  O  N  S #
#############################

# This function checks if the effective sample size (ESS) of each model 
# parameter per chain is at least 400. If not, this means that the 
# parameter space was not adequately explored.
check_ESS <- function(model_fits, n.traits)
{
  
  # Initialise a vector to store the ESS for each model parameter per 
  # chain.
  all_ESS_vals <- c()
  
  # Store ESS values in the vector.
  for ( i in 1:length(model_fits) )
  {
    all_ESS_vals <- c(
      all_ESS_vals,
      effectiveSize(model_fits[[i]]$Sol[,1:n.traits]),
      effectiveSize(model_fits[[i]]$VCV),
      effectiveSize(model_fits[[i]]$CP)
    )
  }
  
  # Remove non-positive ESS values (which can occur if a parameter 
  # is constant).
  all_ESS_vals <- all_ESS_vals[all_ESS_vals > 0]
  
  # Then, check if any parameter has an ESS value below 400.
  low_ESS <- which(all_ESS_vals < 400)
  
  if ( length(low_ESS) > 0 )
  {
    stop('PROBLEM! The ESS is not big enough!')
  } else
  {
    cat('Everything OK!\n', sep = '')
  }
}

# This function checks if the potential scale reduction factor (PSRF; 
# Gelman & Rubin, Stat. Sci., 1992) value of each model parameter is 
# below 1.1. If not, this means that the chains did not converge on 
# statistically indistinguishable posterior distributions.
check_PSRF <- function(model_fits, n.traits)
{
  
  # Initialise lists to store the model parameter values per chain.
  Sols <- 'mcmc.list('
  VCVs <- 'mcmc.list('
  CPs  <- 'mcmc.list('
  
  # Store the values in the lists.
  for ( i in 1:length(model_fits) ) {
    Sols <- paste(
      Sols, 'model_fits[[', i, ']]$Sol[,1:', n.traits, '], ', sep = ''
    )
    VCVs <- paste(
      VCVs, 'model_fits[[', i, ']]$VCV, ', sep = ''
    )
    CPs <- paste(
      CPs, 'model_fits[[', i, ']]$CP, ', sep = ''
    )
  }
  
  Sols <- sub(', $', ')', Sols)
  VCVs <- sub(', $', ')', VCVs)
  CPs  <- sub(', $', ')', CPs)
  
  # Calculate the PSRF for each parameter.
  psrf_vals <- c(
    gelman.diag(
      eval(parse(text = Sols)), multivariate = FALSE
    )$psrf[,1],
    gelman.diag(
      eval(parse(text = VCVs)), multivariate = FALSE
    )$psrf[,1],
    gelman.diag(
      eval(parse(text = CPs)), multivariate = FALSE
    )$psrf[,1]
  )
  
  # Remove NA values (which can occur if a parameter is constant).
  psrf_vals <- psrf_vals[!is.na(psrf_vals)]
  
  # Then, check if any parameter has a PSRF value of 1.1 or higher.
  high_psrf <- which(psrf_vals >= 1.1)
  
  if ( length(high_psrf) > 0 )
  {
    print("parameters having PSRF >= 1.1:")
    print(psrf_vals[high_psrf])
    stop("PROBLEM! The chains have not converged!")
  } else
  {
    cat('Everything OK!\n', sep = '')
  }
}

plot_trace_density = function(chain_id, result.dir, ...){
  pdf(paste0(result.dir, "Trace_Density_chain", chain_id, ".pdf"), height = 3, width = 6)
  
  for(i in 1:ncol(fit$CP)){
    print(plot(fit$CP[,i], main="") + title(main=paste0(gsub("trait", "", colnames(fit$CP)[i]), "\n", "Trace", "_Density")))
  }
  for(j in diets){
    print(plot(fit$VCV[,grep(paste0(paste0("trait", j), ":", paste0("trait", j, ".vertlife_name")), 
                       colnames(fit$VCV))]) + title(main=j))
  }
  dev.off()
}


############################
# M  A  I  N    C  O  D  E #
############################

# Read the working directory provided by the user as a command line argument.
args <- commandArgs(TRUE)
n.traits = as.numeric(args[1])
plot.all.chains = as.character(args[2])

working_dir = "./"
# Look for model fit files in the working directory.
Rda_files_in_working_dir <- list.files(
  path = working_dir, pattern = '^\\d+\\.Rda$'
)

current_model_fits <- list()
for ( i in 1:length(Rda_files_in_working_dir) ) {
  
  # ... load the corresponding .Rda files...
  cat('Now loading ', Rda_files_in_working_dir[i], '...\n', sep ='')
  load(paste(working_dir, '/', Rda_files_in_working_dir[i], sep = ''))
  current_model_fits[[i]] <- fit

  diets = gsub("trait", "", colnames(fit$Sol)[1:n.traits])

  if(plot.all.chains == "T"){
    plot_trace_density(i, working_dir)
  }
}

# ... and check the ESS and PSRF.
cat('Now checking the ESS... ', sep = '')
check_ESS(current_model_fits, n.traits)

cat('Now checking the PSRF... ', sep = '')
check_PSRF(current_model_fits, n.traits)
