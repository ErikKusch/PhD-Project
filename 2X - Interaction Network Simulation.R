# ####################################################################### #
# PROJECT: [PhD - INTERACTION NETWORK TESTING] Data Simulation and Interaction Network Building
# CONTENTS: Simulate Data and Effects, Prepare Data, Specify & Build Interaction Network Models
# AUTHOR: Erik Kusch
# EDIT: 18/11/20
# ####################################################################### #

rm(list=ls())
####### PREAMBLE + SHAPES + REFERENCES ---------------------------------------------------------
source("0 - Preamble.R")
source("0 - ShapeFiles.R")
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

# install.packages(c("coda","mvtnorm","devtools","loo","dagitty"))
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
# devtools::install_github("rmcelreath/rethinking")
library(rethinking)
library(reshape2)

################# SIMULATE DATA ########################################## --------------------------
# Species never go extinct in this simulation and do not have any caps to their growth
Fun_simdata <- function(Sim_Species = c("A", "B", "C", "D", "E"), # species names
                       Sim_Plots = paste("Plot", 1:50, sep = "_"), # plot names
                       Sim_Time = 1:5, # times of resampling of each plot
                       Sim_FitMean = c(10, 100), # species size classes (are assigned to the different species at random)
                       Sim_FitSpread = c(1, 20), # species size class variances (are assigned to the corresponding mean)
                       Sim_TimeEffect = .3, # change in species fitness per time step
                       Sim_PlotEffect = TRUE, # if TRUE has a 2/3 chance to either set species fitness to 50% or 150% for all species at individual plots
                       Init = 42 # to nudge seeds for random processes if I run the same specification multiple times and want to get different data sets
                       ){
  ## Basic data frame ----
  ### basic data frame header which will be filled with simulated data
  sim_df <- data.frame(plotID = NA,
                       fit = NA,
                       focal = NA,
                       focalID = NA)
  for(i in Sim_Species){
    sim_df <- cbind(sim_df, i)
  }
  colnames(sim_df)[-c(1:4)] <- Sim_Species
  
  ## Species-Lookup Fitness ----
  ### randomly identify average fitness proxy for each species
  specfit_df <- data.frame(species = NA,
                           fit = NA)
  for(i in 1:length(Sim_Species)){ # loop over all species
    set.seed(i+Init)
    No_Fit <- sample(1:length(Sim_FitMean), 1) # random number for species
    set.seed(i+Init)
    specfit_df <- rbind(specfit_df, data.frame(species = Sim_Species[i],
                                               fit = rnorm(1, Sim_FitMean[No_Fit], Sim_FitSpread[No_Fit])
                                               )
    )
  } # end of species-loop
  specfit_df <- na.omit(specfit_df)
  
  ## Plot-Level Sampling ----
  ### random generation of plot data
  for(i in 1:length(Sim_Plots)){ # loop over plots
    set.seed(i+Init)
    No_species <- sample(2:length(Sim_Species), 1) # random number of species
    set.seed(i+Init)
    Plot_species <- sample(Sim_Species, No_species, replace = FALSE)
    ### PLOT-LEVEL EFFECTS
    if(Sim_PlotEffect == TRUE){ # plot effect check
      set.seed(i+Init)
      Nudge <- sample(c(-.5, 0, .5),1 , prob = c(1/3, 1/3, 1/3))
    }else{
      Nudge <- 0
    } # end of plot effect check
    ### SPECIES FITNESS AT PLOT
    plot_fit <- as.list(rep(NA, length(Plot_species)))
    counter <- 1
    for(k in Plot_species){ # loop over all species at the current plot
      k_fit <- specfit_df$fit[which(specfit_df$species == k)] # base fitness
      for(t in Sim_Time){ # fitness increase/decrease over time
        k_fit <- c(k_fit, k_fit[1] + k_fit[1] * Sim_TimeEffect * t)
        k_fit <- k_fit + k_fit * Nudge
      }
      plot_fit[[counter]] <- k_fit[-1] # discard the base fitness
      counter <- counter + 1
    } # end of species-loop
    ### DATA FUSING
    plot_df <- data.frame(plotID = rep(Sim_Plots[i], length(Plot_species)*length(Sim_Time)), 
                          fit = unlist(plot_fit),
                          focal = rep(Plot_species, each = length(Sim_Time)),
                          focalID = paste(rep(Sim_Plots[i], length(Plot_species)*length(Sim_Time)), 
                                          rep(Plot_species, each = length(Sim_Time)),
                                          rep(Sim_Time, length(Plot_species)),
                                          sep = "_")
    )
    ### NEIGHBOURHOOD AND SPECIES PRESENCE
    plot_neigh <- matrix(rep(0, dim(plot_df)[1]*length(Sim_Species)), ncol = length(Sim_Species))
    colnames(plot_neigh) <- Sim_Species
    present <- which(colnames(plot_neigh) %in% Plot_species)
    plot_neigh[, present] <- 1
    ### DATA EXPORT
    plot_df <- cbind(plot_df, plot_neigh)
    sim_df <- rbind(sim_df, plot_df)
  } # end of plot-loop
  sim_df <- na.omit(sim_df)
  return(sim_df) 
}

################# PREPARE DATA ########################################### --------------------------
# works on a data frame where the first four columns are plotID, fitness proxy (identified with Fitness argument), focal (species membership), focalID (speciesID and plotID)
Fun_StanList <- function(Fitness = "fit", data = NULL){
  ## Basic List and Data Cleaning ----
  stan.data <- list() # set up the data in list format as preferred by STAN
  data[ , -c(1:4)] <- sapply(data[ , -c(1:4)], as.numeric)
  
  ## Number of Observations ----
  Counts <- data[ , -c(1:4)] # remove first four columns (plot, fitness proxy, focal, focalID)
  Counts[Counts > 0] <- 1 # set all counts/abundances which indicate presence to 1
  Counts <- split(Counts, data$focal) # split into groups by species
  obs <- do.call(rbind, lapply(Counts, colSums))
  
  ## Integers ----
  stan.data$S <- length(unique(data$focal))  # number of species
  stan.data$N <- nrow(data)                  # number of observations
  stan.data$K <- ncol(data[ , -c(1:4)])      # number of neighbours
  stan.data$I <- length(obs[obs > 0])        # number of observed interactions
  # stan.data$Z <- stan.data$S*stan.data$K   # total number of possible interactions         
  
  ## Vectors ----
  stan.data$species_ID <- as.numeric(as.factor(data$focal)) # numeric indices for focal species
  stan.data$fitness <- data[, Fitness] # fitness proxy
  
  ## Indices for Interactions in Alpha-Matrix ----
  stan.data$inter_per_species <- obs # first count the number of interactions observed for each focal species
  stan.data$inter_per_species[stan.data$inter_per_species > 0] <- 1 # set all interactions to one
  stan.data$inter_per_species <- rowSums(stan.data$inter_per_species)
  stan.data$icol <- unlist(apply(ifelse(obs > 0, T, F), 1, which)) # column index in the alpha matrix for each observed interaction
  names(stan.data$icol) <- NULL
  stan.data$icol <- as.vector(stan.data$icol)
  stan.data$irow <- rep(1, stan.data$inter_per_species[[1]]) # begin the row index
  stan.data$istart <- 1 # begin the start and end indices for the vector of interactions per species 
  stan.data$iend <- stan.data$inter_per_species[[1]] # begin the start and end indices for the vector of interactions per species 
  for(s in 2:stan.data$S){ # populate indices for all the other species
    # starting position of a_ij's for i in the vector of observed interactions (ie the 1st 'j')
    stan.data$istart[s] <- sum(stan.data$inter_per_species[s-1:s])+1
    # end position of a_ij's for i in the vector of observed interactions (the last 'j')
    stan.data$iend[s] <-  sum(stan.data$inter_per_species[1:s])
    # row index in the alpha matrix for each observed interaction
    stan.data$irow <- c(stan.data$irow, rep(s, stan.data$inter_per_species[[s]]))
  }
  
  ## Model Matrix ----
  stan.data$X <- as.matrix(data[ , -c(1:4)])  
  
  ## Number of observations per interaction ----
  stan.data$Obs <- as.vector(apply(obs, 1, c)) # vector of the number of observations for each interactions
  stan.data$Obs <- stan.data$Obs[stan.data$Obs > 0] # remove unobserved interactions

  return(stan.data)
}
################# MODEL PRE-CHECKS ######################################## --------------------------
Fun_PreCheck <- function(data = NULL){
  key_speciesID <- unique(data$focal)
  key_neighbourID <- colnames(data[ , -c(1:4)])
  # message(paste0('Dataset selected: ', comm))
  message(paste0('Fitness data dimensions = ', dim(data)[1], ', ', dim(data)[2]))
  message(paste0('Number of focals = ', length(key_speciesID)))
  message(paste0('Number of neighbours = ', length(key_neighbourID)))
}
################# MODEL SPECIFICATION #################################### --------------------------
# model specified in 'Supplement - StanModel.stan'

################# FUNCTION CALLS ######################################### --------------------------
sim_df <- Fun_simdata(Sim_Species = c("A", "B", "C", "D", "E"), # species names
                      Sim_Plots = paste("Plot", 1:50, sep = "_"), # plot names
                      Sim_Time = 1:5, # times of resampling of each plot
                      Sim_FitMean = c(10, 100), # species size classes (are assigned to the different species at random)
                      Sim_FitSpread = c(1, 20), # species size class variances (are assigned to the corresponding mean)
                      Sim_TimeEffect = 0.02, # change in species fitness per time step
                      Sim_PlotEffect = FALSE, # if TRUE has a 2/3 chance to either set species fitness to 50% or 150% for all species at individual plots
                      Init = 42 # to nudge seeds for random processes if I run the same specification multiple times and want to get different data sets
                      )

StanList_ls <- Fun_StanList(Fitness = "fit", # fitness column
                                  data = sim_df # data frame
                                  )

Fun_PreCheck(data = sim_df)

fit <- stan(file = 'Supplement - StanModel.stan',
            data =  FIA_StanList,               # named list of data
            chains = 4,
            warmup = 1000,          # number of warmup iterations per chain
            iter = 5000,            # total number of iterations per chain
            refresh = 100,         # show progress every 'refresh' iterations
            control = list(max_treedepth = 10)
)
