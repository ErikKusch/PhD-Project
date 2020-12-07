# ####################################################################### #
# PROJECT: [PhD - INTERACTION NETWORK TESTING] Data Simulation and Interaction Network Building
# CONTENTS: Simulate Data and Effects, Prepare Data, Specify & Build Interaction Network Models
# AUTHOR: Erik Kusch
# EDIT: 18/11/20
# ####################################################################### #

rm(list=ls())
####### PREAMBLE + SHAPES + REFERENCES ---------------------------------------------------------
source("0 - Preamble.R")
# install.packages(c("coda","mvtnorm","devtools","loo","dagitty"))
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
# devtools::install_github("rmcelreath/rethinking")
library(rethinking)
library(reshape2)
library(igraph)

# source("0 - ShapeFiles.R")
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

################# SIMULATE DATA ########################################## --------------------------
# Species never go extinct in this simulation and do not have any caps to their growth
Fun_simdata <- function(Network_df,
                        Effect = "Fitness",
                        Sim_Plots = paste("Plot", 1:50, sep = "_"), # plot names
                        Sim_Time = 1:5, # times of resampling of each plot
                        Sim_FitMean = c(10, 100), # species size classes (are assigned to the different species at random)
                        Sim_FitSpread = c(1, 20), # species size class variances (are assigned to the corresponding mean)
                        Sim_TimeEffect = .3, # change in species fitness per time step
                        Sim_PlotEffect = TRUE, # if TRUE has a 2/3 chance to either set species fitness to 50% or 150% for all species at individual plots
                        Init = 42 # to nudge seeds for random processes if I run the same specification multiple times and want to get different data sets
                       ){
  ## Basic data frame ----
  Sim_Species <- unique(c(Network_df$Actor, Network_df$Subject)) # species names
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
      Effect <- 1
      if(Effect == "Fitness"){
        Focal_ActorSpecies <- Network_df$Actor[Network_df$Subject == k] # species that act on focal
        Focal_ActorSpecies <- Focal_ActorSpecies[Focal_ActorSpecies %in% Plot_species]
        Effect <- c()
        for(z in 1:length(Focal_ActorSpecies)){
          Effect <- c(Effect, Network_df$Interaction[Network_df$Actor == Focal_ActorSpecies[z] & Network_df$Subject == k])
        }
        Effect <- mean(Effect)
      }
      
      k_fit <- specfit_df$fit[which(specfit_df$species == k)] * Effect # base fitness
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

################# MODEL INSPECTION ####################################### --------------------------
# This function takes the posterior draws for interaction estimates extracted from 
# the STAN model fit object and returns a focal x neighbour x sample array of all 
# interactions, both observed (IFM) and unrealised (RIM)
return_inter_array <- function(joint.post.draws, # posterior draws extracted using the extract.samples() function
                               response = p.samples$response, # samples for the response parameters
                               effect = p.samples$effect,     # samples for the effect parameters
                               focalID,   # vector of focal identities (names)
                               neighbourID,  # vector of neighbour identities (names)
                               ...){
  
  
  # extract the IFM interaction estimates - though observed interactions are already sampled from 
  # the 'beta_ij' parameters, sampling from inter_mat has the benefit of including columns of 0's for 
  # unrealised interactions (which can then be filled with the corresponding RIM estimates)
  ifm_inters <- joint.post.draws$inter_mat
  ifm_inters <- as.data.frame(aperm(ifm_inters, perm = c(1, 3, 2)))
  # there is now one column per interaction - unrealised interactions will simply be a column of 0's
  # sample from the 80% posterior interval for each interaction
  ifm_inters <- apply(ifm_inters, 2, function(x) {
    inter <- x[x > quantile(x, 0.1) & x < quantile(x, 0.9)]
    if (length(inter > 0)) {sample(inter, size = 1000)} 
    else {rep(0, 1000)} # this is for those unobserved interactions (0)
  })
  
  # calculate unrealised interactions 
  rim_inters <- lapply(c(1:length(focalID)), function(x) {
    # randomly re-order samples from the response parameter for each focal 
    r <- sample(response[ , x], dim(response)[1]) 
    # multiply focal's response by all neighbour effect parameters 
    return(r*effect)})
  rim_inters <- do.call(cbind, rim_inters) # this returns RI model estimates for every possible interaction 
  # every column is an interaction
  
  # replace unobserved interactions (columns of 0 in ifm_inters) with the values predicted by the RIM
  all_inters <- ifm_inters
  all_inters[ , apply(all_inters, 2, function(x) {all(x == 0)})] <- 
    rim_inters[ , apply(all_inters, 2, function(x) {all(x == 0)})]
  length(colSums(all_inters)[colSums(all_inters) == 0]) # this should be equal to 0 
  
  # reconstruct species x neighbour x sample interaction arrays  
  inter_mat <- array(data = unlist(all_inters), 
                     dim = c(nrow(all_inters), length(neighbourID), length(focalID)), 
                     dimnames = list('samples' = seq(1, nrow(all_inters)), 
                                     'neighbour' = neighbourID, 
                                     'species' = focalID))
  inter_mat <- aperm(inter_mat, c(3,2,1))
  # inter_mat is now a 3 dimensional array, where rows = focals, columns = neighbours and 3rd dim = samples from the posterior
  # inter_mat[ , , 1] should return a matrix consisting of one sample for every interaction 
  # apply(inter_mat, c(1, 2), mean) will return the mean estimate for every interaction (NB: this is the 
  # mean of the 80% posterior interval, so will be slightly different to the mean value returned from 
  # summary(fit), which is calculated from the full posterior distribution) 
  
}

# This functions spits out diagnostics for convergence
stan_diagnostic <- function(fit, 
                            results_folder,
                            ...) {
  
  # Print diagnostics
  check_hmc_diagnostics(fit)
  
  check_treedepth(fit)
  # if fit saturates the threshold, rerun with larger maximum tree depth and recheck saturation
  check_energy(fit) # check E-BFMI
  check_divergences(fit) # check divergence (validity concern)
  
  fit_summary <- summary(fit)$summary
  # get the range of rhat values
  print(c('rhat range: ', range(na.omit(fit_summary[ , 'Rhat']))))
  # and n_eff
  print(c('n_eff range: ', range(na.omit(fit_summary[ , 'n_eff']))))
  
  # print some diagnostic plots 
  png(paste0(results_folder, '/diag_rhat.png'), width=800, height=800)
  x <- stan_rhat(fit)  # distribution of Rhat
  print(x)
  dev.off()
  png(paste0(results_folder, '/diag_effSS_totalSS.png'), width=800, height=800)
  x <- stan_ess(fit)   # ratio of effective sample size to total sample size
  print(x)
  dev.off()
  png(paste0(results_folder, '/diag_MCse_postsd.png'), width=800, height=800)
  x <- stan_mcse(fit)  # ratio of Monte Carlo standard error to posterior standard deviation
  print(x)
  dev.off() 
  
  #Histograms of lp__ and accept_stat, and their joint distribution
  png(paste0(results_folder, '/diag_sample.png'), width=800, height=800)
  stan_diag(fit, 'sample')
  dev.off()
  # Violin plots showing the distributions of lp__ and accept_stat at each of the sampled step sizes (one per chain)
  png(paste0(results_folder, '/diag_stepsize.png'), width=800, height=800)
  stan_diag(fit, 'stepsize')
  dev.off()
  # Histogram of treedepth and violin plots showing the distributions of lp__ and accept_stat for each value of treedepth
  png(paste0(results_folder, '/diag_treedepth.png'), width=800, height=800)
  stan_diag(fit, 'treedepth')
  dev.off()
  # Violin plots showing the distributions of lp__ and accept_stat for iterations that encountered divergent transitions (divergent=1) and those that did not (divergent=0)
  png(paste0(results_folder, '/diag_divergence.png'), width=800, height=800)
  stan_diag(fit, 'divergence')
  dev.off()
  
  # NB: for further diagnostics, I can explore with
  # - stan_par(fit, 'specific parameter')
  # - stan_ac(fit, 'param') for auto-correlation
}



# This function checks the validity of the RE model
stan_model_check <- function(fit,
                             results_folder,
                             params = c('disp_dev', 'a', 'effect', 'response'),
                             ...){
  
  fit_summary <- summary(fit)$summary
  
  # remove 'mu' from the parameter vector and do it after because it is too long
  params <- params[params != 'mu']
  
  sapply(params, function(x) { 
    
    # this is just to length the plots so they can show all parameters
    N <- length(grep(x, rownames(fit_summary)))
    # some exceptions: 
    if (x == 'a') {N <- 20}
    if (x == 're') {N <- 800}
    if (x == 'sigma_alph') {N <- 10}
    
    # save traceplot of the posterior draws
    png(paste0(results_folder, '/tplot_', x, '.png'), width = 800, height = (N/5*100))
    tplot <- traceplot(fit, pars = x, inc_warmup = F, ncol = 5)
    print(tplot)
    dev.off()
    
    # plot posterior uncertainty intervals
    png(paste0(results_folder, '/postui_', x, '.png'), width = 600, height = (N/5*100))
    postint <- plot(fit, pars = x, show_density = T) 
    print(postint)
    dev.off()
    
  })
  
  
}


# This function plots the posterior predicted seed number vs the actual data
stan_post_pred_check <- function(post.draws,
                                 results_folder,
                                 stan.data,
                                 ...) {
  
  # currently using the loo package, can switch to rethinking
  
  # phi is the overdispersion parameter for the neg binom model
  # mu is the mean for predicted seed number 
  
  # extract mu and phi
  mu <- post.draws$mu # matrix with nrow = draws and ncol = observations
  disp_dev <- post.draws$disp_dev
  phi <- (disp_dev^2)^(-1)
  
  # generating posterior predictions
  seed_pred <- matrix(nrow = dim(mu)[1], ncol = dim(mu)[2])
  for (i in 1:dim(mu)[1]) {     # for each posterior draw
    for (j in 1:dim(mu)[2]) {    # for each observation 
      # draw from the predicted distribution
      seed_pred[i, j] <- rnbinom(1, mu = mu[i, j], size = phi[i])  
    }
  }
  # get maximum density for plot limits
  max.density <- max(c(apply(seed_pred, 1, function(x) {max(density(x)$y)}), 
                       max(density(stan.data$seeds)$y)))
  
  # dev.new(noRStudioGD = T)
  png(paste0(results_folder, '/postpredch.png'), width=800, height=800)
  # start a plot with the first draw 
  ppc.plot <- plot(density(seed_pred[1, ]), ylim = c(0, max.density), col = 'darkgrey',
                   ylab = 'Seed density',
                   main = 'Post. pred. check',
                   sub = '(grey = predicted, black = observed)') 
  for (i in 2:dim(seed_pred)[1]) {
    # add a line for each draw
    ppc.plot <- lines(density(seed_pred[i, ]), col = 'darkgrey')
  }
  # add the actual data
  ppc.plot <- lines(density(stan.data$seeds), col = 'black', lwd = 2)  
  print(ppc.plot)
  dev.off()
  
}

################# FUNCTION CALLS ######################################### --------------------------
## SIMULATING DATA
sim_df <- Fun_simdata(Network_df <- data.frame(Actor = c("A", "A", "B", "C", "D"),
                                               Subject = c("B", "D", "C", "E", "E"),
                                               Interaction = c(1.8, .4, .2, .9, 2)),
                      Effect = "Fitness",
                      Sim_Plots = paste("Plot", 1:50, sep = "_"), # plot names
                      Sim_Time = 1:5, # times of resampling of each plot
                      Sim_FitMean = c(10, 100), # species size classes (are assigned to the different species at random)
                      Sim_FitSpread = c(1, 20), # species size class variances (are assigned to the corresponding mean)
                      Sim_TimeEffect = 0.02, # change in species fitness per time step
                      Sim_PlotEffect = FALSE, # if TRUE has a 2/3 chance to either set species fitness to 50% or 150% for all species at individual plots
                      Init = 42 # to nudge seeds for random processes if I run the same specification multiple times and want to get different data sets
                      )

## REFORMATTING DATA
StanList_ls <- Fun_StanList(Fitness = "fit", # fitness column
                                  data = sim_df # data frame
                                  )
## CHECKING DATA
focalID <- unique(sim_df$focal)
neighbourID <- colnames(sim_df[ , -c(1:4)])
Fun_PreCheck(data = sim_df)

## RUNNING MODEL
fit <- stan(file = 'Supplement - StanModel.stan',
            data =  StanList_ls,               # named list of data
            chains = 4,
            warmup = 500,          # number of warmup iterations per chain
            iter = 2500,            # total number of iterations per chain
            refresh = 100,         # show progress every 'refresh' iterations
            control = list(max_treedepth = 10)
)

## MODEL DIAGNOSTICS
# Get the full posteriors 
joint.post.draws <- extract.samples(fit)
# Select parameters of interest
param.vec <- c('a', 'beta_ij', 'effect', 'response', 're', 'inter_mat', 'mu', 'sigma_alph'
               # 'disp_dev' # including this one leads to "Error in apply(joint.post.draws[[p]], 2, function(x) { : dim(X) must have a positive length"
               )
# Draw 1000 samples from the 80% posterior interval for each parameter of interest
p.samples <- list()
p.samples <- sapply(param.vec[param.vec != 'sigma_alph' & param.vec != 'inter_mat'], function(p) {
  p.samples[[p]] <- apply(joint.post.draws[[p]], 2, function(x){
    sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 1000)
  })  # this only works for parameters which are vectors
})
# there is only one sigma_alph parameter so we must sample differently:
p.samples[['sigma_alph']] <- sample(joint.post.draws$sigma[
  joint.post.draws$sigma > quantile(joint.post.draws$sigma, 0.1) & 
    joint.post.draws$sigma < quantile(joint.post.draws$sigma, 0.9)], size = 1000)
# WARNING: in the STAN model, parameter 'a' lies within a logarithmic, and must thus be logarithmitised to return estimates of intrinsic performance
intrinsic.perf <- log(p.samples$a)
colnames(intrinsic.perf) <- focalID
inter_mat <- return_inter_array(joint.post.draws, 
                                response = p.samples$response,
                                effect = p.samples$effect,
                                focalID,
                                neighbourID)
# inter_mat is now a 3 dimensional array, where rows = focals, columns = neighbours and 3rd dim = samples from the posterior
# inter_mat[ , , 1] should return a matrix consisting of one sample for every interaction 
# apply(inter_mat, c(1, 2), mean) will return the mean estimate for every interaction (NB: this is the 
# mean of the 80% posterior interval, so will be slightly different to the mean value returned from 
# summary(fit), which is calculated from the full posterior distribution) 

# Interactions can now be divided by the appropriate scaling (intrinsic performance, and demographic rates
# if a population dynamics model is used) in order to return per capita interaction strengths. 
Interaction_mean <- apply(inter_mat, c(1, 2), mean)

param.vec <- c('a', 'beta_ij', 'effect', 'response', 're', 'inter_mat', 'mu', 'sigma')
Dir.Simulation <- file.path(Dir.Plots, "Simulation")
dir.create(Dir.Simulation)
stan_model_check(fit = fit,
                 results_folder = Dir.Simulation,
                 params = param.vec)
stan_post_pred_check(post.draws = joint.post.draws, # error here on missing data!!!
                     results_folder = Dir.Simulation,
                     stan.data = StanList_ls)

# !!! need to simulate data with implemented effects of species on species






igraph_Interaction_mean <- graph_from_adjacency_matrix(adjmatrix = Interaction_mean, mode = c("directed"), 
                            weighted = TRUE, diag = TRUE, add.colnames = NULL, add.rownames = NA)

plot(igraph_Interaction_mean)

































