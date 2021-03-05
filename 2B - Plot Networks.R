# ####################################################################### #
# PROJECT: [PhD; 2B - PLOT NETWORKS]
# CONTENTS: Generate plot data base species-association/-interaction networks
# AUTHOR: Erik Kusch
# ####################################################################### #

rm(list=ls())
####### PREAMBLE + SHAPES + REFERENCES ---------------------------------------------------------
source("0 - Preamble.R")
source("0 - ShapeFiles.R")

####### FIA DATA RETRIEVAL ---------------------------------------------------------
FUN_PlotData_FIA <- function(states = c("DE","MD"), ByYear = FALSE, nCores = parallel::detectCores()/2){
  ### EXISTENCE CHECK
  Check_vec <- states %nin% substring(list.files(Dir.Plots.FIA), 1, 2)
  if(length(substring(list.files(Dir.Plots.FIA), 1, 2)) != 0){
    if(unique(substring(list.files(Dir.Plots.FIA), 1, 2) %in% states) != TRUE){stop("Your FIA directory contains data for more than the states you asked to analyse here. Please remove these or change the state argument here to include these files already present.")}
  }
  # might need to run devtools::install_github('hunter-stanke/rFIA') to circumvent "Error in rbindlist(inTables..." as per https://github.com/hunter-stanke/rFIA/issues/7
  if(sum(Check_vec) != 0){
    FIA_df <- rFIA::getFIA(states = states[Check_vec], dir = Dir.Plots.FIA, nCores = nCores) # download FIA state data and save it to the FIA directory
  }else{
    FIA_df <- rFIA::readFIA(dir = Dir.Plots.FIA) # load all of the data in the FIA directory
  }
  
  ## BIOME SHAPE PREPARATION
  FIAMerged_shp <- aggregate(FIA_shp, by = "BIOME") # Merge shapes according to biome type
  FIAMerged_shp@data$Names <- Full_Biomes[match(FIAMerged_shp@data$BIOME, Abbr_Biomes)] # Assign corresponding full text biome names
  FIAPlots_df <- unique(FIA_df$PLOT[, c("LON", "LAT")]) # plot extraction with lat and lon
  FIAPlots_df <- na.omit(FIAPlots_df) # remove any NA rows
  ## plot biomes and FIA sites as mapview object
  Biomes_mv <- mapview(FIAMerged_shp, color = "black", layer.name = "Biomes", zcol = "Names", col.regions = heat.colors)
  FIA_mv <- mapview(FIAPlots_df, xcol = "LON", ycol = "LAT", legend = FALSE, cex = 3.5, layer.name = "FIA Sites", color = "Black", grid = FALSE)
  save(Biomes_mv, FIA_mv, FIAPlots_df, FIAMerged_shp, file = file.path(Dir.Plots, "MapData_FIA.RData"))
  ## may need to run remotes::install_github("r-spatial/mapview") to get this to work
  mapshot(Biomes_mv+FIA_mv, url = paste0(Dir.PlotNets.FIA, "/FIA_SitesBiomes.html"), fgb = FALSE)
  
  ## MASKING BY BIOMES
  FIA_Fitness_ls <- as.list(rep(NA, nrow(FIAMerged_shp@data))) # creating an empty list to hold the fitness calculations for each biome separately
  names(FIA_Fitness_ls) <- FIAMerged_shp@data$BIOME # set names of list elements with biome IDs
  FIA_Fitness_ls <- FIA_Fitness_ls[!(names(FIA_Fitness_ls) %in% c("98", "99"))] # remove ID 98 and 99 which correspond to non-vegetated areas (lakes and barren areas)
  FIA_ls <- FIA_Fitness_ls # empty list to hold the clipped FIA data bases
  
  ## CALCULATION OF FITNESS AS APPROXIMATED BY BIOMASS
  counter2 <- 1
  for(Clip_Iter in as.numeric(names(FIA_ls))){
    print(paste("Handling FIA data for", Full_Biomes[which(Abbr_Biomes == Clip_Iter)]))
    Clip_Shp <- FIAMerged_shp[FIAMerged_shp$BIOME == Clip_Iter, ]
    plot(Clip_Shp, main = Full_Biomes[which(Abbr_Biomes == Clip_Iter)])
    ClipFIA_df <- NULL
    try(ClipFIA_df <- clipFIA(db = FIA_df, mask = Clip_Shp, matchEval = TRUE), silent = TRUE)
    if(is.null(ClipFIA_df)){
      FIA_ls[[counter2]] <- "No Data"
      FIA_Fitness_ls[[counter2]] <- "No Data"
      print("No Data")
      next()
    }
    FIA_ls[[counter2]] <- ClipFIA_df
    FIAIter_df <- biomass(db = ClipFIA_df, # which data base to use
                          bySpecies = TRUE, # group by Species
                          byPlot = TRUE, # group by plot
                          nCores = nCores # use half the machines cores for this
    )
    FIAIter_df <- FIAIter_df[,c("pltID", "BIO_ACRE", "SCIENTIFIC_NAME", "nStems", "YEAR")] # select columns we need
    colnames(FIAIter_df) <- c("plot", "biomass", "focal", "Number at Plot", "Year") # assign new column names
    Species_vec <- unique(FIAIter_df$focal) # identify all species names
    if(ByYear == TRUE){
      FIAIter_df$plot <- with(FIAIter_df, paste(plot, Year, sep="_"))
    }
    FIAIter_df <- FIAIter_df[,-5]
    Plots_vec <- unique(FIAIter_df$plot) # identify all plot IDs
    Interaction_ls <- as.list(rep(NA, length = length(Plots_vec))) # establish empty list with one slot for each plot
    names(Interaction_ls) <- Plots_vec # set name of the list positions to reflect the plotIDs
    counter <- 1 # create counter to index where to put the data frame created in the loop into the list
    for(Iter_plot in Plots_vec){ # plot loop: loop over all plots
      Iter_df <- FIAIter_df[FIAIter_df$plot == Iter_plot, ] # select data for currently focussed plot
      Iter_df <- group_by(.data = Iter_df, .dots="focal") %>% # group by species
        summarise_at(.vars = c("biomass", "Number at Plot"), .funs = median) # summarise biomass and number grouped by species
      Matches_vec <- na.omit(base::match(x = Iter_df$focal, table = Species_vec)) # identify position of matches of species names
      Counts_mat <- matrix(rep(0, length = dim(Iter_df)[1]*length(Species_vec)), nrow = dim(Iter_df)[1]) # create empty matrix
      colnames(Counts_mat) <- Species_vec # assign column names of species names
      Counts_mat[, Matches_vec] <- rep(Iter_df$`Number at Plot`, each = dim(Iter_df)[1]) # save number of individuals to matrix
      Interaction_ls[[counter]] <- cbind(as.data.frame(Iter_df[ ,-3]), as.data.frame(Counts_mat)) # combine matrix with biomass data
      counter <- counter +1 # raise counter
    }
    Interaction_df <- bind_rows(Interaction_ls, .id = "plot") # combine data frames in list elements into one big data frame
    Interaction_df <- cbind(Interaction_df$plot, Interaction_df$biomass, Interaction_df$focal, Interaction_df[,4:dim(Interaction_df)[2]])
    colnames(Interaction_df)[1:3] <- c("plot", "biomass", "focal")
    
    if(ByYear == TRUE){
      focalID <- with(Interaction_df, paste(focal, plot, sep="_"))
      Interaction_df <- cbind(Interaction_df[,1:3], focalID, Interaction_df[,4:dim(Interaction_df)[2]])
    }
    Interaction_df <- Interaction_df[-c(which(Interaction_df$biomass == 0)), ] # remove 0 biomass entries
    print(paste("Data dimensions:", paste(dim(Interaction_df), collapse = " & ")))
    FIA_Fitness_ls[[counter2]] <- Interaction_df
    counter2 <- counter2 + 1
  }
  save(FIA_Fitness_ls, FIA_ls, file = file.path(Dir.Plots, "FIABiomeData_ls.RData"))
}
if(!file.exists(file.path(Dir.Plots, "FIABiomeData_ls.RData"))){
  FUN_PlotData_FIA(states = c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY", "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY"), ByYear = TRUE, nCores = parallel::detectCores()/2)
}
load(file.path(Dir.Plots, "FIABiomeData_ls.RData"))

####### INTERACTION NETWORK FUNCTIONALITY ---------------------------------------------------------
################# PREPARE DATA ###
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
  stan.data$Z <- stan.data$S*stan.data$K   # total number of possible interactions         
  
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

################# MODEL PRE-CHECKS ###
Fun_PreCheck <- function(data = NULL){
  key_speciesID <- unique(data$focal)
  key_neighbourID <- colnames(data[ , -c(1:4)])
  # message(paste0('Dataset selected: ', comm))
  message(paste0('Fitness data dimensions = ', dim(data)[1], ', ', dim(data)[2]))
  message(paste0('Number of focals = ', length(key_speciesID)))
  message(paste0('Number of neighbours = ', length(key_neighbourID)))
}

################# MODEL SPECIFICATION ###
# model specified in 'Supplement - StanModel.stan'

################# MODEL INSPECTION ###
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
    tplot <- rstan::traceplot(fit, pars = x, inc_warmup = F, ncol = 5)
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













################# RUN NETWORK METHOD FOR EACH FIA SUBSET ########################################### --------------------------
nSamples <- 800
nWarmup <- 100
nChains <- 4
#!!! remove "no data" list elements

Iter_df <- FIA_Fitness_ls[[1]]
colnames(Iter_df)[[1]] <- "PlotID"

FIA_StanList <- Fun_StanList(Fitness = "biomass", data = Iter_df)

## CHECKING DATA
focalID <- unique(Iter_df$focal)
neighbourID <- colnames(Iter_df[ , -c(1:4)]) # remove plot, biomass, focal, focalID to retain only neighbours
Fun_PreCheck(data = Iter_df)
## RUNNING MODEL #!!! error here with dimensions declared and found for irow
fit <- stan(file = 'Supplement - StanModel.stan',
            data =  FIA_StanList,               # named list of data
            chains = nChains,
            warmup = nWarmup,          # number of warmup iterations per chain
            iter = nSamples,            # total number of iterations per chain
            refresh = round(nSamples/100, 0),         # show progress every 'refresh' iterations
            control = list(max_treedepth = 10)
)
save(fit, file = file.path(Dir.Malyon, paste0(RUNNAME, "_Model.RData")))
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
                                sort(focalID),
                                neighbourID)
# inter_mat is now a 3 dimensional array, where rows = focals, columns = neighbours and 3rd dim = samples from the posterior; inter_mat[ , , 1] should return a matrix consisting of one sample for every interaction 
param.vec <- c('a', 'beta_ij', 'effect', 'response', 're', 'inter_mat', 'mu', 'sigma')
stan_model_check(fit = fit,
                 results_folder = Dir.Malyon,
                 params = param.vec)
# stan_post_pred_check(post.draws = joint.post.draws, # error here on missing data!!!
#                       results_folder = Dir.Malyon,
#                       stan.data = StanList_ls)
# Interactions can now be divided by the appropriate scaling (intrinsic performance, and demographic rates
# if a population dynamics model is used) in order to return per capita interaction strengths. 
Interaction_mean <- apply(inter_mat, c(1, 2), mean) # will return the mean estimate for every interaction (NB: this is the mean of the 80% posterior interval, so will be slightly different to the mean value returned from summary(fit), which is calculated from the full posterior distribution)  
Interaction_mean <- Interaction_mean[, order(colnames(Interaction_mean))] # sort columns
Interaction_mean <- Interaction_mean*-1 # need to switch sign of results
diag(Interaction_mean) <- NA
Interaction_hpdi <- apply(inter_mat, c(1, 2), HPDI, prob = 0.89)
Interaction_min <- Interaction_hpdi[1,,]
diag(Interaction_min) <- NA
Interaction_max <- Interaction_hpdi[2,,]
diag(Interaction_max) <- NA
Interactions_igraph <- data.frame(Actor = rep(dimnames(Interaction_mean)$neighbour, length(dimnames(Interaction_mean)$species)),
                                  Subject = rep(dimnames(Interaction_mean)$species, each = length(dimnames(Interaction_mean)$neighbour)),
                                  Inter_mean = as.vector(t(Interaction_mean)),
                                  Inter_min = as.vector(t(Interaction_min)),
                                  Inter_max = as.vector(t(Interaction_max))
)
Interactions_Malyon <- Interactions_igraph[order(abs(Interactions_igraph$Inter_mean), decreasing = TRUE), ]
Interactions_Malyon <- na.omit(Interactions_Malyon)

# Interactions_DAG <- Interactions_igraph[order(Interactions_igraph$Actor), ]
# Interactions_DAG <- na.omit(Interactions_DAG)
# ## Graph Plot
# Dag_Paths <- paste(paste(Interactions_DAG$Actor, "->", Interactions_DAG$Subject), collapse = " ")
# dag <- dagitty(x = paste0("dag {", Dag_Paths, "}"))
# tidy_dag <- tidy_dagitty(dag)
# tidy_dag$data$weight <- c(abs(Interactions_DAG$Inter_mean))
# tidy_dag$data$label <- c(Interactions_DAG$Inter_mean)
# tidy_dag$data$colour <- ifelse(c(Interactions_DAG$Inter_mean) > 0, "green", "red")
# WeightMod <- 2/max(tidy_dag$data$weight)
# # SigPos <- sign(Interactions_DAG$Inter_min) != sign(Interactions_DAG$Inter_max)
# # tidy_dag$data[SigPos, "weight"] <- 0
# # tidy_dag$data[SigPos, "colour"] <- "white"
# # tidy_dag$data[SigPos, "label"] <- NA
# 
# Graph_Malyon <- ggplot(tidy_dag, aes(x = x, y = y, xend = xend, yend = yend)) +
#   geom_dag_point(colour = "white") +
#   geom_dag_text(colour = "black", size = 10) +
#   theme_dag() + 
#   geom_dag_edges_arc(aes(edge_colour = colour, edge_width = weight*WeightMod, label = label), 
#                      angle_calc = 'along', label_dodge = grid::unit(8, "points")) + 
#   ggtitle("Realised Interactions (Malyon's Method)")

















































































































####### STAN LIST PREPARING ---------------------------------------------------------
FUN_NetInter_Prep <- function(Fitness = "biomass", data = NULL){
  
  # set up the data in list format as preferred by STAN:
  stan.data <- list()
  # create a matrix which tallies the number of observations for each focal and neighbour
  Counts <- data[ , -c(1:4)] # remove first four columns (plot, fitness proxy, focal, focalID)
  Counts[Counts>0] <- 1
  Counts <- split(Counts, data$focal)
  obs <- do.call(rbind, lapply(Counts, colSums))
  # integers
  stan.data$S <- length(unique(data$focal))  # number of species
  stan.data$N <- nrow(data)                  # number of observations
  stan.data$K <- ncol(data[ , -c(1:4)])      # number of neighbours
  stan.data$I <- length(obs[obs>0])        # number of observed interactions
  # stan.data$Z <- stan.data$S*stan.data$K   # total number of possible interactions         
  # vectors 
  stan.data$species_ID <- as.numeric(as.factor(data$focal))
  stan.data$fitness <- data[,Fitness]
  # set up indices to place observed interaction in the alpha matrix
  # first count the number of interactions observed for each focal species
  stan.data$inter_per_species <- obs
  stan.data$inter_per_species[stan.data$inter_per_species > 0] <- 1
  stan.data$inter_per_species <- rowSums(stan.data$inter_per_species)
  # column index in the alpha matrix for each observed interaction
  stan.data$icol <- unlist(apply(ifelse(obs > 0, T, F), 1, which)) 
  names(stan.data$icol) <- NULL
  stan.data$icol <- as.vector(stan.data$icol)
  # begin the row index
  stan.data$irow <- rep(1, stan.data$inter_per_species[[1]])
  # begin the start and end indices for the vector of interactions per species 
  stan.data$istart <- 1
  stan.data$iend <- stan.data$inter_per_species[[1]] #???
  # populate indices for all the other species
  for(s in 2:stan.data$S) {
    # starting position of a_ij's for i in the vector of observed interactions (ie the 1st 'j')
    stan.data$istart[s] <- sum(stan.data$inter_per_species[s-1:s])+1
    # end position of a_ij's for i in the vector of observed interactions (the last 'j')
    stan.data$iend[s] <-  sum(stan.data$inter_per_species[1:s])
    # row index in the alpha matrix for each observed interaction
    stan.data$irow <- c(stan.data$irow, rep(s, stan.data$inter_per_species[[s]]))
  }
  # model matrix 
  stan.data$X <- as.matrix(data[ , -c(1:4)])  
  # Number of observations per interaction - for model weighting?
  stan.data$Obs <- as.vector(apply(obs, 1, c))   # vector of the number of observations for each interactions
  stan.data$Obs <- stan.data$Obs[stan.data$Obs>0] # remove unobserved interactions
  return(stan.data)
}

FIA_StanList <- FUN_NetInter_Prep(Fitness = "biomass", data = Iter_df) # for some reason, this doesn't work correctly on the server... old R-version? 



# str(FIA_StanList)
# saveRDS(object = FIA_StanList, file = file.path(Dir.Plots, "STANLISTFIA.rds")) # to be executed on home PC
# FIA_StanList <- readRDS(file = file.path(Dir.Plots, "STANLISTFIA.rds")) # to be executed on server

####### EXECUTING MALYON BIMLER'S METHOD ---------------------------------------------------------

data = FIA_df


library(rstan)
rstan_options(auto_write = TRUE)
rstan_options(javascript = FALSE)
options(mc.cores = parallel::detectCores()) 

# install.packages(c("coda","mvtnorm","devtools","loo","dagitty"))
# Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
# devtools::install_github("rmcelreath/rethinking")
library(rethinking)
library(reshape2)

# I need to keep note of how focals and neighbours are indexed
key_speciesID <- unique(data$focal)
key_neighbourID <- colnames(data[ , -c(1:4)])

# message(paste0('Dataset selected: ', comm))
message(paste0('Fecundity data dimensions = ', dim(data)[1], ', ', dim(data)[2]))
message(paste0('Number of focals = ', length(key_speciesID)))
message(paste0('Number of neighbours = ', length(key_neighbourID)))

#---------------------------------------------------
# Estimate interactions with a joint IFM*REM model |
#---------------------------------------------------

fit <- stan(file = 'Supplement - StanModel.stan',
            data =  FIA_StanList,               # named list of data
            chains = 4,
            warmup = 1000,          # number of warmup iterations per chain
            iter = 5000,            # total number of iterations per chain
            refresh = 100,         # show progress every 'refresh' iterations
            control = list(max_treedepth = 10)
)
# You will probably get a few error warnings - including one that says certain chains won't have converged, 
# this is because some the parameters (unobserved interactions) are set to 0 so the Rhat value is NA 
# though there could also be other error codes, of course...

# parameters of interest
param.vec <- c('disp_dev', 'a', 'interactions','effect', 'response', 'sigma_alph', 'mu', 'ifm_alpha', 're')




fit <- readRDS(file.path(Dir.Base, "Y - Shared/Malyon", "FIA_InterNet.rds"))
library(igraph)


source(file.path(Dir.Base, "5 - Collaborators & Contributions/Interaction Network Method (Malyon Bimler)", "return_inter_array.R"))


# Get the full posteriors 
joint.post.draws <- extract.samples(fit)

# Select parameters of interest
param.vec <- c('a', 'beta_ij', 'effect', 'response', 're', 'inter_mat',
               'mu', 'disp_dev', 'sigma_alph')

# Draw 1000 samples from the 80% posterior interval for each parameter of interest
p.samples <- list()
p.samples <- sapply(param.vec[param.vec != 'sigma_alph' & param.vec != 'inter_mat'], function(p) {
  p.samples[[p]] <- apply(joint.post.draws[[p]], 2, function(x){
    sample(x[x > quantile(x, 0.1) & x < quantile(x, 0.9)], size = 1000)
  })  # this only works for parameters which are vectors
})
# there is only one sigma_alph parameter so we must sample differently:
p.samples[['sigma_alph']] <- sample(joint.post.draws$sigma_alph[
  joint.post.draws$sigma_alph > quantile(joint.post.draws$sigma_alph, 0.1) & 
    joint.post.draws$sigma_alph < quantile(joint.post.draws$sigma_alph, 0.9)], size = 1000)

# WARNING: in the STAN model for annual wildflowers, parameter 'a' lies within an exponential,
# 'a' estimates must thus be exponentiated to return estimates of intrinsic performance
intrinsic.perf <- exp(p.samples$a)
colnames(intrinsic.perf) <- focalID

# Build matrices of interaction estimates
#----------------------------------------
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
