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
FUN_PlotData_FIA <- function(states = c("DE","MD"), nCores = parallel::detectCores()/2){
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
  if(!file.exists(paste0(Dir.PlotNets.FIA, "/FIA_SitesBiomes.html"))){
    FIAPlots_df <- unique(FIA_df$PLOT[, c("LON", "LAT")]) # plot extraction with lat and lon
    FIAPlots_df <- na.omit(FIAPlots_df) # remove any NA rows
    ## plot biomes and FIA sites as mapview object
    Biomes_mv <- mapview(FIAMerged_shp, color = "black", layer.name = "Biomes", zcol = "Names", col.regions = heat.colors)
    FIA_mv <- mapview(FIAPlots_df, xcol = "LON", ycol = "LAT", legend = FALSE, cex = 3.5, layer.name = "FIA Sites", color = "Black", grid = FALSE)
    save(Biomes_mv, FIA_mv, FIAPlots_df, FIAMerged_shp, file = file.path(Dir.Plots, "MapData_FIA.RData"))
    ## may need to run remotes::install_github("r-spatial/mapview") to get this to work
    mapshot(Biomes_mv+FIA_mv, url = paste0(Dir.PlotNets.FIA, "/FIA_SitesBiomes.html"), fgb = FALSE)
  }
  
  ## CALCULATION OF FITNESS AS APPROXIMATED BY BIOMASS
  if(!file.exists(file.path(Dir.Plots, "FIABiomes_df.rds"))){
    FIABiomass_df <- biomass(db = FIA_df, # which data base to use
                             polys = FIAMerged_shp, 
                             bySpecies = TRUE, # group by Species
                             byPlot = TRUE, # group by plot
                             nCores = nCores, # use half the machines cores for this
                             treeType = "live" 
    )
    saveRDS(FIABiomass_df, file.path(Dir.Plots, "FIABiomes_df.rds"))
  }else{
    load(file.path(Dir.Plots, "FIABiomes_df.rds"))
  }
  
  ## SPLITTING INTO BIOMES
  FIASplit_ls <- split(FIABiomass_df, FIABiomass_df$polyID) # extract for each biome to separate data frame
  names(FIASplit_ls) <- FIAMerged_shp@data$Names # apply correct names of biomes
  FIASplit_ls <- FIASplit_ls[-c((length(FIASplit_ls)-1):length(FIASplit_ls))] # remove 98 and 99 biome which is barren and limnic
  for(Biome_Iter in 1:length(FIASplit_ls)){
    BiomeName <- names(FIASplit_ls)[Biome_Iter]
    print(BiomeName)
    if(file.exists(file.path(Dir.Plots, paste0("FIABiome", Biome_Iter, ".RData")))){next()}
    FIAIter_df <- FIASplit_ls[[Biome_Iter]]
    FIAIter_df <- FIAIter_df[,c("pltID", "BIO_ACRE", "SCIENTIFIC_NAME", "nStems", "YEAR")] # select columns we need
    colnames(FIAIter_df) <- c("plot", "biomass", "focal", "Number at Plot", "Year") # assign new column names
    Species_vec <- unique(FIAIter_df$focal) # identify all species names
    Plots_vec <- unique(FIAIter_df$plot) # identify all plot IDs
    Interaction_ls <- as.list(rep(NA, length = length(Plots_vec))) # establish empty list with one slot for each plot
    names(Interaction_ls) <- Plots_vec # set name of the list positions to reflect the plotIDs
    counter <- 1 # create counter to index where to put the data frame created in the loop into the list
    for(Iter_plot in Plots_vec){ # plot loop: loop over all plots
      Iter_df <- FIAIter_df[FIAIter_df$plot == Iter_plot, ] # select data for currently focussed plot
      Iter_df <- Iter_df[order(Iter_df$Year, Iter_df$focal), ] # sort by year first and then by focal alphabetically in each year
      # Iter_df <- group_by(.data = Iter_df, .dots=c("focal", "Year")) %>% # group by species
      #   summarise_at(.vars = c("biomass", "Number at Plot"), .funs = median) # summarise biomass and number grouped by species
      
      Counts_mat <- matrix(rep(0, length = dim(Iter_df)[1]*length(Species_vec)), nrow = dim(Iter_df)[1]) # create empty matrix
      colnames(Counts_mat) <- Species_vec # assign column names of species names
      for(k in unique(Iter_df$Year)){
        k_df <- Iter_df[Iter_df$Year == k,] # select data for year at that plot
        Matches_vec <- base::match(x = k_df$focal, table = Species_vec) # identify position of matches of species names
        Counts_k <- rep(Iter_df$`Number at Plot`[Iter_df$Year == k], # repeat counts of each species
                        each = sum(Iter_df$Year == k) # as often as there are observations at that plot in that year
        )
        Counts_mat[Iter_df$Year == k,Matches_vec] <- Counts_k # save counts to right positions in neighbour matrix
      }
      Interaction_ls[[counter]] <- cbind(as.data.frame(Iter_df), as.data.frame(Counts_mat)) # combine matrix with biomass data
      counter <- counter +1 # raise counter
    }
    Interaction_df <- bind_rows(Interaction_ls, .id = "plot") # combine data frames in list elements into one big data frame
    Interaction_df <- cbind(Interaction_df$plot, Interaction_df$biomass, Interaction_df$focal, Interaction_df$Year, Interaction_df[,6:dim(Interaction_df)[2]])
    colnames(Interaction_df)[1:4] <- c("plot", "biomass", "focal", "year")
    # Interaction_df <- Interaction_df[-c(which(Interaction_df$biomass == 0)), ] # remove 0 biomass entries
    print(paste("Data dimensions:", paste(dim(Interaction_df), collapse = " & ")))
    save(BiomeName, FIAIter_df, Interaction_df, file = file.path(Dir.Plots, paste0("FIABiome", Biome_Iter, ".RData")))
  }
}

if(sum(file.exists(file.path(Dir.Plots, paste0("FIABiome", 1:12, ".RData")))) != 12){
  FUN_PlotData_FIA(states = c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA", "HI", "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD", "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY", "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY"), nCores = parallel::detectCores()/2)
}

####### INTERACTION NETWORK FUNCTIONALITY ---------------------------------------------------------
################# PREPARE DATA ###
# works on a data frame where the first four columns are plotID, fitness proxy (identified with Fitness argument), focal (species membership), focalID (speciesID and plotID)
Fun_StanList <- function(Fitness = "fit", data = NULL){
  ## Basic List and Data Cleaning ----
  stan.data <- list() # set up the data in list format as preferred by STAN
  data[ , -c(1:4)] <- sapply(data[ , -c(1:4)], as.numeric)
  data$focal <- as.character(data$focal)
  
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
    png(paste0(results_folder, '/tplot_', x, '.png'), width = 800, height = (N/5*70))
    tplot <- rstan::traceplot(fit, pars = x, inc_warmup = F, ncol = 5)
    print(tplot)
    dev.off()
    
    # plot posterior uncertainty intervals
    png(paste0(results_folder, '/postui_', x, '.png'), width = 600, height = (N/5*70))
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
nSamples <- 10000
nWarmup <- 500
nChains <- 5

FIABiomes_fs <- list.files(path = Dir.Plots, pattern = "FIABiome")

# probably only want to run this for 3, 6, 7, 10, 11
# Model_Iter = 11 # or 7&10&11 for high data

for(Model_Iter in c(3,6,7,11)){ # 10 for biggest data set
  load(file.path(Dir.Plots, FIABiomes_fs[[Model_Iter]]))
  print(BiomeName)
  print(dim(Interaction_df))
  
  Dir.FIABiome <- file.path(Dir.PlotNets.FIA, FIABiomes_fs[[Model_Iter]])
  dir.create(Dir.FIABiome)
  
  Test_df <- Interaction_df[Interaction_df$biomass != 0, ]
  
  tmp <- Test_df[, 1:3]
  tmp$focalID <- with(Test_df, paste(plot, focal, year, sep ="_"))
  Test_df <- cbind(tmp, Test_df[, 4:ncol(Test_df)])
  colnames(Test_df)[1:5] <- c("plotID", "fit", "focal", "focalID", "time")
  
  
  Mal_df <- Test_df[,-5] # data frame without time
  
  Mal_df <- Mal_df[, -(which(colnames(Mal_df)[-c(1:4)] %nin% Mal_df$focal)+4)]
  
  FIA_StanList <- Fun_StanList(Fitness = "fit", data = Mal_df)
  
  #### Preferences
  rstan_options(auto_write = TRUE)
  rstan_options(javascript = FALSE)
  options(mc.cores = nChains) 
  
  
  ## CHECKING DATA
  focalID <- unique(Mal_df$focal)
  neighbourID <- colnames(Mal_df[ , -c(1:4)])
  Fun_PreCheck(data = Mal_df)
  ## RUNNING MODEL
  fit <- stan(file = 'Supplement - StanModel.stan',
              data =  FIA_StanList,               # named list of data
              chains = nChains,
              warmup = nWarmup,          # number of warmup iterations per chain
              iter = nSamples,            # total number of iterations per chain
              refresh = 100,         # show progress every 'refresh' iterations
              control = list(max_treedepth = 10)
  )
  save(fit, file = file.path(Dir.FIABiome, "Model.RData"))
  
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
  try(stan_model_check(fit = fit,
                       results_folder = Dir.FIABiome,
                       params = param.vec))
  # stan_post_pred_check(post.draws = joint.post.draws, # error here on missing data!!!
  #                       results_folder = Dir.Malyon,
  #                       stan.data = StanList_ls)
  # Interactions can now be divided by the appropriate scaling (intrinsic performance, and demographic rates
  # if a population dynamics model is used) in order to return per capita interaction strengths. 
  Interaction_mean <- apply(inter_mat, c(1, 2), mean) # will return the mean estimate for every interaction (NB: this is the mean of the 80% posterior interval, so will be slightly different to the mean value returned from summary(fit), which is calculated from the full posterior distribution)  
  Interaction_mean <- Interaction_mean[, order(colnames(Interaction_mean))] # sort columns
  Interaction_mean <- Interaction_mean*-1 # need to switch sign of results
  diag(Interaction_mean) <- NA
  Interaction_hpdi <- apply(inter_mat, c(1, 2), HPDI, prob = 0.89) *-1 # need to switch sign of results
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
  # Interactions_Malyon <- Interactions_igraph[order(abs(Interactions_igraph$Inter_mean), decreasing = TRUE), ]
  # Interactions_igraph <- na.omit(Interactions_Malyon)
  
  Interactions_DAG <- Interactions_igraph[order(Interactions_igraph$Actor), ]
  Interactions_DAG <- na.omit(Interactions_DAG)
  
  Signs_DAG <- with(Interactions_DAG, data.frame(
    Intrer_mean = sign(Inter_mean),
    Intrer_max = sign(Inter_max),
    Intrer_min = sign(Inter_min)
  ))
  
  Interactions_DAG <- Interactions_DAG[abs(apply(Signs_DAG, 1, sum)) == 3, ]
  
  Interactions_DAG$Actor <- gsub(x = Interactions_DAG$Actor, pattern = " ", replacement = "_")
  Interactions_DAG$Subject <- gsub(x = Interactions_DAG$Subject, pattern = " ", replacement = "_")
  
  
  ## Graph Plot
  Dag_Paths <- paste(paste(Interactions_DAG$Actor, "->", Interactions_DAG$Subject), collapse = " ")
  dag <- dagitty(x = paste0("dag {", Dag_Paths, "}"))
  tidy_dag <- tidy_dagitty(dag, layout = "fr")
  tidy_dag$data$weight <- c(abs(Interactions_DAG$Inter_mean), rep(NA, sum(is.na(tidy_dag$data$direction))))
  tidy_dag$data$label <- c(round(Interactions_DAG$Inter_mean, 2), rep(NA, sum(is.na(tidy_dag$data$direction))))
  tidy_dag$data$colour <- c(ifelse(c(Interactions_DAG$Inter_mean) > 0, "green", "red"), rep("white", sum(is.na(tidy_dag$data$direction))))
  WeightMod <- 1/max(tidy_dag$data$weight, na.rm = TRUE)
  
  DAG_plot <- ggplot(tidy_dag, aes(x = x, y = y, xend = xend, yend = yend), layout = "linear") +
    geom_dag_point(colour = "grey", shape = 18, size = 10) +
    geom_dag_text(colour = "black", size = 2) +
    geom_dag_edges_arc(aes(edge_colour = colour, edge_width = exp(weight)*WeightMod)) +
    # geom_dag_edges_diagonal(aes(edge_colour = colour, edge_width = weight*WeightMod)) + 
    labs(title = BiomeName, subtitle = "Realised Interaction (Bimler et al. Method") + 
    theme_dag()
  ggsave(DAG_plot, units = "cm", width = 32, height = 18, filename = file.path(Dir.FIABiome, "DAG.png"))
  
  dag_data <- na.omit(tidy_dag$data)
  df <- data.frame(
    A = dag_data$name, 
    B = dag_data$to)
  df.g <- graph.data.frame(d = df, directed = TRUE)
  E(df.g)$weight <- dag_data$weight
  E(df.g)$label <- dag_data$label
  E(df.g)$colour <- dag_data$colour
  
  igraph_plot <- ggraph(df.g, layout = 'linear', circular = TRUE) + 
    # geom_edge_arc(aes(col = colour, width = weight*WeightMod, alpha = ..index..)) + 
    # scale_edge_alpha('Edge direction', guide = 'edge_direction') + 
    geom_edge_arc(aes(col = colour, width = exp(dag_data$weight)*WeightMod), show.legend = FALSE) + 
    scale_edge_color_manual(values = rev(unique(dag_data$colour))) +
    geom_node_point(color = "black", size = 2) + 
    geom_node_text(aes(label = name),  repel = TRUE)+
    coord_fixed() + 
    labs(title = BiomeName, subtitle = "Realised Interaction (Bimler et al. Method)") + 
    theme_graph()
  ggsave(igraph_plot, units = "cm", width = 32, height = 32, filename = file.path(Dir.FIABiome, "IGraph.png"))
  
  diag(Interaction_mean) <- 0
  Interaction_mean[abs(sign(Interaction_mean) + sign(Interaction_min) + sign(Interaction_max)) != 3] <- 0
  
  
  
  fv.colors = colorRampPalette(c("red","white","green")) ## define the color ramp
  colorlut = fv.colors(100)[c(1,seq(50,100,length.out=99))] ## select colors to use
  
  seq(sum(Interaction_mean>0))
  
  jpeg(filename = file.path(Dir.FIABiome, "Interactions.jpeg"), units = "cm", width = 32, height = 18, res = 1000)
  pheatmap(Interaction_mean, display_numbers = T, 
           color = c("red", "white","green"), 
           breaks = c(min(Interaction_mean, na.rm = TRUE), -0.01, 0.01, max(Interaction_mean, na.rm = TRUE)), main = "Actors (Columns) x Subject (Rows)",
           fontsize = 5)
  dev.off()
}
