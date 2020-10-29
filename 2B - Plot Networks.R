# ####################################################################### #
# PROJECT: [PhD; 2B - PLOT NETWORKS]
# CONTENTS: Generate plot data base species-association/-interaction networks
# AUTHOR: Erik Kusch
# EDIT: 29/09/2020
# ####################################################################### #

####### FIA DATA RETRIEVAL ---------------------------------------------------------
FUN_PlotData_FIA <- function(states = c("DE","MD"), ByYear = FALSE, nCores = parallel::detectCores()/2){
  
  ### plot numbers with geo-coordinates?!
  
  ### EXISTENCE Check!
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
  
  FIA_df <- biomass(db = FIA_df, # which data base to use
                   bySpecies = TRUE, # group by Species
                   byPlot = TRUE, # group by plot
                   nCores = nCores # use half the machines cores for this
  )
  
  FIA_df <- FIA_df[,c("pltID", "BIO_ACRE", "SCIENTIFIC_NAME", "nStems", "YEAR")] # select columns we need
  colnames(FIA_df) <- c("plot", "biomass", "focal", "Number at Plot", "Year") # assign new column names
  
  Species_vec <- unique(FIA_df$focal) # identify all species names
  
  if(ByYear == TRUE){
    FIA_df$plot <- with(FIA_df, paste(plot, Year, sep="_"))
  }
  
  FIA_df <- FIA_df[,-5]
  
  Plots_vec <- unique(FIA_df$plot) # identify all plot IDs
  
  Interaction_ls <- as.list(rep(NA, length = length(Plots_vec))) # establish empty list with one slot for each plot
  names(Interaction_ls) <- Plots_vec # set name of the list positions to reflect the plotIDs
  counter <- 1 # create counter to index where to put the data frame created in the loop into the list
  for(Iter_plot in Plots_vec){ # plot loop: loop over all plots
    Iter_df <- FIA_df[FIA_df$plot == Iter_plot, ] # select data for currently focussed plot
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
  return(Interaction_df)
}

FIA_df <- FUN_PlotData_FIA(states = c("DE","MD"), ByYear = TRUE, nCores = parallel::detectCores()/2)
FIA_df <- FIA_df[-c(which(FIA_df$biomass == 0)),]# remove 0 biomass entries
plot(log(FIA_df$biomass) ~ rowSums(FIA_df[,-1:-4]))
########################################################################################################

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

FIA_StanList <- FUN_NetInter_Prep(Fitness = "biomass", data = FIA_df) # for some reason, this doesn't work correctly on the server... old R-version? 
# str(FIA_StanList)
# saveRDS(object = FIA_StanList, file = file.path(Dir.Plots, "STANLISTFIA.rds")) # to be executed on home PC
# FIA_StanList <- readRDS(file = file.path(Dir.Plots, "STANLISTFIA.rds")) # to be executed on server

########################################################################################################

data = FIA_df


library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) 

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





