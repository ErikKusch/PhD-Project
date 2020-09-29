# ####################################################################### #
# PROJECT: [PhD; 2B - PLOT NETWORKS]
# CONTENTS: Generate plot data base species-association/-interaction networks
# AUTHOR: Erik Kusch
# EDIT: 29/09/2020
# ####################################################################### #

####### FIA DATA RETRIEVAL ---------------------------------------------------------
PlotData_FIA <- function(states = c("DE","MD"), ByYear = FALSE, nCores = parallel::detectCores()/2){
  
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
  
  if(ByYear = TRUE){
    focalID <- with(Interaction_df, paste(focal, plot, sep="_"))
    Interaction_df <- cbind(Interaction_df[,1:3], focalID, Interaction_df[,4:dim(Interaction_df)[2]])
  }
  return(Interaction_df)
}
























