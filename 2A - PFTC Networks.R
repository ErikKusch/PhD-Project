#' ####################################################################### #
#' PROJECT: [PhD; 2A - PFTC NETWORKS] 
#' CONTENTS: 
#'  - Generate PFTC species-association/-interaction networks
#'  DEPENDENCIES:
#'  - 0 - Preamble.R
#'  - X - Functions_Bayes.R
#'  - X - Functions_Plotting.R
#'  - PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv (available from PFTC group efforts; https://osf.io/hjpwt/)
#'  - PU.10_PFTC3.10_2020_Peru_Coordinates.xlsx (available from PFTC group efforts; https://osf.io/uk85w/)
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

# PREAMBLE ================================================================
rm(list=ls())

## Sourcing ---------------------------------------------------------------
source("0 - Preamble.R")
# source("0 - ShapeFiles.R")
source("X - Functions_Plotting.R")
source("X - Functions_Bayes.R")
source("X - Environment.R")
## Bayes Settings ---------------------------------------------------------
nSamples <- 7000
nWarmup <- 700
nChains <- 4

# DATA ====================================================================
## PFTC META & CLIMATEDATA ------------------------------------------------
## check if bioclimatic data has already been established on hard drive
if(!file.exists(file.path(Dir.PFTC, "Metadata_df.csv"))){ # bioclimatic data not established yet
  if(!file.exists(file.path(Dir.PFTC, "PU.10_PFTC3.10_2020_Peru_Coordinates.xlsx"))){
    download.file(url = "https://osf.io/uk85w/download", 
                  destfile = file.path(Dir.PFTC, "PU.10_PFTC3.10_2020_Peru_Coordinates.xlsx"), mode = "wb")
  }
  Metadata_df <- as.data.frame(readxl::read_xlsx(file.path(Dir.PFTC, "PU.10_PFTC3.10_2020_Peru_Coordinates.xlsx"))) # read metadata file, available here https://osf.io/uk85w/
  Metadata_df$SiteID <- with(Metadata_df, paste(Site, Treatment, PlotID, sep = "_")) # create ID column which combines site, treatment, and plotid
  Metadata_df$Year <- 2019 # set observation years to 2018 for all sites (that was the field season most of the data was collected throughout)
  BioClim_df <- as.data.frame(matrix(rep(NA, nrow(Metadata_df)*19), ncol = 19)) # establish empty data frame for bioclimatic information
  colnames(BioClim_df) <- paste0("BIO", 1:19) # name columns in bioclimatic data frame
  Metadata_df <- cbind(Metadata_df, BioClim_df) # combine data frames to handle just one big object
  ## download and compute bioclimatic data with user-defined function
  Metadata_df <- FUN.ClimData(Data = Metadata_df,
               ID = "SiteID",
               Lat = "Latitude",
               Lon = "Longitude",
               Year = "Year",
               Cores = numberOfCores,
               Dir = Dir.PFTC,
               Shape = NULL,
               FileName = "Metadata_df",
               force = FALSE,
               rawdata = FALSE)
  
}else{ # bioclimatic data has been established
  Metadata_df <- read.csv(file.path(Dir.PFTC, "Metadata_df.csv")) # load the metadata with attached bioclimatic data
}

## PFTC RAW DATA ----------------------------------------------------------
## download PFTC data if needed
if(!file.exists(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"))){
  download.file(url = "https://osf.io/hjpwt/download", 
                destfile = file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv"))
}
Raw_df <- read.csv(file.path(Dir.PFTC, "PFTC3-Puna-PFTC5_Peru_2018-2020_LeafTraits_clean.csv")) # load pftc data
Raw_df <- Raw_df[Raw_df$trait == "dry_mass_g", c("site", "treatment", "plot_id", "taxon", "trait", "value")] # only dry biomass rows and select only relevant columns for speedier data handling 
Raw_df$SiteID <- with(Raw_df, paste(site, treatment, plot_id, sep = "_")) # create SiteID index
Weights_df <- Raw_df[,c("SiteID", "taxon", "value")] # limit data to necessary parts

## MODEL DATA FRAMES -------------------------------------------------------
### INDEX ---------------------------------------------
Index_df <- Weights_df
Index_df <- na.omit(Index_df) # remove NA rows
colnames(Index_df) <- c("SiteID", "Species", "Outcome")
##!!! should probably remove species not recognized by taxonomy here

### NEIGHBOURS ----------------------------------------
## aggregating predictors of individual plants at plots to species-level means
Neigh_df <- aggregate.data.frame(x = Index_df$Outcome, 
                                 by = list(Index_df$SiteID, Index_df$Species), 
                                 FUN = "mean")
colnames(Neigh_df) <- c("SiteID", "Species", "Predictor") # assign more readable column names
## make data frame wide for clearer representation of data
Neigh_df <- reshape(data = Neigh_df, direction = "wide", 
                    idvar = "SiteID", timevar = "Species")
colnames(Neigh_df) <- gsub(colnames(Neigh_df), pattern = "Predictor.", replacement = "") # fix column names which reshape() messed up

### ENVIRONMENT ---------------------------------------
## isolating only climate information
Envir_df <- Metadata_df[,c("SiteID",paste0("BIO", 1:19))]
##!!! should probably remove collinear variables here
#LinCombs <- caret::findLinearCombos(na.omit(Envir_df[, -1]))[["linearCombos"]][[1]] # find which columns to remove according to https://stackoverflow.com/questions/34929208/how-can-one-list-pairs-of-perfectly-collinear-numeric-vectors-in-a-data-frame
#Envir_df <- Envir_df[, -LinCombs-1]

# ANALYSIS ================================================================
## TREATMENT IDENIFICATION ------------------------------------------------
## split SiteID so we can identify observations by Site and Treatment separately 
Treatments_df <- data.frame(Treatment = sapply(str_split(Index_df$SiteID, "_"), "[[", 2),
                            Plot = sapply(str_split(Index_df$SiteID, "_"), "[[", 1)
)
Treatments_vec <- c("ALL", as.character(unique(Treatments_df$Treatment)), as.character(unique(Treatments_df$Plot))) # find treatmend and site identifiers
Treatments_vec <- Treatments_vec[Treatments_vec != "NA"] # remove NA identifier

## TREATMENT LOOP ---------------------------------------------------------
## loop over all treatments/sites and run network methodology
for(Treatment_Iter in Treatments_vec){ 
  ### Directory Creation ----------------------------------------------
  Dir.PlotNets.PFTC.Iter <- file.path(Dir.PlotNets.PFTC, Treatment_Iter)
  if(dir.exists(Dir.PlotNets.PFTC.Iter)){
    print(paste("PFTC models already run for Treatment/Plot, ", Treatment_Iter,". You can find them here:", Dir.PlotNets.PFTC.Iter))
    next()
  }
  dir.create(Dir.PlotNets.PFTC.Iter)
  
  ### Data Preparation and Subsetting ---------------------------------
  ## Create iteration versions of the base data
  Run_Index_df <- Index_df
  Run_Neigh_df <- Neigh_df
  Run_Envir_df <- Envir_df
  ## subset if needed according to iteration
  if(Treatment_Iter != "ALL"){
    ## Identify Treatments in Envir
    EnvirTreatment_df <- data.frame(Treatment = sapply(str_split(Envir_df$SiteID, "_"), "[[", 2),
                                    Plot = sapply(str_split(Envir_df$SiteID, "_"), "[[", 1)
    )
    ## subsetting
    Run_Envir_df <- Run_Envir_df[which(EnvirTreatment_df[,1] == Treatment_Iter | EnvirTreatment_df[,2] == Treatment_Iter), ]
    Run_Envir_df <- Run_Envir_df[which(rowSums(Run_Envir_df[ ,-1], na.rm=TRUE) != 0), ] # remove all sites for which we do not have environmental data
    
    ## Identify Treatments in Neigh
    NeighTreatment_df <- data.frame(Treatment = sapply(str_split(Neigh_df$SiteID, "_"), "[[", 2),
                                   Plot = sapply(str_split(Neigh_df$SiteID, "_"), "[[", 1)
    )
    ## subsetting
    Run_Neigh_df <- Run_Neigh_df[which(NeighTreatment_df[,1] == Treatment_Iter | NeighTreatment_df[,2] == Treatment_Iter), ]
    Run_Neigh_df <- Run_Neigh_df[Run_Neigh_df$SiteID %in% Run_Envir_df$SiteID, ] # keep only sites for which we have climate data
    ## subsetting
    Run_Index_df <- Run_Index_df[which(Treatments_df[,1] == Treatment_Iter | Treatments_df[,2] == Treatment_Iter), ]
    Run_Index_df <- Run_Index_df[Run_Index_df$SiteID %in% Run_Envir_df$SiteID, ] # keep only sites for which we have climate data
  }
  
  ### Data Sanity Checks ----------------------------------------------
  Run_Neigh_df <- Run_Neigh_df[, c(1,which(colSums(Run_Neigh_df[,-1], na.rm = TRUE) != 0)+1)] # retain only neighbours which are observed at least once at each site/treatment
  
  
  # OmitCols <- -(which(colnames(Neigh_df)[-c(1:4)] %nin% Mal_df$focal)+4)
  # if(length(OmitCols) != 0){
  #   Run_df <- Run_df[, OmitCols]
  # }else{
  #   Run_df <- Run_df
  # }
  # Run_df <- Run_df[Run_df$fit != 0, ]
  
  ## MODEL DATA LIST --------------------------------------------------
  Stan_list <- FUN.StanList(
    Outcome = "Outcome", # name of the outcome variable column in Index_df
    Index_df = Run_Index_df, # data frame containing columns for ID, Fitness, Site, and Species
    Neigh_df = Run_Neigh_df, # data frame containing column for site and columns for all species; cells containing predictor variable of each species at each site
    Envir_df = Run_Envir_df # data frame containing column for site and columns for climate variables; cells containing climate variable values at each site
  )
  
  # FIA_StanList <- Fun_StanList(Fitness = "fit", data = Run_df)
  
  #### Preferences
  rstan_options(auto_write = TRUE)
  rstan_options(javascript = FALSE)
  options(mc.cores = nChains) 
  
  ## CHECKING DATA
  # focalID <- unique(Run_df$focal)
  # neighbourID <- colnames(Run_df[ , -c(1:4)])
  FUN.DataDims(data = Stan_list)
  ## RUNNING MODEL
  fit <- stan(file = 'StanModel_Development.stan',
              data =  Stan_list,               # named list of data
              chains = nChains,
              warmup = nWarmup,          # number of warmup iterations per chain
              iter = nSamples,            # total number of iterations per chain
              refresh = 100,         # show progress every 'refresh' iterations
              control = list(max_treedepth = 10)
  )
  save(fit, file = file.path(Dir.PlotNets.PFTC.Iter, "Model.RData"))
  # load(file.path(Dir.PlotNets.PFTC.Iter, "Model.RData")) # for testing
  
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
                       results_folder = Dir.PlotNets.PFTC.Iter,
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
  # Interaction_hpdi <- apply(inter_mat, c(1, 2), HPDI, prob = 0.89) *-1 # need to switch sign of results
  # Interaction_min <- Interaction_hpdi[1,,]
  # diag(Interaction_min) <- NA
  # Interaction_max <- Interaction_hpdi[2,,]
  # diag(Interaction_max) <- NA
  # Interactions_igraph <- data.frame(Actor = rep(dimnames(Interaction_mean)$neighbour, length(dimnames(Interaction_mean)$species)),
  #                                   Subject = rep(dimnames(Interaction_mean)$species, each = length(dimnames(Interaction_mean)$neighbour)),
  #                                   Inter_mean = as.vector(t(Interaction_mean)),
  #                                   Inter_min = as.vector(t(Interaction_min)),
  #                                   Inter_max = as.vector(t(Interaction_max))
  # )
  # # Interactions_Malyon <- Interactions_igraph[order(abs(Interactions_igraph$Inter_mean), decreasing = TRUE), ]
  # # Interactions_igraph <- na.omit(Interactions_Malyon)
  # 
  # Interactions_DAG <- Interactions_igraph[order(Interactions_igraph$Actor), ]
  # Interactions_DAG <- na.omit(Interactions_DAG)
  # 
  # Signs_DAG <- with(Interactions_DAG, data.frame(
  #   Intrer_mean = sign(Inter_mean),
  #   Intrer_max = sign(Inter_max),
  #   Intrer_min = sign(Inter_min)
  # ))
  # 
  # Interactions_DAG <- Interactions_DAG[abs(apply(Signs_DAG, 1, sum)) == 3, ]
  # if(nrow(Interactions_DAG) == 0){
  #   print(paste("No significantly positive or negative interactions for treatment", Treatment))
  #   next()
  # }
  # 
  # Interactions_DAG$Actor <- gsub(x = Interactions_DAG$Actor, pattern = " ", replacement = "_")
  # Interactions_DAG$Subject <- gsub(x = Interactions_DAG$Subject, pattern = " ", replacement = "_")
  # 
  # 
  # ## Graph Plot
  # Dag_Paths <- paste(paste(Interactions_DAG$Actor, "->", Interactions_DAG$Subject), collapse = " ")
  # dag <- dagitty(x = paste0("dag {", Dag_Paths, "}"))
  # tidy_dag <- tidy_dagitty(dag, layout = "fr")
  # tidy_dag$data$weight <- c(abs(Interactions_DAG$Inter_mean), rep(NA, sum(is.na(tidy_dag$data$direction))))
  # tidy_dag$data$label <- c(round(Interactions_DAG$Inter_mean, 2), rep(NA, sum(is.na(tidy_dag$data$direction))))
  # tidy_dag$data$colour <- c(ifelse(c(Interactions_DAG$Inter_mean) > 0, "green", "red"), rep("white", sum(is.na(tidy_dag$data$direction))))
  # WeightMod <- 1/max(tidy_dag$data$weight, na.rm = TRUE)
  # 
  # DAG_plot <- ggplot(tidy_dag, aes(x = x, y = y, xend = xend, yend = yend), layout = "linear") +
  #   geom_dag_point(colour = "grey", shape = 18, size = 10) +
  #   geom_dag_text(colour = "black", size = 2) +
  #   geom_dag_edges_arc(aes(edge_colour = colour, edge_width = exp(weight)*WeightMod)) +
  #   # geom_dag_edges_diagonal(aes(edge_colour = colour, edge_width = weight*WeightMod)) + 
  #   labs(title = Treatment, subtitle = "Realised Interaction (Bimler et al. Method") + 
  #   theme_dag()
  # ggsave(DAG_plot, units = "cm", width = 32, height = 18, filename = file.path(Dir.PlotNets.PFTC.Iter, "DAG.png"))
  # 
  # dag_data <- na.omit(tidy_dag$data)
  # df <- data.frame(
  #   A = dag_data$name, 
  #   B = dag_data$to)
  # df.g <- graph.data.frame(d = df, directed = TRUE)
  # E(df.g)$weight <- dag_data$weight
  # E(df.g)$label <- dag_data$label
  # E(df.g)$colour <- dag_data$colour
  # 
  # igraph_plot <- ggraph(df.g, layout = 'linear', circular = TRUE) + 
  #   # geom_edge_arc(aes(col = colour, width = weight*WeightMod, alpha = ..index..)) + 
  #   # scale_edge_alpha('Edge direction', guide = 'edge_direction') + 
  #   geom_edge_arc(aes(col = colour, width = exp(dag_data$weight)*WeightMod), show.legend = FALSE) + 
  #   scale_edge_color_manual(values = rev(unique(dag_data$colour))) +
  #   geom_node_point(color = "black", size = 2) + 
  #   geom_node_text(aes(label = name),  repel = TRUE)+
  #   coord_fixed() + 
  #   labs(title = Treatment, subtitle = "Realised Interaction (Bimler et al. Method)") + 
  #   theme_graph()
  # ggsave(igraph_plot, units = "cm", width = 32, height = 32, filename = file.path(Dir.PlotNets.PFTC.Iter, "IGraph.png"))
  
  diag(Interaction_mean) <- 0
  Interaction_mean[abs(sign(Interaction_mean) + sign(Interaction_min) + sign(Interaction_max)) != 3] <- 0
  
  
  
  fv.colors = colorRampPalette(c("red","white","green")) ## define the color ramp
  colorlut = fv.colors(100)[c(1,seq(50,100,length.out=99))] ## select colors to use
  
  seq(sum(Interaction_mean>0))
  
  jpeg(filename = file.path(Dir.PlotNets.PFTC.Iter, "Interactions.jpeg"), units = "cm", width = 32, height = 18, res = 1000)
  pheatmap(Interaction_mean, display_numbers = T, 
           color = c("red", "white","green"), 
           breaks = c(min(Interaction_mean, na.rm = TRUE), -0.01, 0.01, max(Interaction_mean, na.rm = TRUE)), main = "Actors (Columns) x Subject (Rows)",
           fontsize = 5)
  dev.off()
  
  FUN.PlotNetUncert(Model = inter_mat, Dir = Dir.PlotNets.PFTC, Name = Treatment_Iter)
}
