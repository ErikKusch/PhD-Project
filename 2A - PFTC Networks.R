# ####################################################################### #
# PROJECT: [PhD; 2A - PFTC NETWORKS]
# CONTENTS: Generate PFTC species-association/-interaction networks
# AUTHOR: Erik Kusch
# EDIT: 30/06/2021
# ####################################################################### #

rm(list=ls())
####### PREAMBLE + SHAPES + REFERENCES ---------------------------------------------------------
source("0 - Preamble.R")
# source("0 - ShapeFiles.R")
source("X - Functions_Plotting.R")
source("X - Functions_Bayes.R")
nSamples <- 7000
nWarmup <- 700
nChains <- 4

####### PFTC METADATA ---------------------------------------------------------
if(file.exists(file.path(Dir.PFTC, "Metadata.csv"))){
  Metadata_df <- read.csv(file.path(Dir.PFTC, "Metadata.csv"))
}else{
  Metadata_df <- data.frame(Site = rep(c("WAY", "TRE", "ACJ", "PIL", "QUE", "OCC"), 
                                       c(10, 10, 10, 13, 5, 5)),
                            Plot = c(rep(1:5, 8), 1:3, rep(1:5, 2)),
                            Treatment = c(rep(c("B", "C", "C", "BB", "B", "C", "B", "C"), each = 5), "BB", "BB", "BB", rep(c("B", "C"), each = 5)),
                            BurnYear = c(rep(c(2007, NA, NA, 2018, 2005, NA, 2007, NA), each = 5), 2013, 2013, 2013, rep(c(NA, NA), each = 5)), 
                            Lat = c(-13.179917, -13.179978, -13.180036, -13.180092, NA, -13.180406, -13.180430, -13.180457, -13.180527, -13.180580, -13.119746, -13.119733, -13.119713, -13.119695, -13.119707, -13.138776, -13.138797, -13.138823, -13.138837, -13.138871, -13.169344, -13.169347, -13.169447, -13.169483, -13.169548, -13.176034, -13.176100, -13.176134, -13.176195, -13.176252, -13.151583, -13.151576, -13.151582, -13.151584, -13.151593, -13.151610, -13.151633, -13.151780, -13.151791, -13.151803, -13.151819, -13.151832, -13.151323, -13.213786, -13.213841, -13.214201, -13.214146, -13.214258, -13.450887, -13.450843, -13.450800, -13.450756, -13.450745),
                            Lon = c(-71.588118, -71.588188, -71.588263, -71.588355, NA, -71.589709, -71.589752, -71.589806, -71.589873, -71.589928, -71.620771, -71.620810, -71.620874, -71.620927, -71.620930, NA, NA, NA, NA, NA, -71.634001, -71.634004, -71.634022, -71.634032, -71.634030, -71.628828, -71.628834, -71.628834, -71.628844, -71.628819, -71.640710, -71.640552, -71.640513, -71.640511, -71.640514, -71.640462, -71.640395, -71.640027, -71.639945, -71.639885, -71.639813, -71.639748, -71.640354, -71.619294, -71.619348, -71.619384, -71.619396, -71.619377, -71.741101, -71.741050, -71.741017, -71.740969, -71.740922),
                            Elevation = c(3071.7, 3073.1, 3077.1, 3078.8, NA, 3118.0, 3120.3, 3121.0, 3123.3, 3125.0, 3714.3, 3714.8, 3714.7, 3714.7, 3714.6, 3618.2, 3619.1, 3619.4, 3620.3, 3621.3, 3444.1, 3444.2, 3446.3, 3447.3, 3447.7, 3487.4, 3487.7, 3489.0, 3490.0, 3492.4, 3626.7, 3687.4, 3687.0, 3686.9, 3686.8, 3685.5, 3683.5, 3673.3, 3671.2, 3669.6, 3665.4, 3665.5, 3694.3, 3882.9, 3884.3, 3890.0, 3889.8, 3893.1, 4384.3, 4382.4, 4380.0, 4382.0, 4385.8)
  )
  write.csv(Metadata_df, file.path(Dir.PFTC, "Metadata.csv"))
}

####### PFTC RAW DATA ---------------------------------------------------------
## load and limit data
# Raw_df <- read.csv(file.path(Dir.PFTC, "PFTC3.1_CommunityCover_2018_Peru.csv")) # cover of species
Raw_df <- read.csv(file.path(Dir.PFTC, "PFTC3.7_Traits_2018_Peru_cleaned.csv")) # trait expressions
Weights_df <- Raw_df[,c("Site", "Treatment", "PlotID", "Taxon", "Wet_Mass_Total_g", "Individual_nr")]
## create backbone for intrinsic fitness data frame
Mal_df <- with(Weights_df, data.frame(
  PlotID = paste(Site, Treatment, PlotID, sep = "_"),
  fit = Wet_Mass_Total_g,
  focal = Taxon,
  focalID = paste(Taxon, Individual_nr, Site, Treatment, PlotID, sep = "_")
  )
)
## create abundance data (to be appended to backbone data frame)
Abund_df <- aggregate(Individual_nr~Site+Treatment+PlotID+Taxon, data = Weights_df, FUN = max)
Abund_df$PlotID <- with(Abund_df, paste(Site, Treatment, PlotID, sep = "_"))
## create columns to hold abundance of species at each plot
Mal_df <- cbind(Mal_df, matrix(rep(NA, dim(Mal_df)[1] * length(unique(Abund_df$Taxon))), nrow = dim(Mal_df)[1]))
colnames(Mal_df)[-1:-4] <- as.character(unique(Abund_df$Taxon))
## add abundances to each plot for each species
for(i in 1:nrow(Abund_df)){
  PlotRows <- which(Mal_df$PlotID == Abund_df$PlotID[i])
  FocalCol <- which(colnames(Mal_df) == Abund_df$Taxon[i])
  Mal_df[PlotRows,FocalCol] <- Abund_df$Individual_nr[i]
}

################# RUN NETWORK METHOD  ########################################### --------------------------
Treatments_df <- data.frame(Treatment = sapply(str_split(Mal_df$PlotID, "_"), "[[", 2),
                            Plot = sapply(str_split(Mal_df$PlotID, "_"), "[[", 1)
)
Treatments_vec <- c("ALL", as.character(unique(Treatments_df$Treatment)), as.character(unique(Treatments_df$Plot)))           
Treatments_vec <- Treatments_vec[Treatments_vec != "NA"]
Mal_df[is.na(Mal_df)] <- 0

for(Treatment in Treatments_vec){
  Run_df <- Mal_df
  if(Treatment != "ALL"){
    Run_df <- Run_df[which(Treatments_df[,1] == Treatment | Treatments_df[,2] == Treatment), ]
  }
  
  Dir.PlotNets.PFTC.Iter <- file.path(Dir.PlotNets.PFTC, Treatment)
  if(dir.exists(Dir.PlotNets.PFTC.Iter)){
    print(paste("PFTC models already run for Treatment/Plot, ", Treatment,". You can find them here:", Dir.PlotNets.PFTC.Iter))
    next()
  }
  dir.create(Dir.PlotNets.PFTC.Iter)
  OmitCols <- -(which(colnames(Mal_df)[-c(1:4)] %nin% Mal_df$focal)+4)
  if(length(OmitCols) != 0){
    Run_df <- Run_df[, OmitCols]
  }else{
    Run_df <- Run_df
  }
  Run_df <- Run_df[Run_df$fit != 0, ]
  
  FIA_StanList <- Fun_StanList(Fitness = "fit", data = Run_df)
  
  #### Preferences
  rstan_options(auto_write = TRUE)
  rstan_options(javascript = FALSE)
  options(mc.cores = nChains) 
  
  ## CHECKING DATA
  focalID <- unique(Run_df$focal)
  neighbourID <- colnames(Run_df[ , -c(1:4)])
  Fun_PreCheck(data = Run_df)
  ## RUNNING MODEL
  fit <- stan(file = 'Supplement - StanModel.stan',
              data =  FIA_StanList,               # named list of data
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
  
  PlotNetUncert(Model = inter_mat, Dir = Dir.PlotNets.PFTC, Name = Treatment)
}
