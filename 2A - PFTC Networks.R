# ####################################################################### #
# PROJECT: [PhD; 2A - PFTC NETWORKS] 
# CONTENTS: Generate PFTC species-association/-interaction networks
# AUTHOR: Erik Kusch
# EDIT: 18/08/2020
# ####################################################################### #
####### DATA LOADING ---------------------------------------------------------
if(!file.exists(file.path(Dir.PFTC, "PFTC3Cover.csv"))){
  Raw_df <- read.csv(file.path(Dir.PFTC, "PFTC3.1_CommunityCover_2018_Peru.csv"))
}else{
  Cover_df <- read.csv(file.path(Dir.PFTC, "PFTC3Cover.csv"))
  Cover_df <- Cover_df[, -1]
  colnames(Cover_df) <- gsub(pattern = "\\.", replacement = " ", x = colnames(Cover_df))
}
if(!file.exists(file.path(Dir.PFTC, "PFTC3Traits.csv"))){
  TRaw_df <- read.csv(file.path(Dir.PFTC, "PFTC3.7_Traits_2018_Peru_cleaned.csv"))
}else{
  Traits_df <- read.csv(file.path(Dir.PFTC, "PFTC3Cover.csv"))
  Traits_df <- Traits_df[, -1]
}
 

####### FIXING PFTC DATA ---------------------------------------------------------
if(exists("Raw_df")){
  Stations_df <- data.frame(Site = rep(c("WAY", "TRE", "ACJ", "PIL", "QUE", "OCC"), 
                                       c(10, 10, 10, 13, 5, 5)),
                            Plot = c(rep(1:5, 8), 1:3, rep(1:5, 2)),
                            Treatment = c(rep(c("B", "C", "C", "BB", "B", "C", "B", "C"), each = 5), "BB", "BB", "BB", rep(c("B", "C"), each = 5)),
                            BurnYear = c(rep(c(2007, NA, NA, 2018, 2005, NA, 2007, NA), each = 5), 2013, 2013, 2013, rep(c(NA, NA), each = 5)), 
                            Lat = c(-13.179917, -13.179978, -13.180036, -13.180092, NA, -13.180406, -13.180430, -13.180457, -13.180527, -13.180580, -13.119746, -13.119733, -13.119713, -13.119695, -13.119707, -13.138776, -13.138797, -13.138823, -13.138837, -13.138871, -13.169344, -13.169347, -13.169447, -13.169483, -13.169548, -13.176034, -13.176100, -13.176134, -13.176195, -13.176252, -13.151583, -13.151576, -13.151582, -13.151584, -13.151593, -13.151610, -13.151633, -13.151780, -13.151791, -13.151803, -13.151819, -13.151832, -13.151323, -13.213786, -13.213841, -13.214201, -13.214146, -13.214258, -13.450887, -13.450843, -13.450800, -13.450756, -13.450745),
                            Lon = c(-71.588118, -71.588188, -71.588263, -71.588355, NA, -71.589709, -71.589752, -71.589806, -71.589873, -71.589928, -71.620771, -71.620810, -71.620874, -71.620927, -71.620930, NA, NA, NA, NA, NA, -71.634001, -71.634004, -71.634022, -71.634032, -71.634030, -71.628828, -71.628834, -71.628834, -71.628844, -71.628819, -71.640710, -71.640552, -71.640513, -71.640511, -71.640514, -71.640462, -71.640395, -71.640027, -71.639945, -71.639885, -71.639813, -71.639748, -71.640354, -71.619294, -71.619348, -71.619384, -71.619396, -71.619377, -71.741101, -71.741050, -71.741017, -71.740969, -71.740922),
                            Elevation = c(3071.7, 3073.1, 3077.1, 3078.8, NA, 3118.0, 3120.3, 3121.0, 3123.3, 3125.0, 3714.3, 3714.8, 3714.7, 3714.7, 3714.6, 3618.2, 3619.1, 3619.4, 3620.3, 3621.3, 3444.1, 3444.2, 3446.3, 3447.3, 3447.7, 3487.4, 3487.7, 3489.0, 3490.0, 3492.4, 3626.7, 3687.4, 3687.0, 3686.9, 3686.8, 3685.5, 3683.5, 3673.3, 3671.2, 3669.6, 3665.4, 3665.5, 3694.3, 3882.9, 3884.3, 3890.0, 3889.8, 3893.1, 4384.3, 4382.4, 4380.0, 4382.0, 4385.8)
  )
  Stations_df <- cbind(Stations_df, matrix(rep(NA, dim(Stations_df)[1] * length(unique(Raw_df$Taxon))), nrow = dim(Stations_df)[1]))
  colnames(Stations_df)[8:dim(Stations_df)[2]] <- sort(as.character(unique(Raw_df$Taxon)))
  TotalSites <- unique(Stations_df$Site)[unique(Stations_df$Site) %in% unique(Raw_df$Site)]
  Counter <- 1 
  for(Iter_Site in 1:length(TotalSites)){
    SiteIter_df <- Raw_df[which(as.character(Raw_df$Site) == as.character(TotalSites[Iter_Site])), ]
    SiteIter_df2 <- Stations_df[which(as.character(Stations_df$Site) == as.character(TotalSites[Iter_Site])), ]
    TotalTreatments <- unique(SiteIter_df2$Treatment)
    for(Iter_Treatment in 1:length(TotalTreatments)){
      TreatmentIter_df <- SiteIter_df[which(SiteIter_df$Treatment == TotalTreatments[Iter_Treatment]), ]
      TreatmentIter_df2 <- SiteIter_df2[which(SiteIter_df2$Treatment == TotalTreatments[Iter_Treatment]), ]
      TotalPlots <- unique(TreatmentIter_df2$Plot)
      for(Iter_Plot in TotalPlots){
        PlotIter_df <- TreatmentIter_df[which(TreatmentIter_df$PlotID == Iter_Plot), ]
        try(PlotIter_df <- aggregate(PlotIter_df$Cover, list(PlotIter_df$Taxon), mean), silent = TRUE)
        colnames(PlotIter_df) <- c("Taxon", "Cover")
        Stations_df[Counter, which(colnames(Stations_df) %in% as.character(PlotIter_df$Taxon))] <- PlotIter_df$Cover
        Counter <- 1 + Counter  
      }
    }
  }
  Discard_rows <- -which(Stations_df$Site %in% unique(Stations_df$Site)[which(unique(Stations_df$Site) %in% unique(Raw_df$Site) == FALSE)])
  Cover_df <- Stations_df[Discard_rows, ]
  rm(TreatmentIter_df, TreatmentIter_df2, SiteIter_df, SiteIter_df2, PlotIter_df, Raw_df, Stations_df, Iter_Plot, Iter_Site, Iter_Treatment, TotalPlots, TotalSites, TotalTreatments, Discard_rows, Counter)
  write.csv(Cover_df, file.path(Dir.PFTC, "PFTC3Cover.csv"))
}

if(exists("TRaw_df")){
Traits_df <- TRaw_df[, c(5,6,11,12,17:20)]
Traits_df <- aggregate(x = Traits_df[-c(1:4)], by = list(Traits_df$Taxon, Traits_df$Site, Traits_df$PlotID, Traits_df$Treatment), mean, na.rm = TRUE)
colnames(Traits_df)[1:4] <- c("Species", "Site", "Plot", "Treatment")
rm(TRaw_df)
write.csv(Traits_df, file.path(Dir.PFTC, "PFTC3Traits.csv"))
}

Cover_df <- Cover_df[, c(1:6, 6+which(colnames(Cover_df)[-c(1:7)] %in% Traits_df$Species))]

####### BUILDING NETWORKS ---------------------------------------------------------

Cover_df$SitePlot <- paste(Cover_df$Site, Cover_df$PlotID, sep="_")

Cover_df$Cover <- Cover_df$Cover * 100
# Sites_ls <- split(Cover_df, Cover_df$Site)
# Site_df <- Sites_ls[[1]]


# library(reshape)
# Community_df <- as.data.frame(cast(Site_df, PlotID ~ Taxon, value='Cover', fun.aggregate = mean))
# rownames(Community_df) <- Community_df[, 1]
# Community_df <- Community_df[, -1]
# Community_df[is.na(Community_df)] <- 0


Sites_vec <- unique(Cover_df$Site)
Species_vec <- unique(Cover_df$Taxon)

Community_df <- data.frame(Site = NA, Species = NA)

for(Iter_Site in 1:length(Sites_vec)){
  Site_rows <- which(Cover_df$Site == Sites_vec[Iter_Site])
  Site_species <- Cover_df$Taxon[Site_rows]
  for(Iter_Species in 1:length(Site_species)){
    SpeciesSite_rows <- which(Cover_df$Taxon[Site_rows] == Site_species[Iter_Species])
    Species <- rep(as.character(Site_species[Iter_Species]), mean(Cover_df$Cover[SpeciesSite_rows]))
    Site <- rep(Sites_vec[Iter_Site], length(Species))
    Community_df_append <- data.frame(Site = Species, Species = Site)
    Community_df <- rbind(Community_df, Community_df_append)
  }
}
Network_df <- table(Community_df)


library(netassoc)
assoc_net <- make_netassoc_network(obs = Network_df)
plot_netassoc_network(assoc_net$network_all)
