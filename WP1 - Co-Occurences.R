# ####################################################################### #
# PROJECT: [BFTP] Identifying Biomes And Their Shifts Using Remote Sensing
# CONTENTS: Functionality to identify clusters of NDVI mean and seasonality
# AUTHOR: Erik Kusch
# EDIT: 29/05/20
# ####################################################################### #
rm(list=ls())
####### PREAMBLE/SOURCING ---------------------------------------------------------
#### Packages
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}

package_vec <- c(
  "sp", # v1.4-0; for handling spatialpolygondataframes
  "rgdal", # v1.4-8; for loading shapefiles of species ranges
  # "rredlist", # v0.6.0; for loading species id's according to IUCN - NOT NEEDED SINCE I AM USING A MANUAL DOWNLOAD OF IUCN RANGES
  # "BIEN", # v1.2.4; for donwloading plant ranges from BIEN - NOT NEEDED SINCE I AM USING A MANUAL DOWNLOAD OF BIEN RANGES
  "raster", # v3.0-12; for storing data as rasters
  "ncdf4", #v4_1.17; for ncdf namespace when loading nertcdf files
  "fasterize", #v1.0.0; for establishing richness maps in a timely manner
  "sf", #v0.8-1; for use of SpatialpolygonsDataFrame objects in fasterize
  "gimms", #v1.1.1; for downloading the reference raster for NDVI data
  "mapview", #v2.7.0; for creating html maps
  "ggplot2", #v2_3.2.1; for plotting various things
  "rasterVis", #v0.47; for plotting rasters
  "gameofthrones", #v1.0.2; for colour palettes
  "dplyr", #v0.8.4; for data manipulation
  "stringr" #v1.4.0; for padding numbers
)
sapply(package_vec, install.load.package)

#### Directories
Dir.Base <- getwd() # read out the project directory
## DATA
Dir.Data <- file.path(Dir.Base, "3 - Data")
Dir.Ranges <- file.path(Dir.Data, "Ranges")
Dir.Shapes <- file.path(Dir.Data, "Shapes")
Dir.Ranges.Amphibians <- file.path(Dir.Ranges, "AMPHIBIANS")
Dir.Ranges.Reptiles <- file.path(Dir.Ranges, "REPTILES")
Dir.Ranges.Mammals <- file.path(Dir.Ranges, "TERRESTRIAL_MAMMALS")
Dir.Ranges.Birds <- file.path(Dir.Ranges, "BIRDS")
Dir.Ranges.Plants <- file.path(Dir.Ranges, "PLANTS")
DataDirs <- c(Dir.Shapes, Dir.Ranges.Amphibians, Dir.Ranges.Birds, Dir.Ranges.Mammals, Dir.Ranges.Plants, Dir.Ranges.Reptiles)
CreateDir <- sapply(DataDirs, function(x) if(!dir.exists(x)) dir.create(x))
## EXPORTS
Dir.Exports <- file.path(Dir.Base, "2 - Exports")
Dir.Richness <- file.path(Dir.Exports, "1 - Richness")
Dir.Richness.Current <- file.path(Dir.Richness, "1 - Current")
ExportDirs <- c(Dir.Exports, Dir.Richness, Dir.Richness.Current)
CreateDir <- sapply(ExportDirs, function(x) if(!dir.exists(x)) dir.create(x))
rm(list = c("CreateDir", "ExportDirs", "DataDirs"))

####### PREPARATIONS ---------------------------------------------------------
#### REFERENCE RASTER (for projection and resampling)
if(!file.exists(file.path(Dir.Data, "ReferenceRaster.nc"))){ # if reference raster doesn't already exist
  gimms_files <- downloadGimms(x = as.Date("2000-01-01"), # download from January 2000
                               y = as.Date("2000-12-31"), # download to December 2000
                               dsn = Dir.Data, # save downloads in data folder
                               quiet = FALSE) # show download progress
  Reference_ras <- rasterizeGimms(x = gimms_files, remove_header = TRUE) # rasterize the data
  Negatives <- which(values(Reference_ras) < 0) # identify all negative values
  values(Reference_ras)[Negatives] <- 0 # set threshold for barren land (NDVI<0)
  Reference_ras <- calc(Reference_ras, fun=mean, na.rm = TRUE) # annual mean
  writeRaster(x = Reference_ras, filename = file.path(Dir.Data, "ReferenceRaster"), format="CDF") # write raster to data directory
  unlink(gimms_files, recursive = TRUE) # delete raw GIMMs data
  rm(list = c("Negatives", "gimms_files")) # delete from environment
}else{ # if reference raster already exists
  Reference_ras <- raster(file.path(Dir.Data, "ReferenceRaster.nc")) # GIMMs NDVI data as baseline raster  
}

#### LAND MASK (for masking species in the sea which are terrestrial and marine)
if(!file.exists(file.path(Dir.Shapes, "LandMask.zip"))){ # if land mask has not been downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_land.zip", destfile = paste(Dir.Shapes, "LandMask.zip", sep="/")) # download cultural vector
  unzip(paste(Dir.Shapes, "LandMask.zip", sep="/"), exdir = Dir.Shapes) # unzip the data
}
LandMask <- readOGR(Dir.Shapes, "ne_10m_land", verbose = FALSE) # read land mask in

#### COUNTRY MASK (for producing maps with national borders)
if(!file.exists(file.path(Dir.Shapes, "CountryMask.zip"))){ # if land mask has not been downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip", destfile = paste(Dir.Shapes, "CountryMask.zip", sep="/")) # download cultural vector
  unzip(paste(Dir.Shapes, "CountryMask.zip", sep="/"), exdir = Dir.Shapes) # unzip the data
}
CountryMask <- readOGR(Dir.Shapes, "ne_10m_admin_0_countries", verbose = FALSE) # read land mask in

#### LAKE MASK
if(!file.exists(file.path(Dir.Shapes, "LakeMask.zip"))){ # if lake mask has not been downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_lakes.zip", destfile = paste(Dir.Shapes, "LakeMask.zip", sep="/")) # download cultural vector
  unzip(paste(Dir.Shapes, "LakeMask.zip", sep="/"), exdir = Dir.Shapes) # unzip the data
}
LakeMask <- readOGR(Dir.Shapes, "ne_10m_lakes", verbose = FALSE) # read lake mask in

#### RIVER MASK
if(!file.exists(file.path(Dir.Shapes, "RiversMask.zip"))){ # if rivers mask has not been downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_rivers_lake_centerlines_scale_rank.zip", destfile = paste(Dir.Shapes, "RiversMask.zip", sep="/")) # download cultural vector
  unzip(paste(Dir.Shapes, "RiversMask.zip", sep="/"), exdir = Dir.Shapes) # unzip the data
}
RiversMask <- readOGR(Dir.Shapes, "ne_10m_rivers_lake_centerlines_scale_rank", verbose = FALSE) # read river mask in

#### PROTECTED AREAS MASK, this data set is just tooo big!!!
# ### downloaded from https://www.protectedplanet.net/ since download through R doesn't seem to work
# ProtectedAreasMask <- readOGR(file.path(Dir.Shapes, "WDPA_May2020-shapefile"), "WDPA_May2020-shapefile-polygons", verbose = FALSE) # read protected areas mask in

#### ECOREGIONS
if(!file.exists(file.path(Dir.Shapes, "WWF_ecoregions"))){
  download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip", destfile = file.path(Dir.Shapes, "wwf_ecoregions.zip"))
  unzip(file.path(Dir.Shapes, "wwf_ecoregions.zip"), exdir = file.path(Dir.Shapes, "WWF_ecoregions"))
}
EcoregionsMask <- readOGR(file.path(Dir.Shapes, "WWF_ecoregions", "official", "wwf_terr_ecos.shp"), verbose = FALSE) # loading shapefile for wwf ecoregions
### vectors for WWF region naming
Abbr_Realms <- levels((EcoregionsMask$REALM))
Full_Realms <- c("Australasia", "Antarctic", "Afrotropics", "IndoMalay", "Nearctic", "Neotropics", "Oceania", "Palearctic")
Abbr_Biomes <- 1:14
Full_Biomes <- c("Tropical & Subtropical Moist Broadleaf Forests",
                 "Tropical & Subtropical Dry Broadleaf Forests",
                 "Tropical & Subtropical Coniferous Forests",
                 "Temperate Broadleaf & Mixed Forests",
                 "Temperate Conifer Forests",
                 "Boreal Forests/Taiga",
                 "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                 "Temperate Grasslands, Savannas & Shrublands",
                 "Flooded Grasslands & Savannas",
                 "Montane Grasslands & Shrublands",
                 "Tundra",
                 "Mediterranean Forests, Woodlands & Scrub",
                 "Deserts & Xeric Shrublands",
                 "Mangroves")

############## RANGE LOADING -------------------------------------------------
if(!file.exists(file.path(Dir.Ranges, "Ranges_ls.RData")) & !file.exists(file.path(Dir.Ranges, "Rasters_ls.RData"))){
  ## IUCN ----
  Dirs_vec <- c(Dir.Ranges.Amphibians, Dir.Ranges.Reptiles, Dir.Ranges.Mammals) # IUCN directories
  Ranges_ls <- as.list(rep(NA, length(Dirs_vec))) # List for range shapefiles in R
  Rasters_ls <- as.list(rep(NA, length(Dirs_vec))) # List for richness rasters in R
  for(Ranges_Iter in 1:length(Dirs_vec)){ # loop over all IUCN directories
    Name <- strsplit(x = Dirs_vec[Ranges_Iter], split = paste0(Dir.Ranges, "/"))[[1]][2] # isolate name of directory
    Shp <- readOGR(dsn = Dirs_vec[Ranges_Iter], verbose = TRUE) # load range shape data
    Shp <- Shp[which(Shp$terrestial == "true"), ] # mask for terrestrial class
    Ranges_ls[[Ranges_Iter]] <- Shp # save shapefiles to list
    if(!file.exists(file.path(Dir.Ranges, paste0("Richness_", Name, ".nc")))){ # if richness raster has not been established yet
      RichRas <- fasterize(st_as_sf(Shp), raster = Reference_ras, field = NULL, fun = "count") # establish richness raster
      RichRas <- mask(RichRas, LandMask) # maks for land mask 
      RichRas <- mask(RichRas, LakeMask, inverse = TRUE) # mask for lakes
      writeRaster(x = RichRas, filename = file.path(Dir.Ranges, paste0("Richness_", Name)), format="CDF") # write raster with directory name to range directory
    }else{ # if richness raster already exists
      RichRas <-  raster(file.path(Dir.Ranges, paste0("Richness_", Name, ".nc"))) # load richness raster
    }
    Rasters_ls[[Ranges_Iter]] <- RichRas # save richness raster to list
  }
  names(Ranges_ls) <- gsub(pattern = paste0(Dir.Ranges, "/"), replacement = "", x = Dirs_vec) # set names of list
  names(Rasters_ls) <- gsub(pattern = paste0(Dir.Ranges, "/"), replacement = "", x = Dirs_vec) # set names of list
  rm(list = c("Dirs_vec", "Name", "Shp", "RichRas", "Ranges_Iter")) # clean environment
  
  ## BIRDLIFE ----
  Birds_shp <- sf::st_read(dsn = file.path(Dir.Ranges.Birds, "BOTW.gdb")) # load BirdLife data
  classes <- class(Birds_shp$Shape[1])[[1]] # read class of first shape object
  for(Classes_Iter in 2:length(Birds_shp$Shape)){ # loop over all shapes in the shapefile
    classes <- c(classes, class(Birds_shp$Shape[Classes_Iter])[[1]]) # append the class of current object
  }
  Birds_shp <- Birds_shp[which(classes == "sfc_MULTIPOLYGON"), ] # retain MULTIPOLYGON objects, this gets rid of MULTISURFACE objects which fasterize can't handle
  if(!file.exists(file.path(Dir.Ranges, "Richness_BIRDS.nc"))){ # if richness raster has not been established yet
    Birdsrich_ras <- fasterize(Birds_shp, raster = Reference_ras, field = NULL, fun = "count") # establish richness raster
    Birdsrich_ras <- mask(Birdsrich_ras, LandMask) # maks for land mask 
    Birdsrich_ras <- mask(Birdsrich_ras, LakeMask, inverse = TRUE) # mask for lakes
    writeRaster(x = Birdsrich_ras, filename = file.path(Dir.Ranges, "Richness_BIRDS"), format="CDF") # write raster with
  }else{ # if richness raster already exists
    Birdsrich_ras <-  raster(file.path(Dir.Ranges, paste0("Richness_BIRDS.nc"))) # load richness raster
  }
  
  ## BIEN ----
  
  
  ## SAVING TO DISK ----
  Ranges_ls[["BIRDS"]] <- Birds_shp # save to list
  save(Ranges_ls, file = file.path(Dir.Ranges, "Ranges_ls.RData"))
  Rasters_ls[["Total"]] <- sum(stack(Rasters_ls), na.rm = TRUE) # calculate total species richness
  Rasters_ls[["Total"]] <- mask(Rasters_ls[["Total"]], LandMask) # maks for land mask 
  Rasters_ls[["Total"]] <- mask(Rasters_ls[["Total"]], LakeMask, inverse = TRUE) # mask for lakes
  Rasters_ls[["BIRDS"]] <- Birdsrich_ras # save to list
  save(Rasters_ls, file = file.path(Dir.Ranges, "Rasters_ls.RData"))
  rm(list = c("Birds_shp", "Birdsrich_ras", "classes", "Classes_Iter")) # cleaning up environment
}else{
  load(file.path(Dir.Ranges, "Ranges_ls.RData"))
  load(file.path(Dir.Ranges, "Rasters_ls.RData"))
}

############## PLOTTING RICHNESS MAPS ----------------------------------------
## MAPVIEW HTML ----

### PICK IT UP HERE!! - SHOW BIOMES AND REALMS AND PROTECTED AREAS ON THE MAP AS WELL AS RIVERS

Ecoregions_sf <- st_as_sf(EcoregionsMask) 


library(dplyr)

Realms_sf <- Ecoregions_sf %>% 
  group_by(REALM) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()
Realms_sf <- Realms_sf[which(Realms_sf$REALM %in% Abbr_Realms), ]
Realms_sf$REALM <- Full_Realms

Biomes_sf <- Ecoregions_sf %>% 
  group_by(BIOME) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()
Biomes_sf <- Biomes_sf[which(Biomes_sf$BIOME %in% Abbr_Biomes), ]
Biomes_sf$BIOME <- Full_Biomes

# Union_sf <- st_union(Realms_sf, Biomes_sf)
RealmBiomes_sf <- st_intersection(Realms_sf, Biomes_sf)
RealmBiomes_sf$REALMBIOMES <- paste(RealmBiomes_sf$REALM, RealmBiomes_sf$BIOME, sep=" - ")
RealmBiomes_sf <- RealmBiomes_sf[, c(-1,-2)]

# ProtectedAreas_mv <- mapview(ProtectedAreasMask, col.regions = "green", color = "black", alpha.regions = 0.3)
# Countries_mv <- mapview(CountryMask, color = "black", alpha.regions = 0)
Rivers_mv <- mapview(RiversMask, color = "blue", alpha.regions = 0)
# RealmBiomes_mv <- mapview(RealmBiomes_sf, zcol = "REALM", labels = paste(Intersec$REALM, Intersec$BIOME, sep=" - "))
RealmBiomes_mv <- mapview(RealmBiomes_sf, color = "black")
# Biomes_mv <- mapview(Biomes_sf, color = "red", legend = FALSE, alpha.regions = 0)
Amphibians_mv <- mapview(layer.name = "Amphibian Species Richness", Rasters_ls$AMPHIBIANS, legend = TRUE, 
                         maxpixels =  ncell(Reference_ras), na.color = "#FFFFFF00") 
Reptiles_mv <- mapview(layer.name = "Reptilian Species Richness", Rasters_ls$REPTILES, legend = TRUE, 
                       maxpixels =  ncell(Reference_ras), na.color = "#FFFFFF00")
Mammals_mv <- mapview(layer.name = "Mammalian Species Richness", Rasters_ls$TERRESTRIAL_MAMMALS, legend = TRUE, 
                      maxpixels =  ncell(Reference_ras), na.color = "#FFFFFF00")
Birds_mv <- mapview(layer.name = "Avian Species Richness", Rasters_ls$BIRDS, legend = TRUE, 
                    maxpixels =  ncell(Reference_ras), na.color = "#FFFFFF00")
# Plants_mv <- mapview(layer.name = "Plant Species Richness", Rasters_ls$PLANTS, legend = TRUE, 
#                     maxpixels =  ncell(Reference_ras), na.color = "#FFFFFF00")
Total_mv <- mapview(layer.name = "Total Species Richness", Rasters_ls$Total, legend = TRUE, 
                    maxpixels =  ncell(Reference_ras), na.color = "#FFFFFF00")
Richness_mv <- Countries_mv + Amphibians_mv + Reptiles_mv + Mammals_mv + Birds_mv + Total_mv # Combine all maps into 1 map
mapshot(Richness_mv, url = paste0(Dir.Richness.Current, "/RICHNESSMaps.html")) # export to disk
rm(list = c("Countries_mv", "Amphibians_mv", "Reptiles_mv", "Mammals_mv", "Birds_mv", "Total_mv", "Richness_mv"))

## JPEGS ----
col.mapview <- got(n = max(unlist(lapply(X = Rasters_ls, FUN = maxValue))), alpha = 1, begin = 0, end = 1, direction = 1, option = "daenerys")
Title_vec <- c("Amphibian", "Reptilian", "Mammalian", "Avian", 
               # "Plant", 
               "Total")
for(Plot_Iter in 1:length(Rasters_ls)){
  jpeg(file=paste(Dir.Richness.Current, "/", Title_vec[[Plot_Iter]], "_Richness.jpeg", sep = ""), width = 16, height = 12, units = "cm", quality = 100, res = 1000)
  plot(LandMask, col = "grey", bg = "black", main = paste(Title_vec[[Plot_Iter]], "Species Richness"))
  plot(Rasters_ls[[Plot_Iter]], col = col.mapview, add = TRUE, horizontal = TRUE, 
       legend.args = list(text='Number of Species'))
  plot(LakeMask, col = "black", add = TRUE)
  dev.off()
}
rm(list = c("Plot_Iter", "Title_vec", "col.mapview"))

####### SPECIES LIST -----------------------------------------------------
if(!file.exists(file.path(Dir.Ranges, "Species_df.RData"))){
  SpeciesNames_vec <- NA
  DataSets_vec <- NA
  for(Names_Iter in 1:length(Ranges_ls)){ # loop over all range data sets
    if(names(Ranges_ls)[[Names_Iter]] != "BIRDS"){ # BIRD data is stored differently than IUCN data
      SpeciesNamesAdd_vec <- as.character(Ranges_ls[[Names_Iter]]$binomial) # extract species names
    }else{
      SpeciesNamesAdd_vec <- as.character(Ranges_ls[[Names_Iter]]$SCINAME) # extract species names
    }
    SpeciesNames_vec <- c(SpeciesNames_vec, SpeciesNamesAdd_vec) # combine species names
    DataSets_vec <- c(DataSets_vec, rep(Names_Iter, length(SpeciesNamesAdd_vec))) # note range data set (i.e. place in Ranges_ls)
  }
  # Make Species data frame listing names, IDs, and data sets
  Species_df <- data.frame(ID = rep(NA, length(DataSets_vec)),
                           Name = SpeciesNames_vec,
                           Group = DataSets_vec)
  Species_df <- Species_df[-1, ] # remove initial NA row
  Species_df <- Species_df[!duplicated(Species_df), ] # remove multiple mentions of the same species
  Species_df <- Species_df[order(Species_df$Name),] # sort by alphabet
  Species_df$ID <- 1:dim(Species_df)[1] # Species ID for later analysis
  save(Species_df, file = file.path(Dir.Ranges, "Species_df.RData"))
  rm(list = c("DataSets_vec", "SpeciesNames_vec", "Names_Iter"))
}else{
  load(file.path(Dir.Ranges, "Species_df.RData"))
}

####### ESTABLISH CO-OCCURRENCES -----------------------------------------------------
### OCCURRENCE DATA FRAME BIOMES IN REALMS ----
if(!file.exists(file.path(Dir.Ranges, "SpeciesCells_df.RData"))){
  if(!file.exists(file.path(Dir.Ranges, "SpeciesCells_ls.RData"))){
    Realms_ras <- fasterize(sf = st_as_sf(EcoregionsMask), raster = Reference_ras, field = "REALM", fun = "first")
    Biomes_ras <- fasterize(sf = st_as_sf(EcoregionsMask), raster = Reference_ras, field = "BIOME", fun = "first")
    SpeciesCells_ls <- as.list(rep(NA, dim(Species_df)[1])) # a list which will hold sparse data frames for each species and it's CellIDs
    names(SpeciesCells_ls) <- Species_df$Name # set names of the list positions
    Cells_pb <- txtProgressBar(min = 0, max = dim(Species_df)[1], style = 3) # make progress bar
    for(Cells_Iter in Species_df$ID){ # loop over all species IDs
      DataSet <- Ranges_ls[[Species_df$Group[[Cells_Iter]]]] # isolate the data set of the current species
      if(names(Ranges_ls)[[Species_df$Group[Cells_Iter]]] != "BIRDS"){ # BIRD data is stored differently than IUCN data
        Polys <- which(as.character(DataSet$binomial) == as.character(Species_df$Name[[Cells_Iter]])) # identify the polygon position(s)
      }else{
        Polys <- which(as.character(DataSet$SCINAME) == as.character(Species_df$Name[[Cells_Iter]])) # identify the polygon position(s)
      }
      Presence_ras <- fasterize(st_as_sf(DataSet[Polys,]), raster = Reference_ras, field = NULL, fun = "max") # raserize species range
      Cells <- which(!is.na(values(Presence_ras))) # identify cells with data
      # create data frame with species ID and CellIDs
      SpeciesCells_df <- data.frame(SpeciesID = rep(Species_df$ID[[Cells_Iter]], length(Cells)), 
                                    CellID = Cells,
                                    RealmID = values(Realms_ras)[Cells],
                                    BiomeID = values(Biomes_ras)[Cells])
      SpeciesCells_ls[[Cells_Iter]] <- SpeciesCells_df # save dta frame to list
      print(as.character(Species_df$Name[[Cells_Iter]]))
      setTxtProgressBar(Cells_pb, Cells_Iter) # update progress bar
    }
    rm(list = c("SpeciesCells_df", "Cells_pb", "Names_Iter", "Presence_ras", "Polys", "SpeciesNamesAdd_vec", "Cells_Iter", "Cells", "DataSet"))
    save(SpeciesCells_ls, file = file.path(Dir.Ranges, "SpeciesCells_ls.RData"))
  }else{
    load(file.path(Dir.Ranges, "SpeciesCells_ls.RData"))
  }
  SpeciesCells_df <- bind_rows(SpeciesCells_ls, .id = 'column_label')
  rm(SpeciesCells_ls)
  SpeciesCells_df <- na.omit(SpeciesCells_df)
  SpeciesCells_df$RealmID <- str_pad(as.character(SpeciesCells_df$RealmID), 2, "left","0")
  SpeciesCells_df$BiomeID <- str_pad(as.character(SpeciesCells_df$BiomeID), 2, "left","0")
  SpeciesCells_df$BiomesInRealms<- with(SpeciesCells_df, paste(RealmID, BiomeID, sep ="_"))
  save(SpeciesCells_df, file = file.path(Dir.Ranges, "SpeciesCells_df.RData"))
}else{
  load(file.path(Dir.Ranges, "SpeciesCells_df.RData"))
}

### COMMUNITY DATA FRAME BY BIOMES IN REALMS ----
if(!file.exists(file.path(Dir.Ranges, "SpeciesCells_df.RData"))){
  BiomeRealms <- sort(unique(SpeciesCells_df$BiomesInRealms))
  SiteCommunities_ls <- as.list(rep(NA, length(BiomeRealms)))
  names(SiteCommunities_ls) <- BiomeRealms
  for(Realm_Iter in 1:length(BiomeRealms)){
    print(BiomeRealms[Realm_Iter])
    CurrRealm <- BiomeRealms[Realm_Iter]
    CurrRows <- which(SpeciesCells_df$BiomesInRealms == CurrRealm)
    CurrCells <- SpeciesCells_df[CurrRows,]
    SiteCommunities_ls[[Realm_Iter]] <- table(CurrCells[,2:3])
  }
  save(SiteCommunities_ls, file = file.path(Dir.Ranges, "SiteCommunities_ls.RData"))
}else{
  load(file.path(Dir.Ranges, "SiteCommunities_ls.RData"))
}
rm(SpeciesCells_df)

### SPECIES ASSOCIATION NETWORKS ----
# install.packages("netassoc")
library(netassoc)

assoc_net <- make_netassoc_network(obs = SiteCommunities_ls[[2]], verbose = TRUE)
plot_netassoc_network(assoc_net$network_all)



### OUTPUT ----


