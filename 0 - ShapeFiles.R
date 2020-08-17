# ####################################################################### #
# PROJECT: [PhD; Preparations] 
# CONTENTS: Reference raster and Shapefiles
# AUTHOR: Erik Kusch
# EDIT: 17/08/2020
# ####################################################################### #
#### REFERENCE RASTER (for projection and resampling) -----------------------------------------------------------
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

#### LAND MASK (for masking species in the sea which are terrestrial and marine) -----------------------------------------------------------
if(!file.exists(file.path(Dir.Shapes, "LandMask.zip"))){ # if land mask has not been downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_land.zip", destfile = paste(Dir.Shapes, "LandMask.zip", sep="/")) # download cultural vector
  unzip(paste(Dir.Shapes, "LandMask.zip", sep="/"), exdir = Dir.Shapes) # unzip the data
}
LandMask <- readOGR(Dir.Shapes, "ne_10m_land", verbose = FALSE) # read land mask in

#### COUNTRY MASK (for producing maps with national borders) -----------------------------------------------------------
if(!file.exists(file.path(Dir.Shapes, "CountryMask.zip"))){ # if land mask has not been downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip", destfile = paste(Dir.Shapes, "CountryMask.zip", sep="/")) # download cultural vector
  unzip(paste(Dir.Shapes, "CountryMask.zip", sep="/"), exdir = Dir.Shapes) # unzip the data
}
CountryMask <- readOGR(Dir.Shapes, "ne_10m_admin_0_countries", verbose = FALSE) # read land mask in

#### LAKE MASK -----------------------------------------------------------
if(!file.exists(file.path(Dir.Shapes, "LakeMask.zip"))){ # if lake mask has not been downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_lakes.zip", destfile = paste(Dir.Shapes, "LakeMask.zip", sep="/")) # download cultural vector
  unzip(paste(Dir.Shapes, "LakeMask.zip", sep="/"), exdir = Dir.Shapes) # unzip the data
}
LakeMask <- readOGR(Dir.Shapes, "ne_10m_lakes", verbose = FALSE) # read lake mask in

#### RIVER MASK -----------------------------------------------------------
if(!file.exists(file.path(Dir.Shapes, "RiversMask.zip"))){ # if rivers mask has not been downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_rivers_lake_centerlines_scale_rank.zip", destfile = paste(Dir.Shapes, "RiversMask.zip", sep="/")) # download cultural vector
  unzip(paste(Dir.Shapes, "RiversMask.zip", sep="/"), exdir = Dir.Shapes) # unzip the data
}
RiversMask <- readOGR(Dir.Shapes, "ne_10m_rivers_lake_centerlines_scale_rank", verbose = FALSE) # read river mask in

#### PROTECTED AREAS MASK, this data set is just tooo big!!!
# ### downloaded from https://www.protectedplanet.net/ since download through R doesn't seem to work
# ProtectedAreasMask <- readOGR(file.path(Dir.Shapes, "WDPA_May2020-shapefile"), "WDPA_May2020-shapefile-polygons", verbose = FALSE) # read protected areas mask in

#### ECOREGIONS -----------------------------------------------------------
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