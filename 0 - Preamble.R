# ####################################################################### #
# PROJECT: [PhD; Preparations] 
# CONTENTS: Packages and Directories
# AUTHOR: Erik Kusch
# EDIT: 17/08/2020
# ####################################################################### #
####### PACKAGES ---------------------------------------------------------
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

####### DIRECTORIES ---------------------------------------------------------
Dir.Base <- getwd() # read out the project directory
## DATA
Dir.Data <- file.path(Dir.Base, "3 - Data")
Dir.PFTC <- file.path(Dir.Data, "PFTC")
Dir.Ranges <- file.path(Dir.Data, "Ranges")
Dir.Shapes <- file.path(Dir.Data, "Shapes")
Dir.Ranges.Amphibians <- file.path(Dir.Ranges, "AMPHIBIANS")
Dir.Ranges.Reptiles <- file.path(Dir.Ranges, "REPTILES")
Dir.Ranges.Mammals <- file.path(Dir.Ranges, "TERRESTRIAL_MAMMALS")
Dir.Ranges.Birds <- file.path(Dir.Ranges, "BIRDS")
Dir.Ranges.Plants <- file.path(Dir.Ranges, "PLANTS")
DataDirs <- c(Dir.Shapes, Dir.PFTC, Dir.Ranges.Amphibians, Dir.Ranges.Birds, Dir.Ranges.Mammals, Dir.Ranges.Plants, Dir.Ranges.Reptiles)
CreateDir <- sapply(DataDirs, function(x) if(!dir.exists(x)) dir.create(x))
## EXPORTS
Dir.Exports <- file.path(Dir.Base, "2 - Exports")
Dir.Richness <- file.path(Dir.Exports, "1 - Richness")
Dir.Richness.Current <- file.path(Dir.Richness, "1 - Current")
ExportDirs <- c(Dir.Exports, Dir.Richness, Dir.Richness.Current)
CreateDir <- sapply(ExportDirs, function(x) if(!dir.exists(x)) dir.create(x))
rm(list = c("CreateDir", "ExportDirs", "DataDirs"))