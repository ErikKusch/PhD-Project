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
  "sp", # v1.4-2; for handling spatialpolygondataframes
  "rgdal", # v1.5-16; for loading shapefiles of species ranges
  # "rredlist", # v0.6.0; for loading species id's according to IUCN - NOT NEEDED SINCE I AM USING A MANUAL DOWNLOAD OF IUCN RANGES
  # "BIEN", # v1.2.4; for donwloading plant ranges from BIEN - NOT NEEDED SINCE I AM USING A MANUAL DOWNLOAD OF BIEN RANGES
  "raster", # v3.3-13; for storing data as rasters
  "ncdf4", #v4_1.17; for ncdf namespace when loading nertcdf files
  "fasterize", #v1.0.3; for establishing richness maps in a timely manner
  "sf", #v0.9-6; for use of SpatialpolygonsDataFrame objects in fasterize
  "gimms", #v1.1.3; for downloading the reference raster for NDVI data
  "mapview", #v2.9.0; for creating html maps
  "ggplot2", #v2_3.3.2; for plotting various things
  "rasterVis", #v0.48; for plotting rasters
  "gameofthrones", #v1.0.2; for colour palettes
  "dplyr", #v1.0.2; for data manipulation
  "stringr", #v1.4.0; for padding numbers
  "rFIA", #0.2.5; for downloading and using Forest Inventory Analysis (FIA) data
  "rethinking", #2.13; for Bayesian methdology and model inspection
  "rstan", #2.12.2; for access to stan
  "dagitty", #0.3-0; for DAG drawing
  "ggdag", # 0.2.3; for DAG drawing in ggplot environment
  "data.table", #1.12.8; for DAG building
  "igraph", #1.2.5; for nicer graph visualisation
  "pheatmap", #1.0.12; for heatmaps of interactions
  "ggraph" #2.0.5 ; for igraph ggplot layouts
)
sapply(package_vec, install.load.package)

####### DIRECTORIES ---------------------------------------------------------
Dir.Base <- getwd() # read out the project directory
## DATA
Dir.Data <- file.path(Dir.Base, "3 - Data")
Dir.PFTC <- file.path(Dir.Data, "PFTC")
Dir.Plots <- file.path(Dir.Data, "Plots")
Dir.Plots.FIA <- file.path(Dir.Plots, "FIA")
Dir.Ranges <- file.path(Dir.Data, "Ranges")
Dir.Shapes <- file.path(Dir.Data, "Shapes")
Dir.Ranges.Amphibians <- file.path(Dir.Ranges, "AMPHIBIANS")
Dir.Ranges.Reptiles <- file.path(Dir.Ranges, "REPTILES")
Dir.Ranges.Mammals <- file.path(Dir.Ranges, "TERRESTRIAL_MAMMALS")
Dir.Ranges.Birds <- file.path(Dir.Ranges, "BIRDS")
Dir.Ranges.Plants <- file.path(Dir.Ranges, "PLANTS")
DataDirs <- c(Dir.Data, Dir.Ranges, Dir.Shapes, Dir.PFTC, Dir.Ranges.Amphibians, Dir.Ranges.Birds, Dir.Ranges.Mammals, Dir.Ranges.Plants, Dir.Ranges.Reptiles, Dir.Plots, Dir.Plots.FIA)
CreateDir <- sapply(DataDirs, function(x) if(!dir.exists(x)) dir.create(x))
## EXPORTS
Dir.Exports <- file.path(Dir.Base, "2 - Exports")
Dir.PlotNets <- file.path(Dir.Exports, "2 - Plot Networks")
Dir.PlotNets.FIA <- file.path(Dir.PlotNets, "FIA")
Dir.PlotNets.PFTC <- file.path(Dir.PlotNets, "PFTC")
Dir.Richness <- file.path(Dir.Exports, "X - Richness")
Dir.Richness.Current <- file.path(Dir.Richness, "1 - Current")
ExportDirs <- c(Dir.Exports, Dir.PlotNets, Dir.PlotNets.FIA, Dir.PlotNets.PFTC, Dir.Richness, Dir.Richness.Current)
CreateDir <- sapply(ExportDirs, function(x) if(!dir.exists(x)) dir.create(x))
rm(list = c("CreateDir", "ExportDirs", "DataDirs"))


####### FUNCTIONALITY ---------------------------------------------------------
`%nin%` <- Negate(`%in%`)
