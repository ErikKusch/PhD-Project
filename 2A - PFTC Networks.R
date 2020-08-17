# ####################################################################### #
# PROJECT: [PhD; 2A - PFTC NETWORKS] 
# CONTENTS: Generate PFTC species-association/-interaction networks
# AUTHOR: Erik Kusch
# EDIT: 17/08/2020
# ####################################################################### #
####### DATA LOADING ---------------------------------------------------------
raster(file.path(Dir.PFTC, "PFTC3.1_CommunityCover_2018_Peru.csv"))
raster(file.path(Dir.PFTC, "PFTC3.7_Traits_2018_Peru_cleaned.csv"))