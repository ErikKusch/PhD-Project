# install.packages("rFIA")


library("rFIA")



FIA <- getFIA(states = "VT")

WP_df <- biomass(db = FIA, # which data base to use
                 bySpecies = TRUE, # group by SPecies
                 byPlot = TRUE, # group by plot
                 nCores = 4
)
View(WP_df)

### plot numbers with geo-coordinates?!











