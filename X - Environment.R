#' ####################################################################### #
#' PROJECT: [PhD; X - ENVIRONMENT] 
#' CONTENTS: 
#'  - Functionality for retrieval of environmental data for observations
#'  DEPENDENCIES:
#'  - 
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #

FUN.ClimData <- function(Data = Metadata_df,
                         ID = "ID",
                         Lat = "Latitude",
                         Lon = "Longitude",
                         Year = "Year",
                         Cores = 1,
                         Shape = NULL,
                         FileName = "Metadata_df", 
                         force = FALSE, 
                         rawdata = FALSE, 
                         Dir = getwd()){
    
  # Download Identification ----
    Data$Download <- 2 # create download query column, 2 indicates insufficient metadata data
    Data$Download[rowSums(!is.na(Data[,c(Lat,Lon, Year)])) == 3] <- 1 # 1 indicates sufficient metadata availability
    Data$Download[which(as.numeric(as.character(Data[,Year])) < 1981 & Data$Download == 1)] <- 3 # 3 indicates study end date predating first available year of era5-land data availability
    if(!is.null(Shape)){
      DownloadQueries <- which(Data$Download == 1)
      CoordinateCheck_df <- Data[DownloadQueries,4:5]
      coordinates(CoordinateCheck_df) <- ~Lon+Lat
      proj4string(CoordinateCheck_df) <- proj4string(Shape)
      OverCheck_df <- sp::over(x = CoordinateCheck_df, y = Shape)
      OutofBounds <- which(is.na(OverCheck_df$featurecla)) ## these are the points which have download query 1, but fall outside of era5-land shape, currently only takes naturalearth land mask data
      Data$Download[DownloadQueries][OutofBounds] <- 4 # 4 indicates locations outside of land mask
    }

    # File Check ----
    if(isTRUE(force)){ # identify which of the entires already have the necessary data
      DownsNeeded <- Data[,ID]
    }else{
      DownsNeeded <- Data[, ID][is.na(Data$BIO1) & Data$Download == 1]
    }
    # DOWNLOAD LOOP ----
    if(length(unique(Data[Data[,ID] %in% DownsNeeded,Year])) == 1){
      ## data
      Y_1 <- 1981
      Y_2 <- unique(Data[Data[,ID] %in% DownsNeeded,Year])
      Down_df <- Data[Data[,ID] %in% DownsNeeded,]
      ## download
      BC_ras <- BioClim(
        Water_Var = "total_precipitation",
        Y_start = Y_1,
        Y_end = Y_2,
        T_res = "day",
        Extent = Down_df,
        DataSet = "era5-land",
        Dir = Dir,
        verbose = TRUE,
        FileName = FileName,
        Keep_Raw = rawdata,
        Keep_Monthly = FALSE,
        Buffer = 0.1,
        ID = ID,
        API_User = API_User,
        API_Key = API_Key,
        Cores = Cores
      )
      
      ###CONTINUE FROM HERE!!!!!!
      ## data extraction and adding to metadata dataframe
      BC_df <- raster::extract(x = BC_ras, y = data.frame(Down_df[, Lon], Down_df[, Lat]), method = "bilinear")
      Data[Data[, ID] %in%  DownsNeeded, which(colnames(Data) == "BIO1"):which(colnames(Data) == "BIO19")] <- BC_df
      ## writing metadata frame
      write.csv(Data, file = file.path(Dir,  paste0(FileName, ".csv")))
    }else{
      for(Down_Iter in DownsNeeded){
        ## data
        Iter_df <- Data[Data[, ID] == Down_Iter, ]
        Y_1 <- 1981
        Y_2 <- as.numeric(as.character(Iter_df[,Year]))
        ## directory
        Dir.Iter <- file.path(Dir, Down_Iter)
        dir.create(Dir.Iter)
        ## bioclimatic
        BC_ras <- BioClim(
          Water_Var = "total_precipitation",
          Y_start = Y_1,
          Y_end = Y_2,
          T_res = "day",
          Extent = Iter_df,
          DataSet = "era5-land",
          Dir = Dir.Iter,
          verbose = TRUE,
          FileName = Down_Iter,
          Keep_Raw = rawdata,
          Keep_Monthly = FALSE,
          Buffer = 0.1,
          ID = ID,
          API_User = API_User,
          API_Key = API_Key,
          Cores = Cores
        )

        ## compilation of raw data if queried
        if(isTRUE(rawdata)){
          ## loading data and extracting values
          setwd(Dir.Iter)
          Temps_ls <- list.files(pattern = "2m_temperature_Temporary")
          Temps_vec <- extract(x = stack(Temps_ls), y = data.frame(Iter_df[, Lon], Iter_df[, Lat]), method = "bilinear")
          Precips_ls <- list.files(pattern = "total_precipitation_Temporary")
          Precips_vec <- extract(x = stack(Precips_ls), y = data.frame(Iter_df[, Lon], Iter_df[, Lat]), method = "bilinear")
          ## figuring out days
          Down_start <- lubridate::date(paste0(Y_1, "-01-01"))
          Down_end <- lubridate::date(paste0(Y_2, "-12-31"))
          T_seq <- seq(Down_start, Down_end, by = "day")
          ## data frame and saving
          Raw_df <- data.frame(
            Date = T_seq,
            Temperature = as.numeric(Temps_vec),
            Precipitation = as.numeric(Precips_vec)
          )
          write.csv(x = Raw_df, file = file.path(Dir, paste0(Down_Iter, "_Raw.csv")))
          rm("Temps_ls", "Temps_vec", "Precips_ls", "Precips_vec")
          setwd(Dir)
        }
        
        ## data extraction and adding to metadata dataframe
        BC_vec <- extract(x = BC_ras, y = data.frame(Iter_df[, Lon], Iter_df[, Lat]), method = "bilinear")
        Data[Data[, ID] == Down_Iter, which(colnames(Data) == "BIO1"):which(colnames(Data) == "BIO19")] <- BC_vec
        ## writing metadata frame
        write.csv(Data, file = file.path(Dir,  paste0(FileName, ".csv")))
    
        ## cleaning of local files
        unlink(Dir.Iter, recursive = TRUE)
      } 
    }
    return(Data)
  }