library(readr)
library(rlang)
library(reshape2)
library(tidyverse)
library(readxl)

##this function is for compiling auditory physiology assessed by
##the Eaton-Peabody Laboratory Cochlear Function Test Suite, where

##the ABRs have been thresholded and peaks selected with
##ABR Peak Analysis Version 1.5.5.56

##the user initials (CB) and the file path will have to be specified
##in this script.

##Input numbers corresponding to physiology files as a vector (x)

##the current output is a tibble in the R Environment for each data type:
##the ABR threshold (thr), distortion product otoacoustic emmissions (dp), 
##dp isoresponse (isodp), and ABR peak amplitudes and latencies (wav)
##with the range of files indicated on the tibble object (e.g., thr_2299_2408)

##the ABR file number can be used as a key to add experimental info as necessary

compile_phys <- function(x){
  
  ##the file path to auditory physiology data that has been post processed
  path <- "demos/physiology/"
  
  ##initials of the operator
  user <- "CB"
  
  file_num <- expr(!!x)
  
  data_list_abrt <- list()
  data_list_abrw <- list()
  data_list_isodp <- list()
  data_list_dp <- list()
  
  for (i in 1:length(file_num)){
    `wav.5.66` <- read_table2(paste(path, user, file_num[i], "/ABR-", file_num[i], "-1-analyzed.txt", sep=""), 
                              col_names = FALSE, skip = 7) %>%
      mutate(
        freq = 5.66
      )
    `thr.5.66` <- read_delim(paste(path, user, file_num[i], "/ABR-", file_num[i], "-1-analyzed.txt", sep=""), 
                             ":", escape_double = FALSE, col_names = FALSE) %>%
      slice(1:3)
    `wav.8` <- read_table2(paste(path, user, file_num[i], "/ABR-", file_num[i], "-2-analyzed.txt", sep=""), 
                           col_names = FALSE, skip = 7) %>%
      mutate(
        freq = 8
      )
    `thr.8` <- read_delim(paste(path, user, file_num[i], "/ABR-", file_num[i], "-2-analyzed.txt", sep=""), 
                          ":", escape_double = FALSE, col_names = FALSE) %>%
      slice(1:3) 
    `wav.11.3` <- read_table2(paste(path, user, file_num[i], "/ABR-", file_num[i], "-3-analyzed.txt", sep=""), 
                              col_names = FALSE, skip = 7) %>%
      mutate(
        freq = 11.3
      )
    `thr.11.3` <- read_delim(paste(path, user, file_num[i], "/ABR-", file_num[i], "-3-analyzed.txt", sep=""), 
                             ":", escape_double = FALSE, col_names = FALSE) %>%
      slice(1:3) 
    `wav.16` <- read_table2(paste(path, user, file_num[i], "/ABR-", file_num[i], "-4-analyzed.txt", sep=""), 
                            col_names = FALSE, skip = 7) %>%
      mutate(
        freq = 16
      )
    `thr.16` <- read_delim(paste(path, user, file_num[i], "/ABR-", file_num[i], "-4-analyzed.txt", sep=""), 
                           ":", escape_double = FALSE, col_names = FALSE) %>%
      slice(1:3) 
    `wav.22.6` <- read_table2(paste(path, user, file_num[i], "/ABR-", file_num[i], "-5-analyzed.txt", sep=""), 
                              col_names = FALSE, skip = 7) %>%
      mutate(
        freq = 22.6
      )
    `thr.22.6` <- read_delim(paste(path, user, file_num[i], "/ABR-", file_num[i], "-5-analyzed.txt", sep=""), 
                             ":", escape_double = FALSE, col_names = FALSE) %>%
      slice(1:3) 
    `wav.32` <- read_table2(paste(path, user, file_num[i], "/ABR-", file_num[i], "-6-analyzed.txt", sep=""), 
                            col_names = FALSE, skip = 7) %>%
      mutate(
        freq = 32
      )
    `thr.32` <- read_delim(paste(path, user, file_num[i], "/ABR-", file_num[i], "-6-analyzed.txt", sep=""), 
                           ":", escape_double = FALSE, col_names = FALSE) %>%
      slice(1:3) 
    `wav.45.25` <- read_table2(paste(path, user, file_num[i], "/ABR-", file_num[i], "-7-analyzed.txt", sep=""), 
                               col_names = FALSE, skip = 7) %>%
      mutate(
        freq = 45.25
      )
    `thr.45.25` <- read_delim(paste(path, user, file_num[i], "/ABR-", file_num[i], "-7-analyzed.txt", sep=""), 
                              ":", escape_double = FALSE, col_names = FALSE) %>%
      slice(1:3) 
    
    wav <- bind_rows(wav.5.66, wav.8, wav.11.3, wav.16, wav.22.6, wav.32, wav.45.25) %>%
      mutate(
        file = file_num[i]
      ) 
    data_list_abrw[[i]] <- wav
    
    thr <- bind_rows(thr.5.66, thr.8, thr.11.3, thr.16, thr.22.6, thr.32, thr.45.25) %>%
      mutate(
        file = file_num[i]
      ) 
    data_list_abrt[[i]] <- thr
    
    isodp <- read_table2(paste(path, user, file_num[i], "/IsoDP-", file_num[i], "-1", sep=""), 
                         col_names = FALSE, 
                         col_types = cols(X1 = col_double(), X2 = col_double(), 
                                          X3 = col_double(), X4 = col_double(), 
                                          X5 = col_double(), X6 = col_double(), 
                                          X7 = col_double(), X8 = col_double()),
                         skip = 1) %>%
      mutate(
        file = file_num[i]
      )
    data_list_isodp[[i]] <- isodp
    
    dp <- read_table2(paste(path, user, file_num[i], "/DP-", file_num[i], "-1", sep=""), 
                      col_types = cols(`(1)` = col_skip(), 
                                       `2f1-f2Phase` = col_skip(), `PostexposureTime(min)` = col_skip(), 
                                       Shox = col_skip(), `f2-f1(dB)` = col_skip(), 
                                       `f2-f1Nse(dB)` = col_skip(),
                                       `:dB` = col_double(), 
                                       `f1(Hz)` = col_double(), `f2(Hz)` = col_double(), 
                                       `f1(dB)` = col_double(), `f2(dB)` = col_double(), 
                                       `2f1-f2(dB)` = col_double(), `2f1-f2Nse(dB)` = col_double()), skip = 5) %>%
      mutate(
        file = file_num[i]
      )
    dp <- dp[-1,]
    data_list_dp[[i]] <- dp
  }
  
  ##data type compilation
  wave_data <- dplyr::bind_rows(data_list_abrw)
  assign(paste('wav_', head(file_num, n=1), "_", tail(file_num, n=1), sep = ""), wave_data, envir = global_env())
  
  isodp_data <- dplyr::bind_rows(data_list_isodp)
  assign(paste('isodp_', head(file_num, n=1), "_", tail(file_num, n=1), sep = ""), isodp_data, envir = global_env())
  
  dp_data <- dplyr::bind_rows(data_list_dp)
  assign(paste('dp_', head(file_num, n=1), "_", tail(file_num, n=1), sep = ""), dp_data, envir = global_env())
  
  
  ##combine threshold data for data cleaning and compilation
  threshold_data <- dplyr::bind_rows(data_list_abrt)
  
  ##extract the three components from the autothresholding data: threshold, frequency, method
  threshold <- filter(threshold_data, X1 == "Threshold (dB SPL)")
  frequency <- filter(threshold_data, X1 == "Frequency (kHz)")
  method <- filter(threshold_data, X1 == "Threshold estimation")
  
  ##coerce the threshold data to numeric values
  thr <- threshold$X2
  threshold$threshold <- as.numeric(thr)
  
  ##add frequency and method info
  threshold$freq <- frequency$X2
  threshold$method <- method$X2
  
  threshold <- threshold %>% dplyr::select(X2, file, freq, method) %>% rename(threshold = X2)
  
  assign(paste('thr_', head(file_num, n=1), "_", tail(file_num, n=1), sep = ""), threshold, envir = global_env())
}

#compile_phys(37:39)
