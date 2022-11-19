library(tidyverse)
library(readxl)
library(reshape2)

##this function takes a tibble/data frame as input, where:
##the first column is 5-min intervals (time in 5 min bins)
##the rest of the columns contain actigraphy measures for each participant

compile_inactivity_data <- function(x) {
  
  cluster <- list()
  data_cluster <- list()
  
  inactivity_data <- melt(expr(!!x), id = "5-min intervals", variable.name = "indiv")
  
  inactivity_data <- inactivity_data %>% group_by(indiv) %>%
    mutate(max = max(value),
           threshold = max*0.1,
           tf_val = ifelse(value<threshold, TRUE, FALSE),
           indiv = as.character(indiv))
  
  names <-  unique(inactivity_data$indiv)
  
  for (i in 1:length(names)) {
    
    ##only TRUE observations
    sample_tb <- inactivity_data %>% filter(tf_val==TRUE)
    
    sample_tb2 <- sample_tb %>% filter(indiv==names[i])
    print(unique(sample_tb2$indiv))
    
    sample_tb2$size <- NA
    sample_tb2$len <- 1:length(sample_tb2$`5-min intervals`)
    sample_tb2$series <- sample_tb2$`5-min intervals`
    ##assign a size value to consecutive observations
    
    for (i in 1:length(sample_tb2$series)) {
      
      print(sample_tb2$series[i] + 1)
      print(sample_tb2$series[i+1])
      
      try(if ((sample_tb2$series[i] + 5) != sample_tb2$series[i+1]) {
        sample_tb2[i,]$size <- 1
      } else if ((sample_tb2$series[i] + 10) != sample_tb2$series[i+2] | is.na(sample_tb2$series[i+2])) {
        sample_tb2[i,]$size <- 2
      } else if ((sample_tb2$series[i] + 15) != sample_tb2$series[i+3] | is.na(sample_tb2$series[i+3])) {
        sample_tb2[i,]$size <- 3
      } else if ((sample_tb2$series[i] + 20) != sample_tb2$series[i+4] | is.na(sample_tb2$series[i+4])) {
        sample_tb2[i,]$size <- 4
      } else
        sample_tb2[i,]$size <- 4, silent = TRUE)
      
    }
    ##make the last item in the size vector equal to 1
    sample_tb2[nrow(sample_tb2),]$size <- 1
    
    ##the items which are size 1 indicate the end of clusters or lone observations
    sample_tb2_1 <- sample_tb2 %>% filter(size==1)
    
    ##use the size = 1 observations to assign clusters
    
    cluster <- list()
    
    for (i in 1:length(sample_tb2_1$series)) {
      try(if (sample_tb2_1[i,]$series==min(sample_tb2_1$series)) {
        clu <- sample_tb2[1:sample_tb2_1[i,]$len,]
        clu <- clu[-1,]
        clu$clu <- i
        clu$clu.n <- nrow(clu)
        cluster[[i]] <- clu 
        print(clu)
      } else {
        clu <- sample_tb2[sample_tb2_1[i,]$len:sample_tb2_1[i+1,]$len,]
        clu <- clu[-1,]
        clu$clu <- i
        clu$clu.n <- nrow(clu)
        cluster[[i]] <- clu 
        print(clu)
      }, silent = TRUE)
    }
    
    data_cluster[[i]] <- bind_rows(cluster)
  }
  
  inact_data <- bind_rows(data_cluster) %>% filter(clu.n>1) %>%
    select(indiv, `5-min intervals`, value, threshold, clu, clu.n)
  
  bouts <- unique(inact_data %>%
                    select(indiv, clu)) %>%
    count(indiv, name = "bouts")
  
  inact_data2 <- left_join(inact_data, bouts) %>% 
    group_by(indiv) %>%
    mutate(ttl.time = n()*5)
  
  first_bout <- inact_data2 %>% group_by(indiv) %>% filter(clu==1) %>% 
    count(indiv, name = "first_bout") %>%
    transmute(first_bout_time = first_bout*5)
  
  inact_data3 <<- left_join(inact_data2, first_bout) %>% 
    group_by(indiv) %>%
    mutate(actual_time = ttl.time - first_bout_time,
           actual_bouts = bouts - 1)
  
  names(inact_data3) <<- c("indiv", "5-min intervals", "value < threshold", "threshold", "cluster_id", "cluster_size (>1)", "bouts", "total_time", "first_bout_time", "actual_time", "actual_bouts")
  
  write.csv(inact_data3, file = "inactivity_data.csv")
}