library(readr)
library(rlang)
library(tidyverse)
library(reshape2)
library(ggthemes)
library(rgdal) 
library(spatstat) 
library(nabor) 

##A function, for images imaged at 40x obj, 3.5x digital zoom, 1-micron z-step size,
##to compile puncta detected by SynapseJ v.1

##input is the label of the image as a string
##the function can be itterated over a string vector using sapply()

##the current output is a tibble in the R Environment, entitled all_puncta_[your image name]

##each row is a puncta observation, where:
##"Label"  (the name of the image with mouse number, cochlear region, and image number)
##"Area"  (the area in microns computed by ImageJ, NA for manual observations)       
##"Mean"   (the mean intensity. intensity of manual observations is just the pixel marked near the puncta)
##"desc"  (observations described as manual or pairing by SynapseJ algorithm)      
##"cluster_id"  (retained from the r function but indicates grouping of manual counts with detected puncta)
##"cluster_size" (retained from the r function but indicates how many components were connected in 3d before collapsing Z and binning by X and Y)
##"dist"  (the distance of a puncta in the pre channel to a puncta in the post channel [pixel units])     
##"channel" (the pre and post channels with puncta)  
##"ihc"  (the puncta corresponding to imaged inner hair cells numbered from left to right) 
##"X"  (pixel scale of puncta X location)
##"Y"  (pixel scale of puncta Y location)
##"Slice" (pixel scale Slice with puncta (Z) location)
##"nuc_X"  (pixel scale of IHC X location)
##"nuc_Y"  (pixel scale of IHC Y location)
##"nuc_Slice" (pixel scale Slice with IHC (Z) location)
##"norm_X" (pixel scale of puncta X location normalized to the corresponding inner hair cell "bin")
##"norm_Y"  (pixel scale of puncta Y location normalized to the inner hair cell "bin")
##"norm_Slice" (pixel scale Slice with puncta (Z) location normalized to the inner hair cell "bin")

##you need to know:
##1) microns per pixel for the images to be processed
##2) the file path to the detection tables from the SynapseJ macro
##3) the file path and name of the table containing nuclear locations for each inner hair cell
##Step 3 was done separately by an observer
##using ImageJ Cell Counter to mark nuclei centroids


puncta_per_hair_cell <- function(x){
  ##micron-to-pixel conversion factor
  mic_pix <- 0.0593047
  ##file path to SynapseJ detection tables
  path1 <- "demos/confocal image/"
  ##file path to the SynapseJ excel subfolder
  path2 <- paste(path1, "excel/", sep = "")
  ##file path to tables with manual counts, such as inner hair cell nuclei locations per image
  path3 <- "demos/confocal image/"
  
  ##import data
  post <- read_delim(paste(path1, expr(!!x), "Post.txt", sep = ""), 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE) %>% dplyr::select(-`...1`)
  
  post_results <- read_delim(paste(path2, expr(!!x), "PostResults.txt", sep = ""), 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE) %>% dplyr::select(-`...1`)
  
  ##manual cell counter counts
  #  manual_post <- read_csv(paste(path3, expr(!!x), "_post.csv", sep = "")) %>%
  #    mutate(Label = expr(!!x),
  #           desc = "manual",
  #           Mean = Value,
  #           X = `X(micron)`/mic_pix,
  #           Y = `Y(micron)`/mic_pix,
  #           Slice = Slice/mic_pix) %>% 
  #    dplyr::select(-Type, -Value, -`C-pos`, -`Z-pos`, -`T-pos`, -`X(micron)`, -`Y(micron)`, -`Z(micron)`)
  
  pre <- read_delim(paste(path1, expr(!!x), "Pre.txt", sep = ""), 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE) %>% dplyr::select(-`...1`)
  
  pre_results <- read_delim(paste(path2, expr(!!x), "PreResults.txt", sep = ""), 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE) %>% dplyr::select(-`...1`)
  
  ##manual cell counter counts
  #  manual_pre <- read_csv(paste(path3, expr(!!x), "_pre.csv", sep = "")) %>%
  #    mutate(Label = expr(!!x),
  #           desc = "manual",
  #           Mean = Value,
  #           X = `X(micron)`/mic_pix,
  #           Y = `Y(micron)`/mic_pix,
  #           Slice = Slice/mic_pix) %>% 
  #    dplyr::select(-Type, -Value, -`C-pos`, -`Z-pos`, -`T-pos`, -`X(micron)`, -`Y(micron)`, -`Z(micron)`)
  
  nuc <- read_csv(paste(path3, expr(!!x), "_nuc.csv", sep = "")) %>%
    mutate(Label = expr(!!x),
           desc = "manual",
           Mean = Value,
           X = `X(micron)`/mic_pix,
           Y = `Y(micron)`/mic_pix,
           Slice = Slice/mic_pix) %>%
    dplyr::select(-Type, -Value, -`C-pos`, -`Z-pos`, -`T-pos`, -`X(micron)`, -`Y(micron)`, -`Z(micron)`)
  
  
  ##prep and clean post data
  synaptic_post <- post_results %>%
    mutate(desc = "pair")
  
  ##non-synaptic post
  post_puncta <- anti_join(post, post_results,
                           by = c("Area", "Mean", "Min", "Max", "X", "Y", "XM", "YM", "Perim.", "Feret", "IntDen", "RawIntDen", "Slice", "FeretX", "FeretY", "FeretAngle", "MinFeret")) %>% 
    mutate(desc = "orphan")
  
  all_post <- bind_rows(synaptic_post, post_puncta)
  ##convert X, Y, and Slice to pixels, 1px = mic_pix microns (for 40x obj, 3.5x digital zoom)
  all_post <- all_post %>%
    mutate(X = X/mic_pix,
           Y = Y/mic_pix,
           Slice = Slice/mic_pix)
  
  #  all_post <- bind_rows(all_post, manual_post) 
  
  all_post <- all_post %>%
    mutate(id = 1:length(Label))
  
  
  ##point cloud clustering with spatstat
  post_pc <- all_post %>% dplyr::select(X, Y, Slice, id) %>%
    mutate(X = as.numeric(X),
           Y = as.numeric(Y),
           Slice = as.numeric(Slice))
  
  class(post_pc) <- 'data.frame'
  
  bb_post <- box3(range(post_pc[,1]), range(post_pc[,2]), range(post_pc[,3]))
  
  pc_post <- pp3(post_pc[,1], post_pc[,2], post_pc[,3], bb_post)
  
  pc_post_labelled <- connected.pp3(pc_post, R = 1.5/mic_pix)
  
  pc_post_df <- data.frame(cluster_id = marks(pc_post_labelled), id = post_pc$id)
  
  all_post <- left_join(all_post, pc_post_df)
  
  
  ##cluster segmentation
  #  cluster_post <- all_post %>% #filter(desc!="manual") %>%
  #    dplyr::select(X, Y, Mean, id, cluster_id)
  #  
  #  cluster_post_id <- unique(cluster_post$cluster_id)
  
  #  cluster_post_list <- list()
  
  #  for (i in 1:length(cluster_post_id)){
  #    clu <- cluster_post %>% filter(cluster_id==cluster_post_id[i])
  #    if (nrow(clu)>1) {
  #      xy <- SpatialPointsDataFrame(matrix(c(clu$X, clu$Y), ncol=2), 
  #                                   data.frame(inter_id=seq(1:length(clu$X))))
  
  #      cxy <- hclust(dist(data.frame(rownames=rownames(xy@data),
  #                                    x=coordinates(xy)[,1],
  #                                    y=coordinates(xy)[,2])), 
  #                    method = "complete")
  
  #      #0.5 micron threshold (8.43 px radius)
  #      cxy_rad <- cutree(cxy, h = 8.43)
  
  #join inter cluster clustering results with connected component point ids
  #      idc <- data.frame(xy@data, inter_clust=cxy_rad)
  #      as_tibble(idc)
  
  ##join with the connected component data and return the
  ##brightest point per inter cluster
  #      clu2 <- bind_cols(clu, idc) %>% group_by(inter_clust) %>%
  #        mutate(max = max(Mean)) %>%
  #        filter(Mean==max) %>% ungroup %>% dplyr::select(cluster_id, id)
  
  #      cluster_list[[i]] <- clu2
  #    } else {
  #      clu2 <- clu %>%
  #        mutate(max = Mean) %>% 
  #        dplyr::select(cluster_id, id)
  
  #      cluster_post_list[[i]] <- clu2
  #    }
  #  }
  
  
  #  cluster_post2 <- dplyr::bind_rows(cluster_post_list)
  
  #  cluster_post2 <- left_join(cluster_post2, all_post)
  
  #  all_post4 <- cluster_post2 %>% ungroup %>% #bind_rows(max_post, man_post) %>% ungroup %>%
  #    dplyr::select(Label, X, Y, Area, Mean, Slice, desc, id, cluster_id) %>%
  #    group_by(cluster_id) %>%
  #    mutate(cluster_size = n())
  
  
  
  ##pre synaptic puncta channel data cleaning and prep
  synaptic_pre <- pre_results %>%
    mutate(desc = "pair")
  
  ##non-synaptic pre
  pre_puncta <- anti_join(pre, pre_results,
                          by = c("Area", "Mean", "Min", "Max", "X", "Y", "XM", "YM", "Perim.", "Feret", "IntDen", "RawIntDen", "Slice", "FeretX", "FeretY", "FeretAngle", "MinFeret")) %>% 
    mutate(desc = "orphan")
  
  all_pre <- bind_rows(synaptic_pre, pre_puncta)
  ##convert X, Y, and Slice to pixels, 1px = mic_pix microns (for 40x obj, 3.5x digital zoom)
  all_pre <- all_pre %>%
    mutate(X = X/mic_pix,
           Y = Y/mic_pix,
           Slice = Slice/mic_pix)
  
  #  all_pre <- bind_rows(all_pre, manual_pre)
  
  all_pre <- all_pre %>%
    mutate(id = 1:length(Label))
  
  
  ##point cloud clustering with spatstat (connected components)
  pre_pc <- all_pre %>% dplyr::select(X, Y, Slice, id) %>%
    mutate(X = as.numeric(X),
           Y = as.numeric(Y),
           Slice = as.numeric(Slice))
  
  class(pre_pc) <- 'data.frame'
  
  bb_pre <- box3(range(pre_pc[,1]), range(pre_pc[,2]), range(pre_pc[,3]))
  
  pc_pre <- pp3(pre_pc[,1], pre_pc[,2], pre_pc[,3], bb_pre)
  
  pc_pre_labelled <- connected.pp3(pc_pre, R = 1.5/mic_pix)
  
  pc_pre_df <- data.frame(cluster_id = marks(pc_pre_labelled), id = pre_pc$id)
  
  all_pre <- left_join(all_pre, pc_pre_df)
  
  
  ##connected component segmentation with sp and rgdal packages
  ##components within a radius of each other are considered the same object and "connected"
  cluster_pre <- all_pre %>% #filter(desc!="manual") %>%
    dplyr::select(X, Y, Mean, id, cluster_id)
  
  cluster_pre_id <- unique(cluster_pre$cluster_id)
  
  cluster_pre_list <- list()
  
  
  for (i in 1:length(cluster_pre_id)){
    clu <- cluster_pre %>% filter(cluster_id==cluster_pre_id[i])
    if (nrow(clu)>1) {
      xy <- SpatialPointsDataFrame(matrix(c(clu$X, clu$Y), ncol=2), 
                                   data.frame(inter_id=seq(1:length(clu$X))))
      
      cxy <- hclust(dist(data.frame(rownames=rownames(xy@data),
                                    x=coordinates(xy)[,1],
                                    y=coordinates(xy)[,2])), 
                    method = "complete")
      
      #0.5 micron threshold (8.43 px radius)
      cxy_rad <- cutree(cxy, h = 0.5/mic_pix)
      
      #join inter-cluster clustering results with connected component point ids
      idc <- data.frame(xy@data, inter_clust=cxy_rad)
      as_tibble(idc)
      
      ##join with the connected component data and return the
      ##brightest point per inter cluster
      clu2 <- bind_cols(clu, idc) %>% group_by(inter_clust) %>%
        mutate(max = max(Mean)) %>%
        filter(Mean==max) %>% ungroup %>% dplyr::select(cluster_id, id)
      
      cluster_pre_list[[i]] <- clu2
    } else {
      clu2 <- clu %>%
        mutate(max = Mean) %>% 
        dplyr::select(cluster_id, id)
      
      cluster_pre_list[[i]] <- clu2
    }
  }
  
  
  cluster_pre2 <- dplyr::bind_rows(cluster_pre_list)
  
  cluster_pre2 <- left_join(cluster_pre2, all_pre)
  
  all_pre4 <- cluster_pre2 %>% ungroup %>% #bind_rows(max_pre, man_pre) %>% ungroup %>%
    dplyr::select(Label, X, Y, Area, Mean, Slice, desc, id, cluster_id) %>%
    group_by(cluster_id) %>%
    mutate(cluster_size = n())
  
  
  
  ##nearest neighbor
  post_nbr <- all_post #%>% filter(desc!="manual") 
  post_nbr_xyz <- post_nbr %>% ungroup %>% dplyr::select(X, Y, Slice)
  
  pre_nbr <- all_pre4 #%>% filter(desc!="manual")
  pre_nbr_xyz <- pre_nbr %>% ungroup %>% dplyr::select(X, Y, Slice)
  
  ##post puncta within 2 microns (equivalent in pixels) of pre puncta
  post_nn <- knn(data = pre_nbr_xyz, query = post_nbr_xyz, k = 1, radius = 2/mic_pix)
  
  post_nbr$dist <- post_nn$nn.dists
  
  all_post <- left_join(all_post, post_nbr) %>%
    mutate(channel = "post")
  
  ##pre puncta within 2 microns of post puncta
  pre_nn <- knn(data = post_nbr_xyz, query = pre_nbr_xyz, k = 1, radius = 2/mic_pix)
  
  pre_nbr$dist <- pre_nn$nn.dists
  
  all_pre4 <- left_join(all_pre4, pre_nbr) %>% 
    mutate(channel = "pre")
  
  all_puncta <- bind_rows(all_pre4, all_post) %>% ungroup %>%
    dplyr::select(Label, X, Y, Area, Mean, Slice, desc, dist, channel, id, cluster_id)
  
  
  
  
  ##now filter and bin the results by IHC nucleus in the image
  
  nuc_coords <- arrange(nuc, X)
  bw <- tibble(n = 2:nrow(nuc_coords)-1)
  for (i in 1:length(nuc_coords$X)){
    bw[i,] = (nuc_coords[i,2]+ nuc_coords[i+1,2])/2  ##indexed column should be X
  }
  bw <- drop_na(bw)
  
  avg_bw <- bw %>%
    mutate(nuc = nuc_coords[1:nrow(bw),2],
           dist_bw = n-nuc$X,
           med_bw = median(dist_bw))
  
  med_val <- unique(avg_bw$med_bw)
  
  
  ##filtering out puncta in the nuclear region
  nuc_xyz <- nuc %>% dplyr::select(X, Y, Slice)
  all_puncta_xyz <- all_puncta %>% ungroup %>% dplyr::select(X, Y, Slice)
  
  ##diameter of DAPI staining is ~7.8 microns
  ##(7.8/2)/mic_pix = 65.76207
  all_puncta_nn <- knn(data = nuc_xyz, query = all_puncta_xyz, k = 1) 
  
  all_puncta$nuc_dist <- all_puncta_nn$nn.dists
  
  ##the working example for this code (770_22.6_1) did not have any of my manual counts within close 
  ##proximity to the nucleus. Puncta detected near the nucleus will be filtered out in this next step.
  all_puncta2 <- all_puncta %>% filter(nuc_dist>(3.9/mic_pix)) %>%
    dplyr::select(-nuc_dist)
  
  all_puncta2.1 <- all_puncta %>% filter(nuc_dist<(3.9/mic_pix)) %>% dplyr::select(channel, cluster_id)
  
  all_puncta2 <-anti_join(all_puncta, all_puncta2.1) 
  
  
  if(length(nuc_coords$X)==1) {
    ihc1 <- all_puncta2 %>% mutate(ihc="1")
    ihcs <- bind_rows(ihc1)
  } else if (length(nuc_coords$X)==2) {
    if(nuc_coords$X[1]-med_val>=0) {
      ihc1 <- all_puncta2 %>% filter(X>=nuc_coords$X[1]-med_val&X<bw$n[1]) %>% mutate(ihc="1")
    } else {
      ihc1 <- all_puncta2 %>% filter(X<bw$n[1]) %>% mutate(ihc="1")
    }
    if(nuc_coords$X[2]-med_val<1024) {
      ihc2 <- all_puncta2 %>% filter(X>=bw$n[1]&X<nuc_coords$X[2]+med_val) %>% mutate(ihc = "2")
    } else {
      ihc2 <- all_puncta2 %>% filter(X>=bw$n[1]) %>% mutate(ihc = "2")
    } 
    ihcs <- bind_rows(ihc1, ihc2)
  } else if(length(nuc_coords$X)==3) {
    if(nuc_coords$X[1]-med_val>=0) {
      ihc1 <- all_puncta2 %>% filter(X>=nuc_coords$X[1]-med_val&X<bw$n[1]) %>% mutate(ihc="1")
    } else {
      ihc1 <- all_puncta2 %>% filter(X<bw$n[1]) %>% mutate(ihc="1")
    }
    ihc2 <- all_puncta2 %>% filter(X>=bw$n[1]&X<bw$n[2]) %>% mutate(ihc = "2")
    if(nuc_coords$X[3]-med_val<1024) {
      ihc3 <- all_puncta2 %>% filter(X>=bw$n[2]&X<nuc_coords$X[3]+med_val) %>% mutate(ihc = "3")
    } else {
      ihc3 <- all_puncta2 %>% filter(X>=bw$n[2]) %>% mutate(ihc = "3")
    }  
    ihcs <- bind_rows(ihc1, ihc2, ihc3)
  } else if(length(nuc_coords$X)==4) {
    if(nuc_coords$X[1]-med_val>=0) {
      ihc1 <- all_puncta2 %>% filter(X>=nuc_coords$X[1]-med_val&X<bw$n[1]) %>% mutate(ihc="1")
    } else {
      ihc1 <- all_puncta2 %>% filter(X<bw$n[1]) %>% mutate(ihc="1")
    }
    ihc2 <- all_puncta2 %>% filter(X>=bw$n[1]&X<bw$n[2]) %>% mutate(ihc = "2")
    ihc3 <- all_puncta2 %>% filter(X>=bw$n[2]&X<bw$n[3]) %>% mutate(ihc = "3")
    if(nuc_coords$X[4]-med_val<1024) {
      ihc4 <- all_puncta2 %>% filter(X>=bw$n[3]&X<nuc_coords$X[4]+med_val) %>% mutate(ihc = "4")
    } else {
      ihc4 <- all_puncta2 %>% filter(X>=bw$n[3]) %>% mutate(ihc = "4")
    }  
    ihcs <- bind_rows(ihc1, ihc2, ihc3, ihc4)
  } else if(length(nuc_coords$X)==5) {
    if(nuc_coords$X[1]-med_val>=0) {
      ihc1 <- all_puncta2 %>% filter(X>=nuc_coords$X[1]-med_val&X<bw$n[1]) %>% mutate(ihc="1")
    } else {
      ihc1 <- all_puncta2 %>% filter(X<bw$n[1]) %>% mutate(ihc="1")
    }
    ihc2 <- all_puncta2 %>% filter(X>=bw$n[1]&X<bw$n[2]) %>% mutate(ihc = "2")
    ihc3 <- all_puncta2 %>% filter(X>=bw$n[2]&X<bw$n[3]) %>% mutate(ihc = "3")
    ihc4 <- all_puncta2 %>% filter(X>=bw$n[3]&X<bw$n[4]) %>% mutate(ihc = "4")
    if(nuc_coords$X[5]-med_val<1024) {
      ihc5 <- all_puncta2 %>% filter(X>=bw$n[4]&X<nuc_coords$X[5]+med_val) %>% mutate(ihc = "5")
    } else {
      ihc5 <- all_puncta2 %>% filter(X>=bw$n[4]) %>% mutate(ihc = "5")
    }  
    ihcs <- bind_rows(ihc1, ihc2, ihc3, ihc4, ihc5)
  } else if(length(nuc_coords$X)==6) {
    if(nuc_coords$X[1]-med_val>=0) {
      ihc1 <- all_puncta2 %>% filter(X>=nuc_coords$X[1]-med_val&X<bw$n[1]) %>% mutate(ihc="1")
    } else {
      ihc1 <- all_puncta2 %>% filter(X<bw$n[1]) %>% mutate(ihc="1")
    }
    ihc2 <- all_puncta2 %>% filter(X>=bw$n[1]&X<bw$n[2]) %>% mutate(ihc = "2")
    ihc3 <- all_puncta2 %>% filter(X>=bw$n[2]&X<bw$n[3]) %>% mutate(ihc = "3")
    ihc4 <- all_puncta2 %>% filter(X>=bw$n[3]&X<bw$n[4]) %>% mutate(ihc = "4")
    ihc5 <- all_puncta2 %>% filter(X>=bw$n[4]&X<bw$n[5]) %>% mutate(ihc = "5")
    if(nuc_coords$X[6]-med_val<1024) {
      ihc6 <- all_puncta2 %>% filter(X>=bw$n[5]&X<nuc_coords$X[6]+med_val) %>% mutate(ihc = "6")
    } else {
      ihc6 <- all_puncta2 %>% filter(X>=bw$n[5]) %>% mutate(ihc = "6")
    }
    ihcs <- bind_rows(ihc1, ihc2, ihc3, ihc4, ihc5, ihc6)
  } else if(length(nuc_coords$X)==7) {
    if(nuc_coords$X[1]-med_val>=0) {
      ihc1 <- all_puncta2 %>% filter(X>=nuc_coords$X[1]-med_val&X<bw$n[1]) %>% mutate(ihc="1")
    } else {
      ihc1 <- all_puncta2 %>% filter(X<bw$n[1]) %>% mutate(ihc="1")
    }
    ihc2 <- all_puncta2 %>% filter(X>=bw$n[1]&X<bw$n[2]) %>% mutate(ihc = "2")
    ihc3 <- all_puncta2 %>% filter(X>=bw$n[2]&X<bw$n[3]) %>% mutate(ihc = "3")
    ihc4 <- all_puncta2 %>% filter(X>=bw$n[3]&X<bw$n[4]) %>% mutate(ihc = "4")
    ihc5 <- all_puncta2 %>% filter(X>=bw$n[4]&X<bw$n[5]) %>% mutate(ihc = "5")
    ihc6 <- all_puncta2 %>% filter(X>=bw$n[5]&X<bw$n[6]) %>% mutate(ihc = "6")
    if(nuc_coords$X[7]-med_val<1024) {
      ihc7 <- all_puncta2 %>% filter(X>=bw$n[6]&X<nuc_coords$X[7]+med_val) %>% mutate(ihc = "7")
    } else {
      ihc7 <- all_puncta2 %>% filter(X>=bw$n[6]) %>% mutate(ihc = "7")
    }
    ihcs <- bind_rows(ihc1, ihc2, ihc3, ihc4, ihc5, ihc6, ihc7)
  } else if(length(nuc_coords$X)==8) {
    if(nuc_coords$X[1]-med_val>=0) {
      ihc1 <- all_puncta2 %>% filter(X>=nuc_coords$X[1]-med_val&X<bw$n[1]) %>% mutate(ihc="1")
    } else {
      ihc1 <- all_puncta2 %>% filter(X<bw$n[1]) %>% mutate(ihc="1")
    }
    ihc2 <- all_puncta2 %>% filter(X>=bw$n[1]&X<bw$n[2]) %>% mutate(ihc = "2")
    ihc3 <- all_puncta2 %>% filter(X>=bw$n[2]&X<bw$n[3]) %>% mutate(ihc = "3")
    ihc4 <- all_puncta2 %>% filter(X>=bw$n[3]&X<bw$n[4]) %>% mutate(ihc = "4")
    ihc5 <- all_puncta2 %>% filter(X>=bw$n[4]&X<bw$n[5]) %>% mutate(ihc = "5")
    ihc6 <- all_puncta2 %>% filter(X>=bw$n[5]&X<bw$n[6]) %>% mutate(ihc = "6")
    ihc7 <- all_puncta2 %>% filter(X>=bw$n[6]&X<bw$n[7]) %>% mutate(ihc = "7")
    if(nuc_coords$X[8]-med_val<1024) {
      ihc8 <- all_puncta2 %>% filter(X>=bw$n[7]&X<nuc_coords$X[8]+med_val) %>% mutate(ihc = "8")
    } else {
      ihc8 <- all_puncta2 %>% filter(X>=bw$n[7]) %>% mutate(ihc = "8")
    }
    ihcs <- bind_rows(ihc1, ihc2, ihc3, ihc4, ihc5, ihc6, ihc7, ihc8)
  } else if(length(nuc_coords$X)==9) {
    if(nuc_coords$X[1]-med_val>=0) {
      ihc1 <- all_puncta2 %>% filter(X>=nuc_coords$X[1]-med_val&X<bw$n[1]) %>% mutate(ihc="1")
    } else {
      ihc1 <- all_puncta2 %>% filter(X<bw$n[1]) %>% mutate(ihc="1")
    }
    ihc2 <- all_puncta2 %>% filter(X>=bw$n[1]&X<bw$n[2]) %>% mutate(ihc = "2")
    ihc3 <- all_puncta2 %>% filter(X>=bw$n[2]&X<bw$n[3]) %>% mutate(ihc = "3")
    ihc4 <- all_puncta2 %>% filter(X>=bw$n[3]&X<bw$n[4]) %>% mutate(ihc = "4")
    ihc5 <- all_puncta2 %>% filter(X>=bw$n[4]&X<bw$n[5]) %>% mutate(ihc = "5")
    ihc6 <- all_puncta2 %>% filter(X>=bw$n[5]&X<bw$n[6]) %>% mutate(ihc = "6")
    ihc7 <- all_puncta2 %>% filter(X>=bw$n[6]&X<bw$n[7]) %>% mutate(ihc = "7")
    ihc8 <- all_puncta2 %>% filter(X>=bw$n[7]&X<bw$n[8]) %>% mutate(ihc = "8")
    if(nuc_coords$X[9]-med_val<1024) {
      ihc9 <- all_puncta2 %>% filter(X>=bw$n[8]&X<nuc_coords$X[9]+med_val) %>% mutate(ihc = "9")
    } else {
      ihc9 <- all_puncta2 %>% filter(X>=bw$n[8]) %>% mutate(ihc = "9")
    }
    ihcs <- bind_rows(ihc1, ihc2, ihc3, ihc4, ihc5, ihc6, ihc7, ihc8, ihc9)
  } else {
    paste("too many hair cells")
  }
  
  ihcs <- ihcs %>% 
    mutate(Label = expr(!!x),
           ihc = as.numeric(ihc))
  
  nuc_xyz <- nuc %>% arrange(X) %>%
    mutate(nuc_X = X,
           nuc_Y = Y,
           nuc_Slice = Slice,
           ihc = 1:length(X)) %>%
    dplyr::select(nuc_X, nuc_Y, nuc_Slice, ihc)
  
  ihcs2 <- left_join(ihcs, nuc_xyz) %>%
    group_by(ihc) %>%
    mutate(norm_X = X - nuc_X,
           norm_Y = Y - nuc_Y,
           norm_Slice = Slice - nuc_Slice) %>%
    dplyr::select(-id, -nuc_dist)
  
  assign(paste('all_puncta_', expr(!!x), sep = ""), ihcs2, envir = global_env())
}

#puncta_per_hair_cell("1000_22.6_2")



                         



