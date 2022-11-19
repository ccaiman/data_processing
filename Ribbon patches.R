library(tidyverse)
library(ggplot2)
library(imager)
library(reshape2)
library(patchwork)
library(magick)
library(raster)
library(figpatch)
library(readr)
library(rlang)


##This function takes as input a tibbles generated from the puncta_per_hair_cell() function
##(SynapseJ data compilation script)

##The output is multiple tibbles in the R Environment corresponding to definitions of puncta:
##synapse_[image label]: presynaptic puncta paired with a postsynaptic puncta
##syanpse_drop[image label]: a small presynaptic puncta (size-filtered from synapse dataset)
##that is paired with a postsynaptic puncta
##orphan_[image label]: presynaptic puncta that is not paired with a postsynaptic puncta

##the regions of interest from each of these tibbles is used to extract a patch from the 
##image and print them to either the "pre" or "post" folder (within the "patches" folder)

##from there, the user is encouraged to determine the ribbon synapses and orphaned ribbons
##by matching the patches from the pre and post channels
##the result is guided, manual counts

##in this script, you need to:
##1)set the path to the tiffs which were processed by the SynapseJ algorithm
##so this function can extract patches from the image
##2)set the micron-to-pixel conversion factor
##3)ensure there is a "patches" folder in your working directory, 
##with a "pre" and "post" subfolder 

SynapseJ_patches <- function(x){
  
  ##micron-to-pixel conversion factor
  mic_pix <- 0.0593047
  ##path to the image
  path = "demos/confocal image/"
  
  ##working with CtBP2 puncta detected by SynapseJ for all figures
  data <- as_tibble(expr(!!x))
  
  detected_3d_xy <- data %>% filter(desc!="manual"&channel=="pre")
  #filter(norm_Y<500) ##my manual counts did not extend so far from the nucleus
  ##though perhaps a more precise filtering will be to group by label and ihc and take everything below the lowest synapse
  
  detected_lob_key <- left_join(detected_3d_xy, detected_3d_xy %>% filter(desc!="manual"&channel=="pre") %>%
                                  #filter(norm_Y<500) ##my manual counts did not extend so far from the nucleus
                                  ##though perhaps a more precise filtering will be to group by label and ihc and take everything below the lowest synapse
                                  group_by(Label, ihc) %>%
                                  filter(dist<40) %>%
                                  #filter(Area>0.35) %>% #obtaining the limit to "search" for synapses based on a typical minimum size (area) exclusion
                                  mutate(lowest_ob = max(norm_Y)) %>%
                                  ungroup)
  
  detected_lob_key <- unique(drop_na(detected_lob_key %>% dplyr::select(Label, ihc, lowest_ob)))
  
  detected_3d_xy <- left_join(detected_3d_xy, detected_lob_key)
  
  detected_3d_xy <- detected_3d_xy %>% 
    group_by(Label, ihc) %>% 
    filter(norm_Y<lowest_ob) %>%
    ungroup %>%
    dplyr::select(!lowest_ob)
  
  labs <- unique(detected_3d_xy %>% ungroup %>% dplyr::select(Label))
  
  labs <- labs$Label
  
  ##size filtering
  detected_3d_xy_szf <- detected_3d_xy %>% filter(Area>0.35)
  
  ##the puncta which were dropped
  detected_3d_xy_filtered <- anti_join(detected_3d_xy, detected_3d_xy_szf)
  
  ##any synapses (within 2 microns of GluA2) that were dropped? We would like to see those.
  synapse_dropped <- detected_3d_xy_filtered %>% filter(dist<40)
  
  
  for (i in 1:length(labs)) {
    ##now, iteration over each Label
    detected <- detected_3d_xy_szf %>% ungroup %>% filter(Label==labs[i])
    
    ##detected synapses with area>0.35 square microns
    synapse <- detected %>% ungroup %>% filter(dist<40) %>% arrange(X,Y)
    synapse_id <- 1:nrow(synapse)
    synapse$id <- synapse_id
    
    assign(paste('synapse_', labs[i], sep = ""), synapse, envir = global_env())
    
    
    ##load the tiff image
    img <- load.image(paste(path, labs[i], "/", labs[i], ".tif", sep = ""))
    
    z_dim <- dim(img)[3] ##return the z stack dimension
    which(1:z_dim%%2!=1) ##used below: split the z dimension where odd vals and even vals each correspond to a channel
    
    ##odd values comprise the pre channel
    imf1 <- frame(img, which(1:z_dim%%2==1))
    
    ##even values comprise the post channel
    imf2 <- frame(img, which(1:z_dim%%2!=1))
    
    
    ##extract patches of synapses
    synapse_extract_pre <- extract_patches3D(imf1, cx = synapse$X, cy = synapse$Y,
                                             cz = synapse$Slice*mic_pix, wx = 1.5/mic_pix, wy = 1.5/mic_pix, wz = 2) ##convert the slice back to integer
    
    synapse_extract_post <- extract_patches3D(imf2, cx = synapse$X, cy = synapse$Y,
                                              cz = synapse$Slice*mic_pix, wx = 1.5/mic_pix, wy = 1.5/mic_pix, wz = 2)
    
    ##convert image list to a data frame
    df_pre <- as.data.frame(synapse_extract_pre)
    
    ##convert to a tibble for tidy transformations
    df_pre_tb <- as_tibble(df_pre)
    
    ##cast the z stack around puncta of interest for MIP calculation
    df_pre_tb2 <- dcast(df_pre_tb, im+x+y~z, value.var = "value") %>%
      mutate(max = pmax(`1`, `2`, `3`)) %>%
      dplyr::select(-`1`, -`2`, -`3`)
    
    ##convert im from the image list to a number
    df_pre_tb2$im2 <- as.numeric(df_pre_tb2$im)
    
    ##make an id column to match the id of the puncta population(synapse, dropped synapse, or orphan)
    df_pre_tb2 <- df_pre_tb2 %>% mutate(id = im2)
    
    ihc_id <- synapse %>% ungroup %>% dplyr::select(id, ihc)
    
    ##join the ihc info to the image table
    df_pre_tb2 <- left_join(df_pre_tb2, ihc_id)
    
    ##make text labels of ihc number to go with each id'd image
    dat_text <- data.frame(
      label = as.character(ihc_id$ihc),
      im2   = as.numeric(ihc_id$id)
    )
    
    ##the facet for the puncta subtype (synapse, dropped synapse, or orphan)
    n_facet <- round(sqrt(max(unique(df_pre_tb2$im2))))
    
    ##make the ggplot and save
    ggplot(data = df_pre_tb2, mapping = aes(x, y)) +
      geom_raster(aes(fill=max)) +
      scale_y_continuous(trans=scales::reverse_trans()) +
      geom_text(data = dat_text, mapping = aes(x = 20, y = 23, label = label), 
                hjust   = 0, vjust   = 0, size = 3, color = "white") +
      facet_wrap(~im2, ncol = n_facet+2) +
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
        strip.text.y = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
        plot.title = element_text(size = 14)
      ) +
      scale_fill_gradient(low="black", high="magenta") +
      coord_fixed(ratio = 1) +
      ggtitle(paste(labs[i], "_pre_synapse (0.35 < x <= 1.4 square microns)", sep = ""))
    
    ggsave(paste(labs[i], "_pre_synapse.tiff", sep = ""), path = "patches/pre/", width = (n_facet+2)/1.7, height = ((n_facet+2)/1.7)-1, units = "in")
    
    ##now working with the post channel
    ##there are less steps because geom assets are reused from above
    df_post <- as.data.frame(synapse_extract_post)
    
    df_post_tb <- as_tibble(df_post)
    
    df_post_tb2 <- dcast(df_post_tb, im+x+y~z, value.var = "value") %>%
      mutate(max = pmax(`1`, `2`, `3`))
    
    df_post_tb2$im2 <- as.numeric(df_post_tb2$im)
    
    ggplot(data = df_post_tb2, mapping = aes(x, y)) +
      geom_raster(aes(fill=max)) +
      scale_y_continuous(trans=scales::reverse_trans()) +
      geom_text(data = dat_text, mapping = aes(x = 20, y = 23, label = label), 
                hjust   = 0, vjust   = 0, size = 3, color = "white") +
      facet_wrap(~im2, ncol = n_facet+2) +
      theme(
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        strip.text.x = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
        strip.text.y = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
        plot.title = element_text(size = 14)
      ) +
      scale_fill_gradient(low="black", high="green") +
      coord_fixed(ratio = 1) +
      ggtitle(paste(labs[i], "_post_synapse (within 2 microns of pre)", sep = ""))
    
    ggsave(paste(labs[i], "_post_synapse.tiff", sep = ""), path = "patches/post/", width = (n_facet+2)/1.7, height = ((n_facet+2)/1.7)-1, units = "in")
    
    
    
    
    if (nrow(detected)>nrow(synapse)){
      ##limit orphan CtBP2 punctas to the size-filtered data set
      orphan <- anti_join(detected, synapse) %>% arrange(X,Y)
      orphan_id <- 1:nrow(orphan)
      orphan$id <- orphan_id
      
      assign(paste('orphan_', labs[i], sep = ""), orphan, envir = global_env())
      
      ##extract patches of orphans
      orphan_extract_pre <- extract_patches3D(imf1, cx = orphan$X, cy = orphan$Y,
                                              cz = orphan$Slice*mic_pix, wx = 1.5/mic_pix, wy = 1.5/mic_pix, wz = 2) ##convert the slice back to numeric
      
      orphan_extract_post <- extract_patches3D(imf2, cx = orphan$X, cy = orphan$Y,
                                               cz = orphan$Slice*mic_pix, wx = 1.5/mic_pix, wy = 1.5/mic_pix, wz = 2)
      
      ##convert image list to a data frame
      df_orph_pre <- as.data.frame(orphan_extract_pre)
      
      ##convert to a tibble for tidy transformations
      df_orph_pre_tb <- as_tibble(df_orph_pre)
      
      ##cast the z stack around puncta of interest for MIP calculation
      df_orph_pre_tb2 <- dcast(df_orph_pre_tb, im+x+y~z, value.var = "value") %>%
        mutate(max = pmax(`1`, `2`, `3`)) %>%
        dplyr::select(-`1`, -`2`, -`3`)
      
      ##convert im from the image list to a number
      df_orph_pre_tb2$im2 <- as.numeric(df_orph_pre_tb2$im)
      
      ##make an id column to match the id of the puncta population(synapse, dropped synapse, or orphan)
      df_orph_pre_tb2 <- df_orph_pre_tb2 %>% mutate(id = im2)
      
      orph_ihc_id <- orphan %>% ungroup %>% dplyr::select(id, ihc)
      
      ##join the ihc info to the image table
      df_orph_pre_tb2 <- left_join(df_orph_pre_tb2, orph_ihc_id)
      
      ##make text labels of ihc number to go with each id'd image
      orph_dat_text <- data.frame(
        label = as.character(orph_ihc_id$ihc),
        im2   = as.numeric(orph_ihc_id$id)
      )
      
      ##the facet for the puncta subtype (synapse, dropped synapse, or orphan)
      n_facet_orph <- round(sqrt(max(unique(df_orph_pre_tb2$im2))))
      
      ##make the ggplot and save
      ggplot(data = df_orph_pre_tb2, mapping = aes(x, y)) +
        geom_raster(aes(fill=max)) +
        scale_y_continuous(trans=scales::reverse_trans()) +
        geom_text(data = orph_dat_text, mapping = aes(x = 20, y = 23, label = label), 
                  hjust   = 0, vjust   = 0, size = 3, color = "white") +
        facet_wrap(~im2, ncol = n_facet_orph+2) +
        theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          strip.text.x = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
          strip.text.y = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
          plot.title = element_text(size = 14)
        ) +
        scale_fill_gradient(low="black", high="magenta") +
        coord_fixed(ratio = 1) +
        ggtitle(paste(labs[i], "_pre_orphan (0.35 <= x <= 1.4 square microns)", sep = ""))
      
      ggsave(paste(labs[i], "_pre_orphan.tiff", sep = ""), path = "patches/pre/", width = (n_facet_orph+2)/1.7, height = ((n_facet_orph+2)/1.7)-1, units = "in")
      
      ##now working with the post channel
      ##there are less steps because geom assessts are reused from above
      df_orph_post <- as.data.frame(orphan_extract_post)
      
      df_orph_post_tb <- as_tibble(df_orph_post)
      
      df_orph_post_tb2 <- dcast(df_orph_post_tb, im+x+y~z, value.var = "value") %>%
        mutate(max = pmax(`1`, `2`, `3`))
      
      df_orph_post_tb2$im2 <- as.numeric(df_orph_post_tb2$im)
      
      ggplot(data = df_orph_post_tb2, mapping = aes(x, y)) +
        geom_raster(aes(fill=max)) +
        scale_y_continuous(trans=scales::reverse_trans()) +
        geom_text(data = orph_dat_text, mapping = aes(x = 20, y = 23, label = label), 
                  hjust   = 0, vjust   = 0, size = 3, color = "white") +
        facet_wrap(~im2, ncol = n_facet_orph+2) +
        theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          strip.text.x = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
          strip.text.y = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
          plot.title = element_text(size = 14)
        ) +
        scale_fill_gradient(low="black", high="green") +
        coord_fixed(ratio = 1) +
        ggtitle(paste(labs[i], "_post_orphan (within 2 microns of pre)", sep = ""))
      
      ggsave(paste(labs[i], "_post_orphan.tiff", sep = ""), path = "patches/post/", width = (n_facet_orph+2)/1.7, height = ((n_facet_orph+2)/1.7)-1, units = "in")
    } else {
      print("there were no orphans")
    }
    
    
    if (nrow(synapse_dropped)>0) {
      ##dropped "synapses"
      synapse_drop <- synapse_dropped %>% ungroup %>% filter(Label==labs[i]) %>% arrange(X,Y)
      synapse_drop_id <- 1:nrow(synapse_drop)
      synapse_drop$id <- synapse_drop_id
      
      assign(paste('synapse_drop_', labs[i], sep = ""), synapse_drop, envir = global_env())
      
      ##extract patches of dropped synapses
      syndrp_extract_pre <- extract_patches3D(imf1, cx = synapse_drop$X, cy = synapse_drop$Y,
                                              cz = synapse_drop$Slice*mic_pix, wx = 1.5/mic_pix, wy = 1.5/mic_pix, wz = 2) ##convert the slice back to integer
      
      syndrp_extract_post <- extract_patches3D(imf2, cx = synapse_drop$X, cy = synapse_drop$Y,
                                               cz = synapse_drop$Slice*mic_pix, wx = 1.5/mic_pix, wy = 1.5/mic_pix, wz = 2)
      
      ##convert image list to a data frame
      df_synd_pre <- as.data.frame(syndrp_extract_pre)
      
      ##convert to a tibble for tidy transformations
      df_synd_pre_tb <- as_tibble(df_synd_pre)
      
      ##cast the z stack around puncta of interest for MIP calculation
      df_synd_pre_tb2 <- dcast(df_synd_pre_tb, im+x+y~z, value.var = "value") %>%
        mutate(max = pmax(`1`, `2`, `3`)) %>%
        dplyr::select(-`1`, -`2`, -`3`)
      
      ##convert im from the image list to a number
      df_synd_pre_tb2$im2 <- as.numeric(df_synd_pre_tb2$im)
      
      ##make an id column to match the id of the puncta population(synapse, dropped synapse, or orphan)
      df_synd_pre_tb2 <- df_synd_pre_tb2 %>% mutate(id = im2)
      
      synd_ihc_id <- synapse_drop %>% ungroup %>% dplyr::select(id, ihc)
      
      ##join the ihc info to the image table
      df_synd_pre_tb2 <- left_join(df_synd_pre_tb2, synd_ihc_id)
      
      ##make text labels of ihc number to go with each id'd image
      synd_dat_text <- data.frame(
        label = as.character(synd_ihc_id$ihc),
        im2   = as.numeric(synd_ihc_id$id)
      )
      
      ##the facet for the puncta subtype (synapse, dropped synapse, or orphan)
      n_facet_synd <- round(sqrt(max(unique(df_synd_pre_tb2$im2))))
      
      ##make the ggplot and save
      ggplot(data = df_synd_pre_tb2, mapping = aes(x, y)) +
        geom_raster(aes(fill=max)) +
        scale_y_continuous(trans=scales::reverse_trans()) +
        geom_text(data = synd_dat_text, mapping = aes(x = 20, y = 23, label = label), 
                  hjust   = 0, vjust   = 0, size = 3, color = "white") +
        facet_wrap(~im2, ncol = n_facet_synd+2) +
        theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          strip.text.x = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
          strip.text.y = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
          plot.title = element_text(size = 14)
        ) +
        scale_fill_gradient(low="black", high="magenta") +
        coord_fixed(ratio = 1) +
        ggtitle(paste(labs[i], "_pre_dropped (x<=0.35 square microns)", sep = ""))
      
      ggsave(paste(labs[i], "_pre_dropped.tiff", sep = ""), path = "patches/pre/", width = (n_facet_synd+2)/1.7, height = ((n_facet_synd+2)/1.7)-1, units = "in")
      
      ##now working with the post channel
      ##there are less steps because geom assessts are reused from above
      df_synd_post <- as.data.frame(syndrp_extract_post)
      
      df_synd_post_tb <- as_tibble(df_synd_post)
      
      df_synd_post_tb2 <- dcast(df_synd_post_tb, im+x+y~z, value.var = "value") %>%
        mutate(max = pmax(`1`, `2`, `3`))
      
      df_synd_post_tb2$im2 <- as.numeric(df_synd_post_tb2$im)
      
      ggplot(data = df_synd_post_tb2, mapping = aes(x, y)) +
        geom_raster(aes(fill=max)) +
        scale_y_continuous(trans=scales::reverse_trans()) +
        geom_text(data = synd_dat_text, mapping = aes(x = 20, y = 23, label = label), 
                  hjust   = 0, vjust   = 0, size = 3, color = "white") +
        facet_wrap(~im2, ncol = n_facet_synd+2) +
        theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          strip.text.x = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
          strip.text.y = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
          plot.title = element_text(size = 14)
        ) +
        scale_fill_gradient(low="black", high="green") +
        coord_fixed(ratio = 1) +
        ggtitle(paste(labs[i], "_post_dropped (within 2 microns of pre)", sep = ""))
      
      ggsave(paste(labs[i], "_post_dropped.tiff", sep = ""), path = "patches/post/", width = (n_facet_synd+2)/1.7, height = ((n_facet_synd+2)/1.7)-1, units = "in")
    } else {
      print("there were no synapses below .35 square microns")
    }
  }
}

#SynapseJ_patches(all_puncta_1000_22.6_2)

##this is the same function but does not write an image
##the function is to regenerate puncta classification tibbles to the R global environment
SynapseJ_files <- function(x){
  ##working with CtBP2 puncta detected by SynapseJ for all figures
  data <- as_tibble(expr(!!x))
  
  detected_3d_xy <- data %>% filter(desc!="manual"&channel=="pre")
  #filter(norm_Y<500) ##my manual counts did not extend so far from the nucleus
  ##though perhaps a more precise filtering will be to group by label and ihc and take everything below the lowest synapse
  
  detected_lob_key <- left_join(detected_3d_xy, detected_3d_xy %>% filter(desc!="manual"&channel=="pre") %>%
                                  #filter(norm_Y<500) ##my manual counts did not extend so far from the nucleus
                                  ##though perhaps a more precise filtering will be to group by label and ihc and take everything below the lowest synapse
                                  group_by(Label, ihc) %>%
                                  filter(dist<40) %>%
                                  #filter(Area>0.35) %>% #obtaining the limit to "search" for synapses based on a typical minimum size (area) exclusion
                                  mutate(lowest_ob = max(norm_Y)) %>%
                                  ungroup)
  
  detected_lob_key <- unique(drop_na(detected_lob_key %>% dplyr::select(Label, ihc, lowest_ob)))
  
  detected_3d_xy <- left_join(detected_3d_xy, detected_lob_key)
  
  detected_3d_xy <- detected_3d_xy %>% 
    group_by(Label, ihc) %>% 
    filter(norm_Y<lowest_ob) %>%
    ungroup %>%
    dplyr::select(!lowest_ob)
  
  labs <- unique(detected_3d_xy %>% ungroup %>% dplyr::select(Label))
  
  labs <- labs$Label
  
  ##size filtering
  detected_3d_xy_szf <- detected_3d_xy %>% filter(Area>0.35)
  
  ##the puncta which were dropped
  detected_3d_xy_filtered <- anti_join(detected_3d_xy, detected_3d_xy_szf)
  
  ##any synapses (within 2 microns of GluA2) that were dropped? We would like to see those.
  synapse_dropped <- detected_3d_xy_filtered %>% filter(dist<40)
  
  
  for (i in 1:length(labs)) {
    ##now, iteration over each Label
    detected <- detected_3d_xy_szf %>% ungroup %>% filter(Label==labs[i])
    
    ##detected synapses with area>0.35 square microns
    synapse <- detected %>% ungroup %>% filter(dist<40) %>% arrange(X,Y)
    synapse_id <- 1:nrow(synapse)
    synapse$id <- synapse_id
    
    assign(paste('synapse_', labs[i], sep = ""), synapse, envir = global_env())
    
    
    # ##load the tiff image
    # img <- load.image(paste(path, labs[i], "/", labs[i], ".tif", sep = ""))
    # 
    # z_dim <- dim(img)[3] ##return the z stack dimension
    # which(1:z_dim%%2!=1) ##used below: split the z dimension where odd vals and even vals each correspond to a channel
    # 
    # ##odd values comprise the pre channel
    # imf1 <- frame(img, which(1:z_dim%%2==1))
    # 
    # ##even values comprise the post channel
    # imf2 <- frame(img, which(1:z_dim%%2!=1))
    # 
    # 
    # ##extract patches of synapses
    # synapse_extract_pre <- extract_patches3D(imf1, cx = synapse$X, cy = synapse$Y,
    #                                          cz = synapse$Slice*0.0593047, wx = 25, wy = 25, wz = 2) ##convert the slice back to integer
    # 
    # synapse_extract_post <- extract_patches3D(imf2, cx = synapse$X, cy = synapse$Y,
    #                                           cz = synapse$Slice*0.0593047, wx = 25, wy = 25, wz = 2)
    # 
    # ##convert image list to a data frame
    # df_pre <- as.data.frame(synapse_extract_pre)
    # 
    # ##convert to a tibble for tidy transformations
    # df_pre_tb <- as_tibble(df_pre)
    # 
    # ##cast the z stack around puncta of interest for MIP calculation
    # df_pre_tb2 <- dcast(df_pre_tb, im+x+y~z, value.var = "value") %>%
    #   mutate(max = pmax(`1`, `2`, `3`)) %>%
    #   dplyr::select(-`1`, -`2`, -`3`)
    # 
    # ##convert im from the image list to a number
    # df_pre_tb2$im2 <- as.numeric(df_pre_tb2$im)
    # 
    # ##make an id column to match the id of the puncta population(synapse, dropped synapse, or orphan)
    # df_pre_tb2 <- df_pre_tb2 %>% mutate(id = im2)
    # 
    # ihc_id <- synapse %>% ungroup %>% dplyr::select(id, ihc)
    # 
    # ##join the ihc info to the image table
    # df_pre_tb2 <- left_join(df_pre_tb2, ihc_id)
    # 
    # ##make text labels of ihc number to go with each id'd image
    # dat_text <- data.frame(
    #   label = as.character(ihc_id$ihc),
    #   im2   = as.numeric(ihc_id$id)
    # )
    # 
    # ##the facet for the puncta subtype (synapse, dropped synapse, or orphan)
    # n_facet <- round(sqrt(max(unique(df_pre_tb2$im2))))
    # 
    # ##make the ggplot and save
    # ggplot(data = df_pre_tb2, mapping = aes(x, y)) +
    #   geom_raster(aes(fill=max)) +
    #   scale_y_continuous(trans=scales::reverse_trans()) +
    #   geom_text(data = dat_text, mapping = aes(x = 20, y = 23, label = label), 
    #             hjust   = 0, vjust   = 0, size = 3, color = "white") +
    #   facet_wrap(~im2, ncol = n_facet+2) +
    #   theme(
    #     axis.text = element_blank(),
    #     axis.title = element_blank(),
    #     axis.ticks = element_blank(),
    #     strip.text.x = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
    #     strip.text.y = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
    #     plot.title = element_text(size = 14)
    #   ) +
    #   scale_fill_gradient(low="black", high="magenta") +
    #   coord_fixed(ratio = 1) +
    #   ggtitle(paste(labs[i], "_pre_synapse (0.35 < x <= 1.4 square microns)", sep = ""))
    # 
    # ggsave(paste(labs[i], "_pre_synapse.tiff", sep = ""), path = "patches/pre/", width = (n_facet+2)/1.7, height = ((n_facet+2)/1.7)-1, units = "in")
    # 
    ##now working with the post channel
    ##there are less steps because geom assets are reused from above
    # df_post <- as.data.frame(synapse_extract_post)
    # 
    # df_post_tb <- as_tibble(df_post)
    # 
    # df_post_tb2 <- dcast(df_post_tb, im+x+y~z, value.var = "value") %>%
    #   mutate(max = pmax(`1`, `2`, `3`))
    # 
    # df_post_tb2$im2 <- as.numeric(df_post_tb2$im)
    # 
    # ggplot(data = df_post_tb2, mapping = aes(x, y)) +
    #   geom_raster(aes(fill=max)) +
    #   scale_y_continuous(trans=scales::reverse_trans()) +
    #   geom_text(data = dat_text, mapping = aes(x = 20, y = 23, label = label), 
    #             hjust   = 0, vjust   = 0, size = 3, color = "white") +
    #   facet_wrap(~im2, ncol = n_facet+2) +
    #   theme(
    #     axis.text = element_blank(),
    #     axis.title = element_blank(),
    #     axis.ticks = element_blank(),
    #     strip.text.x = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
    #     strip.text.y = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
    #     plot.title = element_text(size = 14)
    #   ) +
    #   scale_fill_gradient(low="black", high="green") +
    #   coord_fixed(ratio = 1) +
    #   ggtitle(paste(labs[i], "_post_synapse (within 2 microns of pre)", sep = ""))
    # 
    # ggsave(paste(labs[i], "_post_synapse.tiff", sep = ""), path = "patches/post/", width = (n_facet+2)/1.7, height = ((n_facet+2)/1.7)-1, units = "in")
    # 
    # 
    
    
    if (nrow(detected)>nrow(synapse)){
      ##limit orphan CtBP2 punctas to the size-filtered data set
      orphan <- anti_join(detected, synapse) %>% arrange(X,Y)
      orphan_id <- 1:nrow(orphan)
      orphan$id <- orphan_id
      
      assign(paste('orphan_', labs[i], sep = ""), orphan, envir = global_env())
      
      # ##extract patches of orphans
      # orphan_extract_pre <- extract_patches3D(imf1, cx = orphan$X, cy = orphan$Y,
      #                                         cz = orphan$Slice*0.0593047, wx = 25, wy = 25, wz = 2) ##convert the slice back to integer
      # 
      # orphan_extract_post <- extract_patches3D(imf2, cx = orphan$X, cy = orphan$Y,
      #                                          cz = orphan$Slice*0.0593047, wx = 25, wy = 25, wz = 2)
      # 
      # ##convert image list to a data frame
      # df_orph_pre <- as.data.frame(orphan_extract_pre)
      # 
      # ##convert to a tibble for tidy transformations
      # df_orph_pre_tb <- as_tibble(df_orph_pre)
      # 
      # ##cast the z stack around puncta of interest for MIP calculation
      # df_orph_pre_tb2 <- dcast(df_orph_pre_tb, im+x+y~z, value.var = "value") %>%
      #   mutate(max = pmax(`1`, `2`, `3`)) %>%
      #   dplyr::select(-`1`, -`2`, -`3`)
      # 
      # ##convert im from the image list to a number
      # df_orph_pre_tb2$im2 <- as.numeric(df_orph_pre_tb2$im)
      # 
      # ##make an id column to match the id of the puncta population(synapse, dropped synapse, or orphan)
      # df_orph_pre_tb2 <- df_orph_pre_tb2 %>% mutate(id = im2)
      # 
      # orph_ihc_id <- orphan %>% ungroup %>% dplyr::select(id, ihc)
      # 
      # ##join the ihc info to the image table
      # df_orph_pre_tb2 <- left_join(df_orph_pre_tb2, orph_ihc_id)
      # 
      # ##make text labels of ihc number to go with each id'd image
      # orph_dat_text <- data.frame(
      #   label = as.character(orph_ihc_id$ihc),
      #   im2   = as.numeric(orph_ihc_id$id)
      # )
      # 
      # ##the facet for the puncta subtype (synapse, dropped synapse, or orphan)
      # n_facet_orph <- round(sqrt(max(unique(df_orph_pre_tb2$im2))))
      # 
      # ##make the ggplot and save
      # ggplot(data = df_orph_pre_tb2, mapping = aes(x, y)) +
      #   geom_raster(aes(fill=max)) +
      #   scale_y_continuous(trans=scales::reverse_trans()) +
      #   geom_text(data = orph_dat_text, mapping = aes(x = 20, y = 23, label = label), 
      #             hjust   = 0, vjust   = 0, size = 3, color = "white") +
      #   facet_wrap(~im2, ncol = n_facet_orph+2) +
      #   theme(
      #     axis.text = element_blank(),
      #     axis.title = element_blank(),
      #     axis.ticks = element_blank(),
      #     strip.text.x = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
      #     strip.text.y = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
      #     plot.title = element_text(size = 14)
      #   ) +
      #   scale_fill_gradient(low="black", high="magenta") +
      #   coord_fixed(ratio = 1) +
      #   ggtitle(paste(labs[i], "_pre_orphan (0.35 <= x <= 1.4 square microns)", sep = ""))
      # 
      # ggsave(paste(labs[i], "_pre_orphan.tiff", sep = ""), path = "patches/pre/", width = (n_facet_orph+2)/1.7, height = ((n_facet_orph+2)/1.7)-1, units = "in")
      # 
      # ##now working with the post channel
      ##there are less steps because geom assessts are reused from above
      #   df_orph_post <- as.data.frame(orphan_extract_post)
      #   
      #   df_orph_post_tb <- as_tibble(df_orph_post)
      #   
      #   df_orph_post_tb2 <- dcast(df_orph_post_tb, im+x+y~z, value.var = "value") %>%
      #     mutate(max = pmax(`1`, `2`, `3`))
      #   
      #   df_orph_post_tb2$im2 <- as.numeric(df_orph_post_tb2$im)
      #   
      #   ggplot(data = df_orph_post_tb2, mapping = aes(x, y)) +
      #     geom_raster(aes(fill=max)) +
      #     scale_y_continuous(trans=scales::reverse_trans()) +
      #     geom_text(data = orph_dat_text, mapping = aes(x = 20, y = 23, label = label), 
      #               hjust   = 0, vjust   = 0, size = 3, color = "white") +
      #     facet_wrap(~im2, ncol = n_facet_orph+2) +
      #     theme(
      #       axis.text = element_blank(),
      #       axis.title = element_blank(),
      #       axis.ticks = element_blank(),
      #       strip.text.x = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
      #       strip.text.y = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
      #       plot.title = element_text(size = 14)
      #     ) +
      #     scale_fill_gradient(low="black", high="green") +
      #     coord_fixed(ratio = 1) +
      #     ggtitle(paste(labs[i], "_post_orphan (within 2 microns of pre)", sep = ""))
      #   
      #   ggsave(paste(labs[i], "_post_orphan.tiff", sep = ""), path = "patches/post/", width = (n_facet_orph+2)/1.7, height = ((n_facet_orph+2)/1.7)-1, units = "in")
    } else {
      print("there were no orphans")
    }
    
    
    if (nrow(synapse_dropped)>0) {
      ##dropped "synapses"
      synapse_drop <- synapse_dropped %>% ungroup %>% filter(Label==labs[i]) %>% arrange(X,Y)
      synapse_drop_id <- 1:nrow(synapse_drop)
      synapse_drop$id <- synapse_drop_id
      
      assign(paste('synapse_drop_', labs[i], sep = ""), synapse_drop, envir = global_env())
      
      #   ##extract patches of dropped synapses
      #   syndrp_extract_pre <- extract_patches3D(imf1, cx = synapse_drop$X, cy = synapse_drop$Y,
      #                                           cz = synapse_drop$Slice*0.0593047, wx = 25, wy = 25, wz = 2) ##convert the slice back to integer
      #   
      #   syndrp_extract_post <- extract_patches3D(imf2, cx = synapse_drop$X, cy = synapse_drop$Y,
      #                                            cz = synapse_drop$Slice*0.0593047, wx = 25, wy = 25, wz = 2)
      #   
      #   ##convert image list to a data frame
      #   df_synd_pre <- as.data.frame(syndrp_extract_pre)
      #   
      #   ##convert to a tibble for tidy transformations
      #   df_synd_pre_tb <- as_tibble(df_synd_pre)
      #   
      #   ##cast the z stack around puncta of interest for MIP calculation
      #   df_synd_pre_tb2 <- dcast(df_synd_pre_tb, im+x+y~z, value.var = "value") %>%
      #     mutate(max = pmax(`1`, `2`, `3`)) %>%
      #     dplyr::select(-`1`, -`2`, -`3`)
      #   
      #   ##convert im from the image list to a number
      #   df_synd_pre_tb2$im2 <- as.numeric(df_synd_pre_tb2$im)
      #   
      #   ##make an id column to match the id of the puncta population(synapse, dropped synapse, or orphan)
      #   df_synd_pre_tb2 <- df_synd_pre_tb2 %>% mutate(id = im2)
      #   
      #   synd_ihc_id <- synapse_drop %>% ungroup %>% dplyr::select(id, ihc)
      #   
      #   ##join the ihc info to the image table
      #   df_synd_pre_tb2 <- left_join(df_synd_pre_tb2, synd_ihc_id)
      #   
      #   ##make text labels of ihc number to go with each id'd image
      #   synd_dat_text <- data.frame(
      #     label = as.character(synd_ihc_id$ihc),
      #     im2   = as.numeric(synd_ihc_id$id)
      #   )
      #   
      #   ##the facet for the puncta subtype (synapse, dropped synapse, or orphan)
      #   n_facet_synd <- round(sqrt(max(unique(df_synd_pre_tb2$im2))))
      #   
      #   ##make the ggplot and save
      #   ggplot(data = df_synd_pre_tb2, mapping = aes(x, y)) +
      #     geom_raster(aes(fill=max)) +
      #     scale_y_continuous(trans=scales::reverse_trans()) +
      #     geom_text(data = synd_dat_text, mapping = aes(x = 20, y = 23, label = label), 
      #               hjust   = 0, vjust   = 0, size = 3, color = "white") +
      #     facet_wrap(~im2, ncol = n_facet_synd+2) +
      #     theme(
      #       axis.text = element_blank(),
      #       axis.title = element_blank(),
      #       axis.ticks = element_blank(),
      #       strip.text.x = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
      #       strip.text.y = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
      #       plot.title = element_text(size = 14)
      #     ) +
      #     scale_fill_gradient(low="black", high="magenta") +
      #     coord_fixed(ratio = 1) +
      #     ggtitle(paste(labs[i], "_pre_dropped (x<=0.35 square microns)", sep = ""))
      #   
      #   ggsave(paste(labs[i], "_pre_dropped.tiff", sep = ""), path = "patches/pre/", width = (n_facet_synd+2)/1.7, height = ((n_facet_synd+2)/1.7)-1, units = "in")
      #   
      #   ##now working with the post channel
      #   ##there are less steps because geom assessts are reused from above
      #   df_synd_post <- as.data.frame(syndrp_extract_post)
      #   
      #   df_synd_post_tb <- as_tibble(df_synd_post)
      #   
      #   df_synd_post_tb2 <- dcast(df_synd_post_tb, im+x+y~z, value.var = "value") %>%
      #     mutate(max = pmax(`1`, `2`, `3`))
      #   
      #   df_synd_post_tb2$im2 <- as.numeric(df_synd_post_tb2$im)
      #   
      #   ggplot(data = df_synd_post_tb2, mapping = aes(x, y)) +
      #     geom_raster(aes(fill=max)) +
      #     scale_y_continuous(trans=scales::reverse_trans()) +
      #     geom_text(data = synd_dat_text, mapping = aes(x = 20, y = 23, label = label), 
      #               hjust   = 0, vjust   = 0, size = 3, color = "white") +
      #     facet_wrap(~im2, ncol = n_facet_synd+2) +
      #     theme(
      #       axis.text = element_blank(),
      #       axis.title = element_blank(),
      #       axis.ticks = element_blank(),
      #       strip.text.x = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
      #       strip.text.y = element_text(size = 8, margin = margin(.08,0,.08,0, "mm")),
      #       plot.title = element_text(size = 14)
      #     ) +
      #     scale_fill_gradient(low="black", high="green") +
      #     coord_fixed(ratio = 1) +
      #     ggtitle(paste(labs[i], "_post_dropped (within 2 microns of pre)", sep = ""))
      #   
      #   ggsave(paste(labs[i], "_post_dropped.tiff", sep = ""), path = "patches/post/", width = (n_facet_synd+2)/1.7, height = ((n_facet_synd+2)/1.7)-1, units = "in")
    } else {
      print("there were no synapses below .35 square microns")
    }
  }
}




