library(tidyverse)
library(sf)
library(rlang)
library(reshape2)
library(ggthemes)
library(patchwork)

##function for data processing using imagej cell counter over an image z stack

##you need to know:
#the pixels per micron
#the file path to the folder with 1) manual counts with the ImageJ Cell Counter plugin,
  ##2) the line drawn through the IHC region, 3) the line drawn through the OSC region

##IMPORTANT: each line must be started, then ended, on the same side of the image for each image
##so, either left-to-right or top-to-bottom
##the line direction can change between images but must be the the same within each image

##the input is a mouse number

##the current output is a tibble (X[mouse#]_mid_counts) in the R environment
##where each manual count was defined as within the sensory epithelium (location="sensep")
##or outside the sensory epithelium (location="other")


imcell_region_counts <- function(x) {
  ##import data
  ##notice the x and y pixels are given in the cell counter dataset are not precise (no decimal values)
  ##so, we can go ahead and replace those X and Y columns with calculated pixel values
  
  pix_mic <- 1.4453
  
  path <- "demos/immune cells/"
  
  mid_c1_counts <- read_csv(paste(path , expr(!!x), "_mid_c1_counts.csv", sep=""),
                            col_types = cols(`C-pos` = col_skip(), 
                                             `Z-pos` = col_skip(), `T-pos` = col_skip())) %>%
    mutate(Type = "Iba1",
           Xp = round(`X(micron)`*pix_mic, digits = 3), #pix_mic pixels per micron
           Yp = round(`Y(micron)`*pix_mic, digits = 3))
  
  # mid_c2_counts <- read_csv(paste(path , expr(!!x), "_mid_c2_counts.csv", sep=""),
  #                           col_types = cols(`C-pos` = col_skip(), 
  #                                            `Z-pos` = col_skip(), `T-pos` = col_skip())) %>%
  #   mutate(Type = "CD45",
  #          Xp = round(`X(micron)`*pix_mic, digits = 3), #pix_mic pixels per micron
  #          Yp = round(`Y(micron)`*pix_mic, digits = 3))
  
  mid_counts <- bind_rows(mid_c1_counts)
  mid_ihc <- read_csv(paste(path, expr(!!x), "_mid_ihc.csv", sep=""))
  mid_osc <- read_csv(paste(path, expr(!!x), "_mid_osc.csv", sep=""))
  
  
  ##convert data to sf style
  mid_ihc_matrix <- as.matrix( mid_ihc)
  mid_ihc_ls <- st_linestring(x =  mid_ihc_matrix, dim = "XY")
  
  mid_osc_matrix <- as.matrix( mid_osc)
  mid_osc_ls <- st_linestring(x =  mid_osc_matrix, dim = "XY")
  
  mid_counts_sf <-  mid_counts %>% st_as_sf(coords = c("Xp", "Yp"))
  
  ##create edges to connect ihc and osc regions
  lead_edge <- bind_rows(head( mid_ihc, n = 1), head( mid_osc, n = 1))
  lead_edge_matrix <- as.matrix(lead_edge)
  lead_edge_ls <- st_linestring(x = lead_edge_matrix, dim = "XY")
  
  trail_edge <- bind_rows(tail( mid_ihc, n = 1), tail( mid_osc, n = 1))
  trail_edge_matrix <- as.matrix(trail_edge)
  trail_edge_ls <- st_linestring(x = trail_edge_matrix, dim = "XY")
  
  ##Make a sensory epithelium polygon with all four linestring edges
  sens_ep_polygon2 <- st_multilinestring(list(lead_edge_ls,  mid_ihc_ls, 
                                              trail_edge_ls,  mid_osc_ls))
  
  ##Convert to a polygon
  sens_ep_polygon3 <- st_polygonize(sens_ep_polygon2)
  
  ##now filter counts through the polygon
  within_result <- st_within(  mid_counts_sf, sens_ep_polygon3)
  within_result2 <- melt(within_result)
  points_within <-  mid_counts_sf[within_result2$L1, ]
  
  
  ##but we can make the region crisp by joining ihc and osc regions at their 
  ##common ends
  ##finding the nearest-neighbor endpoints of the OSC region to the IHC region
  ## to get new edges of the IHC region
  
  osc_lead <- head( mid_osc, n = 1)
  
  osc_lead2 <- osc_lead %>% st_as_sf(coords = c("X", "Y"))
  
  mid_ihc_mp <- st_cast(mid_ihc_ls, "MULTIPOINT")
  lead_edge <- st_nearest_points( mid_ihc_mp, osc_lead2)
  
  ihc_lead_point <- st_cast(lead_edge[[1]], "POINT")
  ihc_lead_point_dbl <- as.double(ihc_lead_point)
  ihc_lead_point_tbl <- tibble(
    X = as.double(ihc_lead_point_dbl[1]),
    Y = ihc_lead_point_dbl[2]
  )
  
  mid_ihc_tb <- tibble(
    X = mid_ihc_matrix[,1],
    Y = mid_ihc_matrix[,2],
    index = 1:length(mid_ihc_matrix[,1])
  )
  
  start <- inner_join(mid_ihc_tb, ihc_lead_point_tbl)
  
  ihc_lead_edge <-  mid_ihc_ls[start$index:nrow( mid_ihc_ls), ]
  
  ihc_lead_edge_ls <- st_linestring(ihc_lead_edge, dim = "X,Y")
  
  ##now for the new trailing edge of the IHC region
  ##a better way to specify the end of the ihc line is by index instead of X location
  ##because the line can double back (multiple line points for a single X value)
  osc_trail <- tail( mid_osc, n = 1)
  
  osc_trail2 <- osc_trail %>% st_as_sf(coords = c("X", "Y"))
  
  mid_ihc_mp <- st_cast(mid_ihc_ls, "MULTIPOINT")
  trail_edge <- st_nearest_points(mid_ihc_mp, osc_trail2)
  
  ihc_trail_point <- st_cast(trail_edge[[1]], "POINT")
  ihc_trail_point_dbl <- as.double(ihc_trail_point)
  ihc_trail_point_tbl <- tibble(
    X = ihc_trail_point_dbl[1],
    Y = ihc_trail_point_dbl[2]
  )
  
  ihc_lead_edge_tb <- tibble(
    X = ihc_lead_edge[,1],
    Y = ihc_lead_edge[,2],
    index = 1:length(ihc_lead_edge[,1])
  )
  
  end <- inner_join(ihc_lead_edge_tb, ihc_trail_point_tbl)
  
  ihc_trail_edge <- ihc_lead_edge[1:end$index, ]
  
  ihc_trail_edge_ls <- st_linestring(ihc_trail_edge, dim = "X,Y")
  
  ##make new edges to connect the osc region with the new ihc region
  ##ihc_trail_edge is a combination of lead and trail modifications to the ihc region
  new_ihc_lead_edge <- as_tibble(head(ihc_trail_edge, n=1))
  
  new_lead_edge <- bind_rows(new_ihc_lead_edge, head( mid_osc, n = 1))
  new_lead_edge_matrix <- as.matrix(new_lead_edge)
  new_lead_edge_ls <- st_linestring(x = new_lead_edge_matrix, dim = "XY")
  
  new_ihc_trail_edge <- as_tibble(tail(ihc_trail_edge, n=1))
  
  new_trail_edge <- bind_rows(new_ihc_trail_edge, tail( mid_osc, n = 1))
  new_trail_edge_matrix <- as.matrix(new_trail_edge)
  new_trail_edge_ls <- st_linestring(x = new_trail_edge_matrix, dim = "XY")
  
  
  sens_ep_polygon4 <- st_multilinestring(list(new_lead_edge_ls, ihc_trail_edge_ls, 
                                              new_trail_edge_ls,  mid_osc_ls))
  sens_ep_polygon5 <- st_polygonize(sens_ep_polygon4)
  
  plot(sens_ep_polygon3)
  plot(sens_ep_polygon5)
  
  assign(paste('X', expr(!!x), "_sens_ep", sep = ""), sens_ep_polygon5, envir = global_env())
  
  ##now filter counts through the new polygon
  within_result3 <- st_within(  mid_counts_sf, sens_ep_polygon5)
  within_result4 <- melt(within_result3)
  points_within2 <-  mid_counts[within_result4$L1, ]
  sensep_points <- inner_join(mid_counts, points_within2) %>%
    mutate(location = "sens_ep")
  other_points <- anti_join(mid_counts, points_within2) %>% 
    mutate(location = "other")
  mid_counts2 <- bind_rows(sensep_points, other_points) %>%
    mutate(mouse = expr(!!x),
           region = "mid",
           length = st_length(ihc_trail_edge_ls)/pix_mic) ##length along the inner hair cell region (convert to microns)
  
  assign(paste('X', expr(!!x), "_mid_counts", sep = ""), mid_counts2, envir = global_env())
}

#nums <- c(842, 843, 844, 845, 846, 848, 849, 850, 852, 853, 854, 855, 857, 858, 861, 847, 851, 856, 859)

#lapply(nums, imcell_region_counts)

