# TD Aim 2 BTV Americas Cost Affinity Matrix MDS Plots
## Originally Made by TD 07.09.25
## Edited for final version on 10.14.25
## Edited to include upper and lower on 10.16.26
#######################################################

#Load Required Libraries
library(cluster)
library(stringr)
library(viridis)
library(ape)
# library(treespace)
library(tidyverse)
library(ggplot2)
library(writexl)

#######################################################

# CLUSTER SUPPORT COST
#######################################################
# code ran in terminal
# cluster_matrix --average BTV_Americas_s1_aln_midpoint_root_nexus.tre BTV_Americas_s2_aln_midpoint_root_nexus.tre BTV_Americas_s3_aln_midpoint_root_nexus.tre BTV_Americas_s4_aln_midpoint_root_nexus.tre BTV_Americas_s5_aln_midpoint_root_nexus.tre BTV_Americas_s6_aln_midpoint_root_nexus.tre BTV_Americas_s7_aln_midpoint_root_nexus.tre BTV_Americas_s8_aln_midpoint_root_nexus.tre BTV_Americas_s9_aln_midpoint_root_nexus.tre BTV_Americas_s10_aln_midpoint_root_nexus.tre BTV_Americas_cluster_support_matrix_062526.png --cost cluster_support --csv_output BTV_Americas_cluster_support_matrix_062526.csv
####################################################### 


# SYMMETRICAL MATRIX
#######################################################
# read in data and turn it into a matrix. 
segment_data <- read.csv("BTV_Americas_cluster_support_matrix_062526_CLEAN.csv", 
                         row.names = 1, check.names = FALSE)
segment_matrix <- as.matrix(segment_data)

# get the averages to create a symmetrical matrix
segment_matrix_symmetrical <- ((segment_matrix + t(segment_matrix))/2)

# perform MDS analysis
mds.out <- cmdscale(segment_matrix_symmetrical)
save(mds.out,segment_matrix_symmetrical,file="MDS_objects.RData")

# turn mds.upper.out into a data frame
mds.out <- as.data.frame(mds.out)

# make segment number an additional row
mds.out$Segment <- rownames(mds.out)

# define mapping to unique segments
segment_levels <- unique(mds.out$Segment)

# include colors for all 10 segments
my_colors <- c("#88CCEE", "#44A499", "#117633", "#332288", "#DDCC77",
               "#999933", "#CC6677", "#862254", "#AA4499", "#A6A6A6")

# define which segment numbers get shape 16 and which get 17
shape_16_segments <- c(1, 2, 3, 4, 6, 7, 9)
shape_17_segments <- setdiff(1:10, shape_16_segments)

# function to get shape by segment name
get_shape <- function(segment_name) {
  # Extract the segment number from the name, e.g. "Segment 3" -> 3
  seg_num <- as.integer(sub("Segment ", "", segment_name))
  if (seg_num %in% shape_16_segments) {
    return(16)
  } else if (seg_num %in% shape_17_segments) {
    return(17)
  } else {
    return(16)  # default shape if not matched
  }
}

# apply shape function to all segments in data
shape_map <- sapply(segment_levels, get_shape)
names(shape_map) <- segment_levels

# create named lookup vectors
color_map <- setNames(my_colors, segment_levels)

# plot MDS data (based off code from Sean C.)
# if (save.plot) pdf("./mds_plot_upper.pdf")
cols <- plasma(10)
n.samples <- nrow(mds.out)/10
plot(mds.out[,1], mds.out[,2],
     xlab = "x", ylab = "y",
     col = color_map[mds.out$Segment],
     pch = shape_map[mds.out$Segment],
     cex = 2.5)
legend("topleft",
       legend = segment_levels,
       col = color_map[segment_levels],
       pch = shape_map[segment_levels],
       pt.cex = 2,
       bty = "n",
       ncol=2)
legend("topright",
       legend = c("Structural", "Non-structural"),
       pch = c(16, 17),
       pt.cex = 2,
       bty = "n")

# SAVED AS PNG IMAGE WITH THE FOLLOWING: WIDTH = 628, HEIGHT = 532 # 

####################################################### END SYMMETRICAL MATRIX


















####################################################### CLUSTER SUPPORT
####################################################### UPPER MATRIX
# read in data and turn it into a matrix. 
segment_data <- read.csv("BTV_Americas_cluster_support_matrix_062526_CLEAN.csv", 
                         row.names = 1, check.names = FALSE)
segment_matrix <- as.matrix(segment_data)

# pull out the upper triangle of the matrix and mirror the upper to the lower
upper_tri <- segment_matrix
upper_tri[lower.tri(upper_tri)] <- t(upper_tri)[lower.tri(upper_tri)]

# rename the triangles for clarity
matrix_upper_sym <- upper_tri

# perform MDS analysis
mds.upper.out <- cmdscale(matrix_upper_sym)
save(mds.upper.out,matrix_upper_sym,file="MDS_objects_upper.RData")

# turn mds.upper.out into a data frame
mds.upper.out <- as.data.frame(mds.upper.out)

# make segment number an additional row
mds.upper.out$Segment <- rownames(mds.upper.out)

# define mapping to unique segments
segment_levels <- unique(mds.upper.out$Segment)

# include colors for all 10 segments
my_colors <- c("#88CCEE", "#44A499", "#117633", "#332288", "#DDCC77",
               "#999933", "#CC6677", "#862254", "#AA4499", "#A6A6A6")

# define which segment numbers get shape 16 and which get 17
shape_16_segments <- c(1, 2, 3, 4, 6, 7, 9)
shape_17_segments <- setdiff(1:10, shape_16_segments)

# function to get shape by segment name
get_shape <- function(segment_name) {
  # Extract the segment number from the name, e.g. "Segment 3" -> 3
  seg_num <- as.integer(sub("Segment ", "", segment_name))
  if (seg_num %in% shape_16_segments) {
    return(16)
  } else if (seg_num %in% shape_17_segments) {
    return(17)
  } else {
    return(16)  # default shape if not matched
  }
}

# apply shape function to all segments in data
shape_map <- sapply(segment_levels, get_shape)
names(shape_map) <- segment_levels

# create named lookup vectors
color_map <- setNames(my_colors, segment_levels)

# plot MDS data (based off code from Sean C.)
# if (save.plot) pdf("./mds_plot_upper.pdf")
cols <- plasma(10)
n.samples <- nrow(mds.upper.out)/10
plot(mds.upper.out[,1], mds.upper.out[,2],
     xlab = "x", ylab = "y",
     col = color_map[mds.upper.out$Segment],
     pch = shape_map[mds.upper.out$Segment],
     cex = 2.5)
legend("topleft",
       legend = segment_levels,
       col = color_map[segment_levels],
       pch = shape_map[segment_levels],
       pt.cex = 2,
       bty = "n",
       ncol=2)
legend("right",
       legend = c("Structural", "Non-structural"),
       pch = c(16, 17),
       pt.cex = 2,
       bty = "n")

# SAVED AS PNG IMAGE WITH THE FOLLOWING: WIDTH = 628, HEIGHT = 532 # 
####################################################### END UPPER MATRIX





####################################################### LOWER MATRIX
# pull out the lower triangle of the matrix and mirror the lower to the upper
lower_tri <- segment_matrix
lower_tri[upper.tri(lower_tri)] <- t(lower_tri)[upper.tri(lower_tri)]

# rename the triangles for clarity
matrix_lower_sym <- lower_tri

# perform MDS analysis
mds.lower.out <- cmdscale(matrix_lower_sym)
save(mds.lower.out,matrix_lower_sym,file="MDS_objects_lower.RData")


# turn mds.upper.out into a data frame
mds.lower.out <- as.data.frame(mds.lower.out)

# make segment number an additional row
mds.lower.out$Segment <- rownames(mds.lower.out)

# define mapping to unique segments
segment_levels <- unique(mds.lower.out$Segment)

# include colors for all 10 segments
# my_colors <- c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b",
#                "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2")

# my_colors <- c("#dd75ce", "#422759", "#2444a2", "#85d6ff", "#79e2a8",
#                "#9b9f7e", "#dad30f", "#e8a237", "#da7d4a", "#ff3934")

my_colors <- c("#88CCEE", "#44A499", "#117633", "#332288", "#DDCC77",
               "#999933", "#CC6677", "#862254", "#AA4499", "#A6A6A6")

# define which segment numbers get shape 16 and which get 17
shape_16_segments <- c(1, 2, 3, 4, 6, 7, 9)
shape_17_segments <- setdiff(1:10, shape_16_segments)

# function to get shape by segment name
get_shape <- function(segment_name) {
  # Extract the segment number from the name, e.g. "Segment 3" -> 3
  seg_num <- as.integer(sub("Segment ", "", segment_name))
  if (seg_num %in% shape_16_segments) {
    return(16)
  } else if (seg_num %in% shape_17_segments) {
    return(17)
  } else {
    return(16)  # default shape if not matched
  }
}

# apply shape function to all segments in data
shape_map <- sapply(segment_levels, get_shape)
names(shape_map) <- segment_levels

# create named lookup vectors
color_map <- setNames(my_colors, segment_levels)

# plot MDS data (based off code from Sean C.)
# if (save.plot) pdf("./mds_plot_upper.pdf")
cols <- plasma(10)
n.samples <- nrow(mds.lower.out)/10
plot(mds.lower.out[,1], mds.lower.out[,2],
     xlab = "x", ylab = "y",
     col = color_map[mds.lower.out$Segment],
     pch = shape_map[mds.lower.out$Segment],
     cex = 2.5)
legend("topright",
       legend = segment_levels,
       col = color_map[segment_levels],
       pch = shape_map[segment_levels],
       pt.cex = 2,
       bty = "n",
       ncol=2)
legend("topleft",
       legend = c("Structural", "Non-structural"),
       pch = c(16, 17),
       pt.cex = 2,
       bty = "n")

# SAVED AS PNG IMAGE WITH THE FOLLOWING: WIDTH = 628, HEIGHT = 532 # 
####################################################### END LOWER MATRIX
