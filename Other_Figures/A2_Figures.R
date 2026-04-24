#######################################################
# BTV Americas - Figures: Maps, Pie Graphs, Bar Charts 
# TD 03.25.25
# TD Cleaned on 02.23.26
#######################################################
library(ggplot2)
library(tidyverse)
library(dplyr)
library(readxl)
library(tidyr)

library(ggrepel)  
library(ggmap)
library(maps)
library(mapdata)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)
#######################################################


### Unique Genotype-code Heatmap ### 
#######################################################
# read in excel spreadsheet
back_df <- read_excel("BTV_Americas_Metadata_030425.xlsx", sheet = "All Metadata")

# turn genogroup_code into factor, then order alphabetically by genogroup and then numerically by serotype
back_df$genogroup_code <- factor(back_df$genogroup_code)

# then order numerically by serotype and then alphabetically by genogroup
back_df$genogroup_code <- factor(back_df$genogroup_code, levels = c('AH01AAADAAAA','AH01AAADAAAC','AI02AAADAAAA',
                                                                    'AI02AAADAAAC', 'AI02AAADCAAA','AB03AAAEAAAA',
                                                                    'AB03AAAEAAAC','AB03AAAEBAAA','AB03AAAEBAAC',
                                                                    'AE05AAAEAAAC','AE05AAAEBAAC','AC06AAAEBAAA',
                                                                    'AC06AAAEBAAC','AE09AAADAAAA','AE09AAADAAAC',
                                                                    'AA10AAAFAAAA','AA10AAAFAAAC','AA11AAAFAAAA',
                                                                    'AA11AAAFAAAC','AG12AAAAAAAA','AG12AAAAAAAC',
                                                                    'AB13AAAEBAAA','AB13AAAEBAAC','AC14AAAEBAAA',
                                                                    'AC14AAAEBAAC','AJ15AAABEAAC','AA17AAAFAAAA',
                                                                    'AA17AAAFAAAC','AD18AAADAAAA','AD18AAADAAAC',
                                                                    'AF19AAACDAAA','AF19AAACDAAC','AG22AAAAAAAA',
                                                                    'AG22AAAABAAA','AA24AAAFAAAA','AA24AAAFAAAC'), ordered = TRUE)

# turn country into factor, then order by the place with the most sequences (highest diversity)
back_df$country <- factor(back_df$country)
back_df$country <- factor(back_df$country, levels = c('United States','Cuba','Costa Rica','Guatemala','Panama','Honduras','Jamaica',
                                                      'Barbados','Ecuador','El Salvador','Nicaragua','Belize','Dominican Republic',
                                                      'Puerto Rico','Trinidad','Tobago','Brazil','Martinique'), ordered = TRUE)

# count the number of genotypes per country
counts_df <- back_df %>%
  group_by(genogroup_code, country) %>%
  summarise(count = n(), .groups = "drop")

# join the counts into the main df
back_df_counts <- back_df %>%
  left_join(counts_df, by = c("genogroup_code", "country"))

# turn serotype into a factor
back_df$serotype <- factor(back_df$serotype)

# total the counts to figure out which has the most
genogroup_totals <- counts_df %>%
  group_by(genogroup_code) %>%
  summarise(total_count = sum(count)) %>%
  arrange(desc(total_count))

# order the graph by the genogroup totals
back_df_counts$genogroup_code <- factor(back_df_counts$genogroup_code, 
                                        levels = genogroup_totals$genogroup_code)

# plot with genogroup totals order and NO COLOR
ggplot(back_df_counts, aes(x = genogroup_code, y = country)) + 
  geom_tile() + 
  geom_text(aes(label = count), color = "white", size = 3) +
  theme_bw()  +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Genotype") + ylab("Country") + labs(fill = "Count") + 
  ggtitle("Genotypes per Country")

# ggsave("Genogroup_country_heatmap_with_counts_ordered_NOCOLOR_122625.pdf", width=10, height=8, units="in")

#######################################################


### Genotype-code Heatmap - s2 order ### 
#######################################################

# read in excel spreadsheet
backbone_df <- read_excel("BTV_Americas_Metadata_030425.xlsx", sheet = "All Metadata")

# parse out genogroup into its own segment column 
parsed_df <- backbone_df %>%
  mutate(genogroup = str_split(genogroup_code, "")) %>%
  unnest_longer(genogroup) %>%
  filter(str_detect(genogroup, "[A-Za-z]")) %>%
  group_by(name, serotype, genogroup_code) %>%
  mutate(segment_number = row_number()) %>%
  ungroup()

# read in data from tempest organized by how tree is
s2_midroot_data <- read.csv("BTV_Americas_s2_aln_midpoint_root_color_by_genogroup_2.fasta_tempest_081525.txt", header = TRUE, sep = "\t")

# subset data - names only
s2_midroot_data_subset <- s2_midroot_data[,1, drop=FALSE]

# reverse order to match tree
s2_midroot_name_vector <- rev(s2_midroot_data_subset$Name)

# rename dataset
parsed_df_renamed_s2 <- parsed_df %>%
  mutate(name = str_replace_all(name, "[/|]", "_"))

# reorder dataset by imported data from tree
parsed_df_renamed_s2$name <- factor(parsed_df_renamed_s2$name, levels = s2_midroot_name_vector, ordered = TRUE)

# make heat map ordered by data from tree and vertical so it can get put next to tree
parsed_df_renamed_s2 %>% 
  count(name, segment_number, genogroup) %>%  
  ggplot(mapping = aes(x = segment_number, y = name)) +
  geom_tile(mapping = aes(fill = genogroup)) +
  geom_tile(color = "black", fill = NA, linewidth = 0.5) +
  scale_x_continuous(breaks = 1:10) + 
  scale_fill_discrete(name = "BTV-Glue Genogroup") + 
  xlab("Segment") + ylab("Isolate")

# ggsave("BTV_Americas_genogroup_heatmap_s2_treeorder_050225.pdf", height=19, width=15, units="in")
# ggsave("BTV_Americas_genogroup_heatmap_s2_treeorder_081525.pdf", height=19, width=15, units="in")

#######################################################


### MDS Visualization ###
#######################################################
#Load Required Libraries
library(cluster)
library(stringr)
library(viridis)
library(ape)
library(tidyverse)
library(ggplot2)
library(writexl)

# read in data and turn it into a matrix. 
segment_data <- read.csv("BTV_Americas_cluster_support_matrix_100725_tidy.csv", 
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
legend("bottomleft",
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

# export symmetrical matrix
# write_xlsx(as.data.frame(segment_matrix_symmetrical), "TD_CAM_Values_BTV_Americas_all_segments_midpoint_root_nexus_cluster_matrix_082225_tidy_SYMMETRICAL.xlsx")

#######################################################


### Map of collection locations ### 
#######################################################
# Get world and americas map from map_data package
world <- map_data("world")

americas <- subset(world, region %in% c("USA","Brazil","Mexico", "Colombia", "Argentina", "Canada",
                                        "Peru","Venezuela","Chile","Guatemala","Ecuador", "Bolivia", "Cuba",
                                        "Honduras", "Paraguay", "Nicaragua","El Salvador", "Costa Rica", "Panama",
                                        "Uruguay",  "Jamaica",  "Trinidad and Tobago", "Guyana", "Suriname", "Belize",
                                        "Barbados", "Saint Lucia", "Grenada", "Saint Vincent and the Grenadines", 
                                        "Antigua and Barbuda", "Saint Kitts and Nevis", "Dominican Republic"))

country_colors <- c("USA" = "blue", "Barbados" = "red", "Belize" = "red", "Brazil" = "red",
                    "Costa Rica" = "red", "Cuba" = "red", "Dominican Republic" = "red", "Ecuador" = "red", 
                    "El Salvador" = "red", "Guatemala" = "red", "Honduras" = "red", "Jamaica" = "red", 
                    "Martinique" = "red", "Nicaragua" = "red", "Panama" = "red", "Puerto Rico" = "red",
                    "Trinidad and Tobago" = "red")

# read in lat long data
All_Points <- read_excel("BTV_Americas_Metadata_030425.xlsx", sheet = "All lat long")

# this defines limits of the world map in long/lat that will be plotted
xlim = c(-180,180)
ylim = c(-55,75)

# Set fill color but don't include Alaska and Hawaii
americas_colored <- americas %>%
  mutate(fill_color = ifelse(region %in% names(country_colors) & !(region == "USA" & long < -125), 
                             country_colors[region], "grey90"))

# Plot the map with country-specific colors and single black points
americas_map_colored_points <- ggplot() + 
  geom_polygon(data = americas_colored, aes(x = long, y = lat, group = group, fill = fill_color), color = "grey70",  size = 0.2) +
  geom_point(data = All_Points, aes(x = long, y = lat), color = "black", size = 2, alpha = 0.8, inherit.aes = FALSE) +
  geom_text(data = All_Points, aes(x = long, y = lat, label = country), hjust = -0.1, vjust = 0, size = 3) +  
  scale_fill_identity() +
  coord_fixed(ratio = 1.1, xlim = c(-180, -35)) +
  theme_void()

americas_map_colored_points

# plot with multiple points for each sample.

# include sample metadata
All_Metadata <- read_excel("BTV_Americas_Metadata_030425.xlsx", sheet = "All Metadata")

# left join with lat long data
All_Metadata_and_Points <- All_Metadata %>%
  left_join(All_Points %>% select(county, state, country, continent, lat, long), 
            by = c("county", "state", "country", "continent"))

#create plot
ggplot() + 
  geom_polygon(data = americas_colored, aes(x = long, y = lat, group = group, fill = fill_color), color = "grey70",  size = 0.2) +
  geom_point(data = All_Metadata_and_Points, aes(x = long, y = lat), shape = 21, fill = NA, color = "black", size = 1, 
             alpha = 0.8,position = position_jitter(width = 0.8, height = 0.8)) +
  # geom_text(data = All_Metadata_and_Points, aes(x = long, y = lat, label = country), hjust = -0.1, vjust = 0, size = 3, check_overlap = TRUE) +  
  scale_fill_identity() +
  coord_fixed(ratio = 1.1, xlim = c(-180, -35)) +
  theme_void()

# ggsave("Americas_map_colored_with_points.pdf", width=7, height=7, units="in")

#######################################################