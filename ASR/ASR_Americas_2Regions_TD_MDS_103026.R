################################################################
# Aim 2 - Ancestral State Reconstruction (ASR)
# Using correct cluster files
# Americas Region Categories: Inside USA & Outside USA
# TD 01.08.26
# TD 06.23.26 Repeated with trees made with bootstrapping
# PREVIOUSLY: n = 163 (removed samples causing polytomies)
# TD 10.14.26 Repeated with trees with identical sequences
# ^ following reviewer comments from Oxford Academic
# ^ n = 169 (kept samples causing polytomies)
# ^ adding small branch lengths to resolve issues with this
# Updated model selection
# ^ Use SYM/ER unless delta AIC is greater than 2, then use ARD
################################################################

# load required libraries
library(tidyverse)
library(phytools)
library(geiger)
library(plotrix)
library(corHMM)

#######################################################

#### Segment 1 ####

########################
# SIMPLE REGIONS 
# AMERICAS
# SEGMENT 1
########################

# plot base tree
#######################################################
# read in tree
s1.tree<-read.tree(file="./trees/BTV_Americas_s1_aln_REPEAT.treefile")
s1.tree

# does this tree have any zero-length branches?
any(s1.tree$edge.length == 0)
# yes it does

plot_tree_with_colored_zero_length_branches <- function(tree) {
  
  # how many tips and nodes does this tree have?
  number_tips  = Ntip(tree)
  number_nodes = tree$Nnode
  
  # get indexes of all the zero-length edges
  zero_length_edges <- which(tree$edge.length == 0) 
  
  # all the downstream nodes/tips of zero-length branches
  # all of them
  zero_length_children <- tree$edge[tree$edge.length == 0, 2]
  # just the tips
  zero_length_tips  <- zero_length_children[zero_length_children <= number_tips]
  # just the nodes
  zero_length_nodes <- zero_length_children[zero_length_children > number_tips]
  
  # create color schemes
  tip_colors  = rep("black", number_tips)
  edge_colors = rep("black", length(tree$edge.length))
  edge_widths = rep(0.5, length(tree$edge.length))
  
  # color any 0-length tips or branches red, and thick
  tip_colors[zero_length_tips] <- "red"
  edge_colors[zero_length_edges] <- "red"
  edge_widths[zero_length_edges] <- 2
  
  # make the plot
  plot(tree, 
       tip.color = tip_colors, 
       edge.color = edge_colors, 
       edge.width = edge_widths, 
       cex = 0.25)
}

# plot s1.tree with tips downstream of 0-length branches colored red
plot_tree_with_colored_zero_length_branches(s1.tree)

# now convert 0-length branches into branches of length 1e-6
# using base R indexing
s1.tree.no_zero <- s1.tree
s1.tree.no_zero$edge.length[s1.tree.no_zero$edge.length == 0] <- 1e-6

# does this have any 0-length branches?
any(s1.tree.no_zero$edge.length == 0 )

# plot it
plot_tree_with_colored_zero_length_branches(s1.tree.no_zero)

# midpoint root tree
s1.tree.midpoint_rooted <- midpoint_root(s1.tree.no_zero)

# write out midpoint rooted trees in .nexus format
write.tree(s1.tree.midpoint_rooted, file="./trees/rooted/s1_tree_R_midpoint_rooted.newick")
 
# make sure that tip label name structure matches data
s1.tree.midpoint_rooted$tip.label <- gsub("/", "_", s1.tree.midpoint_rooted$tip.label)
s1.tree.midpoint_rooted$tip.label <- gsub("\\|", "_", s1.tree.midpoint_rooted$tip.label)

# read in metadata
s1.data<-read.csv(file="A2_Americas_Metadata.csv", row.names=1)

# keep only the tips that are found in the tree
s1.data.subset <- s1.data %>% filter(rownames(s1.data) %in% s1.tree.midpoint_rooted$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s1.tree.midpoint_rooted, s1.data.subset)
summary(chk)

# create a named vector
s1.trait <- setNames(s1.data.subset[,"trait"], rownames(s1.data.subset))
unique(s1.trait)

plotTree(s1.tree.midpoint_rooted,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s1.trait),names(s1.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s1.tree.midpoint_rooted),2,
             dimnames=list(s1.tree.midpoint_rooted$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s1.tree.midpoint_rooted)) 
  pies[s1.tree.midpoint_rooted$tip.label[i],
       X[[s1.tree.midpoint_rooted$tip.label[i]]]]<-
  rep(1/length(X[[s1.tree.midpoint_rooted$tip.label[i]]]),
      length(X[[s1.tree.midpoint_rooted$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s1.tree.midpoint_rooted$tip.label, tip = 1:Ntip(s1.tree.midpoint_rooted), adj = -0.05, 
          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(cols),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)

#######################################################
# fit model
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
s1.ard<-fitMk(s1.tree.midpoint_rooted,s1.trait,model="ARD")
s1.ard

### SYM
s1.sym<-fitMk(s1.tree.midpoint_rooted,s1.trait,model="SYM")
s1.sym

### ER
s1.er<-fitMk(s1.tree.midpoint_rooted,s1.trait,model="ER")
s1.er

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s1.ard),
                    logLik(s1.sym),
                    logLik(s1.er)),
           AIC=c(AIC(s1.ard),
                 AIC(s1.sym),
                 AIC(s1.er)))

# calculate delta AIC
aic_vals <- c(
  ER  = AIC(s1.er),
  SYM = AIC(s1.sym),
  ARD = AIC(s1.ard))
delta_aic <- aic_vals - min(aic_vals)
delta_aic

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit
# SYM/ER are equally supported so plot both

####################################################### 

# Marginal ASR SYM
####################################################### 
# marginal ancestral reconstruction
s1_marginal_sym<-ancr(s1.sym, type="marginal")
s1_marginal_sym
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s1.sym)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s1_QMatrix_americas_SYM_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s1.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s1_marginal_sym,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols

# SAVE AS: 15 x 18 pdf.
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s1.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s1_marginal_sym,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols

# SAVE AS: 10 x 8 pdf.
####################################################### 


# Marginal ASR ER
####################################################### 
# marginal ancestral reconstruction
s1_marginal_er<-ancr(s1.er, type="marginal")
s1_marginal_er
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s1.er)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s1_QMatrix_americas_ER_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s1.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s1_marginal_er,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s1.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s1_marginal_er,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

########################
#
# END SEGMENT 1
#
########################



















#### Segment 2 ####

########################
# SIMPLE REGIONS 
# AMERICAS
# SEGMENT 2
########################

# plot base tree
#######################################################
# read in tree
s2.tree<-read.tree(file="./trees/BTV_Americas_s2_aln.fasta.treefile")
s2.tree

# does this tree have any zero-length branches?
any(s2.tree$edge.length == 0)
# yes it does

plot_tree_with_colored_zero_length_branches <- function(tree) {
  
  # how many tips and nodes does this tree have?
  number_tips  = Ntip(tree)
  number_nodes = tree$Nnode
  
  # get indexes of all the zero-length edges
  zero_length_edges <- which(tree$edge.length == 0) 
  
  # all the downstream nodes/tips of zero-length branches
  # all of them
  zero_length_children <- tree$edge[tree$edge.length == 0, 2]
  # just the tips
  zero_length_tips  <- zero_length_children[zero_length_children <= number_tips]
  # just the nodes
  zero_length_nodes <- zero_length_children[zero_length_children > number_tips]
  
  # create color schemes
  tip_colors  = rep("black", number_tips)
  edge_colors = rep("black", length(tree$edge.length))
  edge_widths = rep(0.5, length(tree$edge.length))
  
  # color any 0-length tips or branches red, and thick
  tip_colors[zero_length_tips] <- "red"
  edge_colors[zero_length_edges] <- "red"
  edge_widths[zero_length_edges] <- 2
  
  # make the plot
  plot(tree, 
       tip.color = tip_colors, 
       edge.color = edge_colors, 
       edge.width = edge_widths, 
       cex = 0.25)
}

# plot s2.tree with tips downstream of 0-length branches colored red
plot_tree_with_colored_zero_length_branches(s2.tree)

# now convert 0-length branches into branches of length 1e-6
# using base R indexing
s2.tree.no_zero <- s2.tree
s2.tree.no_zero$edge.length[s2.tree.no_zero$edge.length == 0] <- 1e-6

# does this have any 0-length branches?
any(s2.tree.no_zero$edge.length == 0 )

# plot it
plot_tree_with_colored_zero_length_branches(s2.tree.no_zero)

# midpoint root tree
s2.tree.midpoint_rooted <- midpoint_root(s2.tree)

# write out midpoint rooted trees in .nexus format
write.tree(s2.tree.midpoint_rooted, file="./trees/rooted/s2_tree_R_midpoint_rooted.newick")

# make sure that tip label name structure matches data
s2.tree.midpoint_rooted$tip.label <- gsub("/", "_", s2.tree.midpoint_rooted$tip.label)
s2.tree.midpoint_rooted$tip.label <- gsub("\\|", "_", s2.tree.midpoint_rooted$tip.label)

# read in metadata
s2.data<-read.csv(file="A2_Americas_Metadata.csv", row.names=1)

# keep only the tips that are found in the tree
s2.data.subset <- s2.data %>% filter(rownames(s2.data) %in% s2.tree.midpoint_rooted$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s2.tree.midpoint_rooted, s2.data.subset)
summary(chk)

# create a named vector
s2.trait <- setNames(s2.data.subset[,"trait"], rownames(s2.data.subset))
unique(s2.trait)

plotTree(s2.tree.midpoint_rooted,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s2.trait),names(s2.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s2.tree.midpoint_rooted),2,
             dimnames=list(s2.tree.midpoint_rooted$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s2.tree.midpoint_rooted)) 
  pies[s2.tree.midpoint_rooted$tip.label[i],
       X[[s2.tree.midpoint_rooted$tip.label[i]]]]<-
  rep(1/length(X[[s2.tree.midpoint_rooted$tip.label[i]]]),
      length(X[[s2.tree.midpoint_rooted$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s2.tree.midpoint_rooted$tip.label, tip = 1:Ntip(s2.tree.midpoint_rooted), adj = -0.05, 
          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(cols),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)

#######################################################
# fit model
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
s2.ard<-fitMk(s2.tree.midpoint_rooted,s2.trait,model="ARD")
s2.ard

### SYM
s2.sym<-fitMk(s2.tree.midpoint_rooted,s2.trait,model="SYM")
s2.sym

### ER
s2.er<-fitMk(s2.tree.midpoint_rooted,s2.trait,model="ER")
s2.er

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s2.ard),
                    logLik(s2.sym),
                    logLik(s2.er)),
           AIC=c(AIC(s2.ard),
                 AIC(s2.sym),
                 AIC(s2.er)))


# calculate delta AIC
aic_vals <- c(
  ER  = AIC(s2.er),
  SYM = AIC(s2.sym),
  ARD = AIC(s2.ard))
delta_aic <- aic_vals - min(aic_vals)
delta_aic

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit
# ARD fits best, only plot that 

####################################################### 

# Marginal ASR
####################################################### 
# marginal ancestral reconstruction
s2_marginal<-ancr(s2.ard, type="marginal")
s2_marginal
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s2.ard)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s2_QMatrix_americas_ARD_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s2.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s2_marginal,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s2.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s2_marginal,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

########################
#
# END SEGMENT 2
#
########################



















#### Segment 3 ####

########################
# SIMPLE REGIONS 
# AMERICAS
# SEGMENT 3
########################

# plot base tree
#######################################################
# read in tree
s3.tree<-read.tree(file="./trees/BTV_Americas_s3_aln_REPEAT.fasta.treefile")
s3.tree

# does this tree have any zero-length branches?
any(s3.tree$edge.length == 0)
# yes it does

plot_tree_with_colored_zero_length_branches <- function(tree) {
  
  # how many tips and nodes does this tree have?
  number_tips  = Ntip(tree)
  number_nodes = tree$Nnode
  
  # get indexes of all the zero-length edges
  zero_length_edges <- which(tree$edge.length == 0) 
  
  # all the downstream nodes/tips of zero-length branches
  # all of them
  zero_length_children <- tree$edge[tree$edge.length == 0, 2]
  # just the tips
  zero_length_tips  <- zero_length_children[zero_length_children <= number_tips]
  # just the nodes
  zero_length_nodes <- zero_length_children[zero_length_children > number_tips]
  
  # create color schemes
  tip_colors  = rep("black", number_tips)
  edge_colors = rep("black", length(tree$edge.length))
  edge_widths = rep(0.5, length(tree$edge.length))
  
  # color any 0-length tips or branches red, and thick
  tip_colors[zero_length_tips] <- "red"
  edge_colors[zero_length_edges] <- "red"
  edge_widths[zero_length_edges] <- 2
  
  # make the plot
  plot(tree, 
       tip.color = tip_colors, 
       edge.color = edge_colors, 
       edge.width = edge_widths, 
       cex = 0.25)
}

# plot s3.tree with tips downstream of 0-length branches colored red
plot_tree_with_colored_zero_length_branches(s3.tree)

# now convert 0-length branches into branches of length 1e-6
# using base R indexing
s3.tree.no_zero <- s3.tree
s3.tree.no_zero$edge.length[s3.tree.no_zero$edge.length == 0] <- 1e-6

# does this have any 0-length branches?
any(s3.tree.no_zero$edge.length == 0 )

# plot it
plot_tree_with_colored_zero_length_branches(s3.tree.no_zero)

# midpoint root tree
s3.tree.midpoint_rooted <- midpoint_root(s3.tree)

# write out midpoint rooted trees in .nexus format
write.tree(s3.tree.midpoint_rooted, file="./trees/rooted/s3_tree_R_midpoint_rooted.newick")

# make sure that tip label name structure matches data
s3.tree.midpoint_rooted$tip.label <- gsub("/", "_", s3.tree.midpoint_rooted$tip.label)
s3.tree.midpoint_rooted$tip.label <- gsub("\\|", "_", s3.tree.midpoint_rooted$tip.label)

# read in metadata
s3.data<-read.csv(file="A2_Americas_Metadata.csv", row.names=1)

# keep only the tips that are found in the tree
s3.data.subset <- s3.data %>% filter(rownames(s3.data) %in% s3.tree.midpoint_rooted$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s3.tree.midpoint_rooted, s3.data.subset)
summary(chk)

# create a named vector
s3.trait <- setNames(s3.data.subset[,"trait"], rownames(s3.data.subset))
unique(s3.trait)

plotTree(s3.tree.midpoint_rooted,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s3.trait),names(s3.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s3.tree.midpoint_rooted),2,
             dimnames=list(s3.tree.midpoint_rooted$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s3.tree.midpoint_rooted)) 
  pies[s3.tree.midpoint_rooted$tip.label[i],
       X[[s3.tree.midpoint_rooted$tip.label[i]]]]<-
  rep(1/length(X[[s3.tree.midpoint_rooted$tip.label[i]]]),
      length(X[[s3.tree.midpoint_rooted$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s3.tree.midpoint_rooted$tip.label, tip = 1:Ntip(s3.tree.midpoint_rooted), adj = -0.05, 
          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(cols),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)

#######################################################
# fit model
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
s3.ard<-fitMk(s3.tree.midpoint_rooted,s3.trait,model="ARD")
s3.ard

### SYM
s3.sym<-fitMk(s3.tree.midpoint_rooted,s3.trait,model="SYM")
s3.sym

### ER
s3.er<-fitMk(s3.tree.midpoint_rooted,s3.trait,model="ER")
s3.er

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s3.ard),
                    logLik(s3.sym),
                    logLik(s3.er)),
           AIC=c(AIC(s3.ard),
                 AIC(s3.sym),
                 AIC(s3.er)))


# calculate delta AIC
aic_vals <- c(
  ER  = AIC(s3.er),
  SYM = AIC(s3.sym),
  ARD = AIC(s3.ard))
delta_aic <- aic_vals - min(aic_vals)
delta_aic

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit
# SYM/ER are equally supported so plot both

####################################################### 

# Marginal ASR SYM
####################################################### 
# marginal ancestral reconstruction
s3_marginal_sym<-ancr(s3.sym, type="marginal")
s3_marginal_sym
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s3.sym)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s3_QMatrix_americas_SYM_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s3.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s3_marginal_sym,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topright",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s3.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s3_marginal_sym,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topright",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 


# Marginal ASR ER
####################################################### 
# marginal ancestral reconstruction
s3_marginal_er<-ancr(s3.er, type="marginal")
s3_marginal_er
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s3.er)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s3_QMatrix_americas_ER_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s3.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s3_marginal_er,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topright",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s3.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s3_marginal_er,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topright",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

########################
#
# END SEGMENT 3
#
########################



















#### Segment 4 ####

########################
# SIMPLE REGIONS 
# AMERICAS
# SEGMENT 4
########################

# plot base tree
#######################################################
# read in tree
s4.tree<-read.tree(file="./trees/BTV_Americas_s4_aln.fasta.treefile")
s4.tree

# does this tree have any zero-length branches?
any(s4.tree$edge.length == 0)
# yes it does

plot_tree_with_colored_zero_length_branches <- function(tree) {
  
  # how many tips and nodes does this tree have?
  number_tips  = Ntip(tree)
  number_nodes = tree$Nnode
  
  # get indexes of all the zero-length edges
  zero_length_edges <- which(tree$edge.length == 0) 
  
  # all the downstream nodes/tips of zero-length branches
  # all of them
  zero_length_children <- tree$edge[tree$edge.length == 0, 2]
  # just the tips
  zero_length_tips  <- zero_length_children[zero_length_children <= number_tips]
  # just the nodes
  zero_length_nodes <- zero_length_children[zero_length_children > number_tips]
  
  # create color schemes
  tip_colors  = rep("black", number_tips)
  edge_colors = rep("black", length(tree$edge.length))
  edge_widths = rep(0.5, length(tree$edge.length))
  
  # color any 0-length tips or branches red, and thick
  tip_colors[zero_length_tips] <- "red"
  edge_colors[zero_length_edges] <- "red"
  edge_widths[zero_length_edges] <- 2
  
  # make the plot
  plot(tree, 
       tip.color = tip_colors, 
       edge.color = edge_colors, 
       edge.width = edge_widths, 
       cex = 0.25)
}

# plot s4.tree with tips downstream of 0-length branches colored red
plot_tree_with_colored_zero_length_branches(s4.tree)

# now convert 0-length branches into branches of length 1e-6
# using base R indexing
s4.tree.no_zero <- s4.tree
s4.tree.no_zero$edge.length[s4.tree.no_zero$edge.length == 0] <- 1e-6

# does this have any 0-length branches?
any(s4.tree.no_zero$edge.length == 0 )

# plot it
plot_tree_with_colored_zero_length_branches(s4.tree.no_zero)

# midpoint root tree
s4.tree.midpoint_rooted <- midpoint_root(s4.tree)

# write out midpoint rooted trees in .nexus format
write.tree(s4.tree.midpoint_rooted, file="./trees/rooted/s4_tree_R_midpoint_rooted.newick")

# make sure that tip label name structure matches data
s4.tree.midpoint_rooted$tip.label <- gsub("/", "_", s4.tree.midpoint_rooted$tip.label)
s4.tree.midpoint_rooted$tip.label <- gsub("\\|", "_", s4.tree.midpoint_rooted$tip.label)

# read in metadata
s4.data<-read.csv(file="A2_Americas_Metadata.csv", row.names=1)

# keep only the tips that are found in the tree
s4.data.subset <- s4.data %>% filter(rownames(s4.data) %in% s4.tree.midpoint_rooted$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s4.tree.midpoint_rooted, s4.data.subset)
summary(chk)

# create a named vector
s4.trait <- setNames(s4.data.subset[,"trait"], rownames(s4.data.subset))
unique(s4.trait)

plotTree(s4.tree.midpoint_rooted,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s4.trait),names(s4.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s4.tree.midpoint_rooted),2,
             dimnames=list(s4.tree.midpoint_rooted$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s4.tree.midpoint_rooted)) 
  pies[s4.tree.midpoint_rooted$tip.label[i],
       X[[s4.tree.midpoint_rooted$tip.label[i]]]]<-
  rep(1/length(X[[s4.tree.midpoint_rooted$tip.label[i]]]),
      length(X[[s4.tree.midpoint_rooted$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s4.tree.midpoint_rooted$tip.label, tip = 1:Ntip(s4.tree.midpoint_rooted), adj = -0.05, 
          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(cols),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)

#######################################################
# fit model
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
s4.ard<-fitMk(s4.tree.midpoint_rooted,s4.trait,model="ARD")
s4.ard

### SYM
s4.sym<-fitMk(s4.tree.midpoint_rooted,s4.trait,model="SYM")
s4.sym

### ER
s4.er<-fitMk(s4.tree.midpoint_rooted,s4.trait,model="ER")
s4.er

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s4.ard),
                    logLik(s4.sym),
                    logLik(s4.er)),
           AIC=c(AIC(s4.ard),
                 AIC(s4.sym),
                 AIC(s4.er)))


# calculate delta AIC
aic_vals <- c(
  ER  = AIC(s4.er),
  SYM = AIC(s4.sym),
  ARD = AIC(s4.ard))
delta_aic <- aic_vals - min(aic_vals)
delta_aic

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit
# SYM/ER are equally supported so plot both

####################################################### 

# Marginal ASR SYM
####################################################### 
# marginal ancestral reconstruction
s4_marginal_sym<-ancr(s4.sym, type="marginal")
s4_marginal_sym
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s4.sym)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s4_QMatrix_americas_SYM_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s4.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s4_marginal_sym,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s4.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s4_marginal_sym,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 


# Marginal ASR ER
####################################################### 
# marginal ancestral reconstruction
s4_marginal_er<-ancr(s4.er, type="marginal")
s4_marginal_er
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s4.er)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s4_QMatrix_americas_ER_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s4.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s4_marginal_er,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s4.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s4_marginal_er,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

########################
#
# END SEGMENT 4
#
########################



















#### Segment 5 ####

########################
# SIMPLE REGIONS 
# AMERICAS
# SEGMENT 5
########################

# plot base tree
#######################################################
# read in tree
s5.tree<-read.tree(file="./trees/BTV_Americas_s5_aln_REPEAT.fasta.treefile")
s5.tree

# does this tree have any zero-length branches?
any(s5.tree$edge.length == 0)
# yes it does

plot_tree_with_colored_zero_length_branches <- function(tree) {
  
  # how many tips and nodes does this tree have?
  number_tips  = Ntip(tree)
  number_nodes = tree$Nnode
  
  # get indexes of all the zero-length edges
  zero_length_edges <- which(tree$edge.length == 0) 
  
  # all the downstream nodes/tips of zero-length branches
  # all of them
  zero_length_children <- tree$edge[tree$edge.length == 0, 2]
  # just the tips
  zero_length_tips  <- zero_length_children[zero_length_children <= number_tips]
  # just the nodes
  zero_length_nodes <- zero_length_children[zero_length_children > number_tips]
  
  # create color schemes
  tip_colors  = rep("black", number_tips)
  edge_colors = rep("black", length(tree$edge.length))
  edge_widths = rep(0.5, length(tree$edge.length))
  
  # color any 0-length tips or branches red, and thick
  tip_colors[zero_length_tips] <- "red"
  edge_colors[zero_length_edges] <- "red"
  edge_widths[zero_length_edges] <- 2
  
  # make the plot
  plot(tree, 
       tip.color = tip_colors, 
       edge.color = edge_colors, 
       edge.width = edge_widths, 
       cex = 0.25)
}

# plot s5.tree with tips downstream of 0-length branches colored red
plot_tree_with_colored_zero_length_branches(s5.tree)

# now convert 0-length branches into branches of length 1e-6
# using base R indexing
s5.tree.no_zero <- s5.tree
s5.tree.no_zero$edge.length[s5.tree.no_zero$edge.length == 0] <- 1e-6

# does this have any 0-length branches?
any(s5.tree.no_zero$edge.length == 0 )

# plot it
plot_tree_with_colored_zero_length_branches(s5.tree.no_zero)

# midpoint root tree
s5.tree.midpoint_rooted <- midpoint_root(s5.tree)

# write out midpoint rooted trees in .nexus format
write.tree(s5.tree.midpoint_rooted, file="./trees/rooted/s5_tree_R_midpoint_rooted.newick")

# make sure that tip label name structure matches data
s5.tree.midpoint_rooted$tip.label <- gsub("/", "_", s5.tree.midpoint_rooted$tip.label)
s5.tree.midpoint_rooted$tip.label <- gsub("\\|", "_", s5.tree.midpoint_rooted$tip.label)

# read in metadata
s5.data<-read.csv(file="A2_Americas_Metadata.csv", row.names=1)

# keep only the tips that are found in the tree
s5.data.subset <- s5.data %>% filter(rownames(s5.data) %in% s5.tree.midpoint_rooted$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s5.tree.midpoint_rooted, s5.data.subset)
summary(chk)

# create a named vector
s5.trait <- setNames(s5.data.subset[,"trait"], rownames(s5.data.subset))
unique(s5.trait)

plotTree(s5.tree.midpoint_rooted,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s5.trait),names(s5.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s5.tree.midpoint_rooted),2,
             dimnames=list(s5.tree.midpoint_rooted$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s5.tree.midpoint_rooted)) 
  pies[s5.tree.midpoint_rooted$tip.label[i],
       X[[s5.tree.midpoint_rooted$tip.label[i]]]]<-
  rep(1/length(X[[s5.tree.midpoint_rooted$tip.label[i]]]),
      length(X[[s5.tree.midpoint_rooted$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s5.tree.midpoint_rooted$tip.label, tip = 1:Ntip(s5.tree.midpoint_rooted), adj = -0.05, 
          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(cols),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)

#######################################################
# fit model
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
s5.ard<-fitMk(s5.tree.midpoint_rooted,s5.trait,model="ARD")
s5.ard

### SYM
s5.sym<-fitMk(s5.tree.midpoint_rooted,s5.trait,model="SYM")
s5.sym

### ER
s5.er<-fitMk(s5.tree.midpoint_rooted,s5.trait,model="ER")
s5.er

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s5.ard),
                    logLik(s5.sym),
                    logLik(s5.er)),
           AIC=c(AIC(s5.ard),
                 AIC(s5.sym),
                 AIC(s5.er)))


# calculate delta AIC
aic_vals <- c(
  ER  = AIC(s5.er),
  SYM = AIC(s5.sym),
  ARD = AIC(s5.ard))
delta_aic <- aic_vals - min(aic_vals)
delta_aic

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit
# SYM/ER are equally supported so plot both

####################################################### 

# Marginal ASR SYM
####################################################### 
# marginal ancestral reconstruction
s5_marginal_sym<-ancr(s5.sym, type="marginal")
s5_marginal_sym
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s5.sym)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s5_QMatrix_americas_SYM_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s5.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s5_marginal_sym,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s5.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s5_marginal_sym,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 


# Marginal ASR ER
####################################################### 
# marginal ancestral reconstruction
s5_marginal_er<-ancr(s5.er, type="marginal")
s5_marginal_er
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s5.er)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s5_QMatrix_americas_ER_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s5.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s5_marginal_er,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s5.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s5_marginal_er,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

########################
#
# END SEGMENT 5
#
########################



















#### Segment 6 ####

########################
# SIMPLE REGIONS 
# AMERICAS
# SEGMENT 6
########################

# plot base tree
#######################################################
# read in tree
s6.tree<-read.tree(file="./trees/BTV_Americas_s6_aln.fasta.treefile")
s6.tree

# does this tree have any zero-length branches?
any(s6.tree$edge.length == 0)
# yes it does

plot_tree_with_colored_zero_length_branches <- function(tree) {
  
  # how many tips and nodes does this tree have?
  number_tips  = Ntip(tree)
  number_nodes = tree$Nnode
  
  # get indexes of all the zero-length edges
  zero_length_edges <- which(tree$edge.length == 0) 
  
  # all the downstream nodes/tips of zero-length branches
  # all of them
  zero_length_children <- tree$edge[tree$edge.length == 0, 2]
  # just the tips
  zero_length_tips  <- zero_length_children[zero_length_children <= number_tips]
  # just the nodes
  zero_length_nodes <- zero_length_children[zero_length_children > number_tips]
  
  # create color schemes
  tip_colors  = rep("black", number_tips)
  edge_colors = rep("black", length(tree$edge.length))
  edge_widths = rep(0.5, length(tree$edge.length))
  
  # color any 0-length tips or branches red, and thick
  tip_colors[zero_length_tips] <- "red"
  edge_colors[zero_length_edges] <- "red"
  edge_widths[zero_length_edges] <- 2
  
  # make the plot
  plot(tree, 
       tip.color = tip_colors, 
       edge.color = edge_colors, 
       edge.width = edge_widths, 
       cex = 0.25)
}

# plot s6.tree with tips downstream of 0-length branches colored red
plot_tree_with_colored_zero_length_branches(s6.tree)

# now convert 0-length branches into branches of length 1e-6
# using base R indexing
s6.tree.no_zero <- s6.tree
s6.tree.no_zero$edge.length[s6.tree.no_zero$edge.length == 0] <- 1e-6

# does this have any 0-length branches?
any(s6.tree.no_zero$edge.length == 0 )

# plot it
plot_tree_with_colored_zero_length_branches(s6.tree.no_zero)

# midpoint root tree
s6.tree.midpoint_rooted <- midpoint_root(s6.tree)

# write out midpoint rooted trees in .nexus format
write.tree(s6.tree.midpoint_rooted, file="./trees/rooted/s6_tree_R_midpoint_rooted.newick")

# make sure that tip label name structure matches data
s6.tree.midpoint_rooted$tip.label <- gsub("/", "_", s6.tree.midpoint_rooted$tip.label)
s6.tree.midpoint_rooted$tip.label <- gsub("\\|", "_", s6.tree.midpoint_rooted$tip.label)

# read in metadata
s6.data<-read.csv(file="A2_Americas_Metadata.csv", row.names=1)

# keep only the tips that are found in the tree
s6.data.subset <- s6.data %>% filter(rownames(s6.data) %in% s6.tree.midpoint_rooted$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s6.tree.midpoint_rooted, s6.data.subset)
summary(chk)

# create a named vector
s6.trait <- setNames(s6.data.subset[,"trait"], rownames(s6.data.subset))
unique(s6.trait)

plotTree(s6.tree.midpoint_rooted,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s6.trait),names(s6.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s6.tree.midpoint_rooted),2,
             dimnames=list(s6.tree.midpoint_rooted$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s6.tree.midpoint_rooted)) 
  pies[s6.tree.midpoint_rooted$tip.label[i],
       X[[s6.tree.midpoint_rooted$tip.label[i]]]]<-
  rep(1/length(X[[s6.tree.midpoint_rooted$tip.label[i]]]),
      length(X[[s6.tree.midpoint_rooted$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s6.tree.midpoint_rooted$tip.label, tip = 1:Ntip(s6.tree.midpoint_rooted), adj = -0.05, 
          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(cols),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)

#######################################################
# fit model
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
s6.ard<-fitMk(s6.tree.midpoint_rooted,s6.trait,model="ARD")
s6.ard

### SYM
s6.sym<-fitMk(s6.tree.midpoint_rooted,s6.trait,model="SYM")
s6.sym

### ER
s6.er<-fitMk(s6.tree.midpoint_rooted,s6.trait,model="ER")
s6.er

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s6.ard),
                    logLik(s6.sym),
                    logLik(s6.er)),
           AIC=c(AIC(s6.ard),
                 AIC(s6.sym),
                 AIC(s6.er)))


# calculate delta AIC
aic_vals <- c(
  ER  = AIC(s6.er),
  SYM = AIC(s6.sym),
  ARD = AIC(s6.ard))
delta_aic <- aic_vals - min(aic_vals)
delta_aic

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit
# ARD fits best, only plot that 

####################################################### 

# Marginal ASR
####################################################### 
# marginal ancestral reconstruction
s6_marginal<-ancr(s6.ard, type="marginal")
s6_marginal
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s6.ard)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s6_QMatrix_americas_ARD_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s6.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s6_marginal,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s6.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s6_marginal,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

########################
#
# END SEGMENT 6
#
########################



















#### Segment 7 ####

########################
# SIMPLE REGIONS 
# AMERICAS
# SEGMENT 7
########################

# plot base tree
#######################################################
# read in tree
s7.tree<-read.tree(file="./trees/BTV_Americas_s7_aln.fasta.treefile")
s7.tree

# does this tree have any zero-length branches?
any(s7.tree$edge.length == 0)
# yes it does

plot_tree_with_colored_zero_length_branches <- function(tree) {
  
  # how many tips and nodes does this tree have?
  number_tips  = Ntip(tree)
  number_nodes = tree$Nnode
  
  # get indexes of all the zero-length edges
  zero_length_edges <- which(tree$edge.length == 0) 
  
  # all the downstream nodes/tips of zero-length branches
  # all of them
  zero_length_children <- tree$edge[tree$edge.length == 0, 2]
  # just the tips
  zero_length_tips  <- zero_length_children[zero_length_children <= number_tips]
  # just the nodes
  zero_length_nodes <- zero_length_children[zero_length_children > number_tips]
  
  # create color schemes
  tip_colors  = rep("black", number_tips)
  edge_colors = rep("black", length(tree$edge.length))
  edge_widths = rep(0.5, length(tree$edge.length))
  
  # color any 0-length tips or branches red, and thick
  tip_colors[zero_length_tips] <- "red"
  edge_colors[zero_length_edges] <- "red"
  edge_widths[zero_length_edges] <- 2
  
  # make the plot
  plot(tree, 
       tip.color = tip_colors, 
       edge.color = edge_colors, 
       edge.width = edge_widths, 
       cex = 0.25)
}

# plot s7.tree with tips downstream of 0-length branches colored red
plot_tree_with_colored_zero_length_branches(s7.tree)

# now convert 0-length branches into branches of length 1e-6
# using base R indexing
s7.tree.no_zero <- s7.tree
s7.tree.no_zero$edge.length[s7.tree.no_zero$edge.length == 0] <- 1e-6

# does this have any 0-length branches?
any(s7.tree.no_zero$edge.length == 0 )

# plot it
plot_tree_with_colored_zero_length_branches(s7.tree.no_zero)

# midpoint root tree
s7.tree.midpoint_rooted <- midpoint_root(s7.tree)

# write out midpoint rooted trees in .nexus format
write.tree(s7.tree.midpoint_rooted, file="./trees/rooted/s7_tree_R_midpoint_rooted.newick")

# make sure that tip label name structure matches data
s7.tree.midpoint_rooted$tip.label <- gsub("/", "_", s7.tree.midpoint_rooted$tip.label)
s7.tree.midpoint_rooted$tip.label <- gsub("\\|", "_", s7.tree.midpoint_rooted$tip.label)

# read in metadata
s7.data<-read.csv(file="A2_Americas_Metadata.csv", row.names=1)

# keep only the tips that are found in the tree
s7.data.subset <- s7.data %>% filter(rownames(s7.data) %in% s7.tree.midpoint_rooted$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s7.tree.midpoint_rooted, s7.data.subset)
summary(chk)

# create a named vector
s7.trait <- setNames(s7.data.subset[,"trait"], rownames(s7.data.subset))
unique(s7.trait)

plotTree(s7.tree.midpoint_rooted,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s7.trait),names(s7.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s7.tree.midpoint_rooted),2,
             dimnames=list(s7.tree.midpoint_rooted$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s7.tree.midpoint_rooted)) 
  pies[s7.tree.midpoint_rooted$tip.label[i],
       X[[s7.tree.midpoint_rooted$tip.label[i]]]]<-
  rep(1/length(X[[s7.tree.midpoint_rooted$tip.label[i]]]),
      length(X[[s7.tree.midpoint_rooted$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s7.tree.midpoint_rooted$tip.label, tip = 1:Ntip(s7.tree.midpoint_rooted), adj = -0.05, 
          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(cols),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)

#######################################################
# fit model
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
s7.ard<-fitMk(s7.tree.midpoint_rooted,s7.trait,model="ARD")
s7.ard

### SYM
s7.sym<-fitMk(s7.tree.midpoint_rooted,s7.trait,model="SYM")
s7.sym

### ER
s7.er<-fitMk(s7.tree.midpoint_rooted,s7.trait,model="ER")
s7.er

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s7.ard),
                    logLik(s7.sym),
                    logLik(s7.er)),
           AIC=c(AIC(s7.ard),
                 AIC(s7.sym),
                 AIC(s7.er)))


# calculate delta AIC
aic_vals <- c(
  ER  = AIC(s7.er),
  SYM = AIC(s7.sym),
  ARD = AIC(s7.ard))
delta_aic <- aic_vals - min(aic_vals)
delta_aic

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit
# ARD fits best but does not beat out ER/SYM. Delta AIC < 2. 
# Plot ER/SYM

####################################################### 

# Marginal ASR SYM
####################################################### 
# marginal ancestral reconstruction
s7_marginal_sym<-ancr(s7.sym, type="marginal")
s7_marginal_sym
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s7.sym)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s7_QMatrix_americas_SYM_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s7.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s7_marginal_sym,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s7.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s7_marginal_sym,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 


# Marginal ASR ER
####################################################### 
# marginal ancestral reconstruction
s7_marginal_er<-ancr(s7.er, type="marginal")
s7_marginal_er
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s7.er)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s7_QMatrix_americas_ER_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s7.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s7_marginal_er,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s7.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s7_marginal_er,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

########################
#
# END SEGMENT 7
#
########################



















#### Segment 8 ####

########################
# SIMPLE REGIONS 
# AMERICAS
# SEGMENT 8
########################

# plot base tree
#######################################################
# read in tree
s8.tree<-read.tree(file="./trees/BTV_Americas_s8_aln.fasta.treefile")
s8.tree

# does this tree have any zero-length branches?
any(s8.tree$edge.length == 0)
# yes it does

plot_tree_with_colored_zero_length_branches <- function(tree) {
  
  # how many tips and nodes does this tree have?
  number_tips  = Ntip(tree)
  number_nodes = tree$Nnode
  
  # get indexes of all the zero-length edges
  zero_length_edges <- which(tree$edge.length == 0) 
  
  # all the downstream nodes/tips of zero-length branches
  # all of them
  zero_length_children <- tree$edge[tree$edge.length == 0, 2]
  # just the tips
  zero_length_tips  <- zero_length_children[zero_length_children <= number_tips]
  # just the nodes
  zero_length_nodes <- zero_length_children[zero_length_children > number_tips]
  
  # create color schemes
  tip_colors  = rep("black", number_tips)
  edge_colors = rep("black", length(tree$edge.length))
  edge_widths = rep(0.5, length(tree$edge.length))
  
  # color any 0-length tips or branches red, and thick
  tip_colors[zero_length_tips] <- "red"
  edge_colors[zero_length_edges] <- "red"
  edge_widths[zero_length_edges] <- 2
  
  # make the plot
  plot(tree, 
       tip.color = tip_colors, 
       edge.color = edge_colors, 
       edge.width = edge_widths, 
       cex = 0.25)
}

# plot s8.tree with tips downstream of 0-length branches colored red
plot_tree_with_colored_zero_length_branches(s8.tree)

# now convert 0-length branches into branches of length 1e-6
# using base R indexing
s8.tree.no_zero <- s8.tree
s8.tree.no_zero$edge.length[s8.tree.no_zero$edge.length == 0] <- 1e-6

# does this have any 0-length branches?
any(s8.tree.no_zero$edge.length == 0 )

# plot it
plot_tree_with_colored_zero_length_branches(s8.tree.no_zero)

# midpoint root tree
s8.tree.midpoint_rooted <- midpoint_root(s8.tree)

# write out midpoint rooted trees in .nexus format
write.tree(s8.tree.midpoint_rooted, file="./trees/rooted/s8_tree_R_midpoint_rooted.newick")

# make sure that tip label name structure matches data
s8.tree.midpoint_rooted$tip.label <- gsub("/", "_", s8.tree.midpoint_rooted$tip.label)
s8.tree.midpoint_rooted$tip.label <- gsub("\\|", "_", s8.tree.midpoint_rooted$tip.label)

# read in metadata
s8.data<-read.csv(file="A2_Americas_Metadata.csv", row.names=1)

# keep only the tips that are found in the tree
s8.data.subset <- s8.data %>% filter(rownames(s8.data) %in% s8.tree.midpoint_rooted$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s8.tree.midpoint_rooted, s8.data.subset)
summary(chk)

# create a named vector
s8.trait <- setNames(s8.data.subset[,"trait"], rownames(s8.data.subset))
unique(s8.trait)

plotTree(s8.tree.midpoint_rooted,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s8.trait),names(s8.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s8.tree.midpoint_rooted),2,
             dimnames=list(s8.tree.midpoint_rooted$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s8.tree.midpoint_rooted)) 
  pies[s8.tree.midpoint_rooted$tip.label[i],
       X[[s8.tree.midpoint_rooted$tip.label[i]]]]<-
  rep(1/length(X[[s8.tree.midpoint_rooted$tip.label[i]]]),
      length(X[[s8.tree.midpoint_rooted$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s8.tree.midpoint_rooted$tip.label, tip = 1:Ntip(s8.tree.midpoint_rooted), adj = -0.05, 
          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(cols),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)

#######################################################
# fit model
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
s8.ard<-fitMk(s8.tree.midpoint_rooted,s8.trait,model="ARD")
s8.ard

### SYM
s8.sym<-fitMk(s8.tree.midpoint_rooted,s8.trait,model="SYM")
s8.sym

### ER
s8.er<-fitMk(s8.tree.midpoint_rooted,s8.trait,model="ER")
s8.er

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s8.ard),
                    logLik(s8.sym),
                    logLik(s8.er)),
           AIC=c(AIC(s8.ard),
                 AIC(s8.sym),
                 AIC(s8.er)))


# calculate delta AIC
aic_vals <- c(
  ER  = AIC(s8.er),
  SYM = AIC(s8.sym),
  ARD = AIC(s8.ard))
delta_aic <- aic_vals - min(aic_vals)
delta_aic

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit
# ARD fits best but does not beat out ER/SYM. Delta AIC < 2. 
# Plot ER/SYM

####################################################### 

# Marginal ASR SYM
####################################################### 
# marginal ancestral reconstruction
s8_marginal_sym<-ancr(s8.sym, type="marginal")
s8_marginal_sym
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s8.sym)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s8_QMatrix_americas_SYM_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s8.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s8_marginal_sym,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s8.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s8_marginal_sym,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 


# Marginal ASR ER
####################################################### 
# marginal ancestral reconstruction
s8_marginal_er<-ancr(s8.er, type="marginal")
s8_marginal_er
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s8.er)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s8_QMatrix_americas_ER_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s8.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s8_marginal_er,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s8.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s8_marginal_er,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 


########################
#
# END SEGMENT 8
#
########################



















#### Segment 9 ####

########################
# SIMPLE REGIONS 
# AMERICAS
# SEGMENT 9
########################

# plot base tree
#######################################################
# read in tree
s9.tree<-read.tree(file="./trees/BTV_Americas_s9_aln.fasta.treefile")
s9.tree

# does this tree have any zero-length branches?
any(s9.tree$edge.length == 0)
# yes it does

plot_tree_with_colored_zero_length_branches <- function(tree) {
  
  # how many tips and nodes does this tree have?
  number_tips  = Ntip(tree)
  number_nodes = tree$Nnode
  
  # get indexes of all the zero-length edges
  zero_length_edges <- which(tree$edge.length == 0) 
  
  # all the downstream nodes/tips of zero-length branches
  # all of them
  zero_length_children <- tree$edge[tree$edge.length == 0, 2]
  # just the tips
  zero_length_tips  <- zero_length_children[zero_length_children <= number_tips]
  # just the nodes
  zero_length_nodes <- zero_length_children[zero_length_children > number_tips]
  
  # create color schemes
  tip_colors  = rep("black", number_tips)
  edge_colors = rep("black", length(tree$edge.length))
  edge_widths = rep(0.5, length(tree$edge.length))
  
  # color any 0-length tips or branches red, and thick
  tip_colors[zero_length_tips] <- "red"
  edge_colors[zero_length_edges] <- "red"
  edge_widths[zero_length_edges] <- 2
  
  # make the plot
  plot(tree, 
       tip.color = tip_colors, 
       edge.color = edge_colors, 
       edge.width = edge_widths, 
       cex = 0.25)
}

# plot s9.tree with tips downstream of 0-length branches colored red
plot_tree_with_colored_zero_length_branches(s9.tree)

# now convert 0-length branches into branches of length 1e-6
# using base R indexing
s9.tree.no_zero <- s9.tree
s9.tree.no_zero$edge.length[s9.tree.no_zero$edge.length == 0] <- 1e-6

# does this have any 0-length branches?
any(s9.tree.no_zero$edge.length == 0 )

# plot it
plot_tree_with_colored_zero_length_branches(s9.tree.no_zero)

# midpoint root tree
s9.tree.midpoint_rooted <- midpoint_root(s9.tree)

# write out midpoint rooted trees in .nexus format
write.tree(s9.tree.midpoint_rooted, file="./trees/rooted/s9_tree_R_midpoint_rooted.newick")

# make sure that tip label name structure matches data
s9.tree.midpoint_rooted$tip.label <- gsub("/", "_", s9.tree.midpoint_rooted$tip.label)
s9.tree.midpoint_rooted$tip.label <- gsub("\\|", "_", s9.tree.midpoint_rooted$tip.label)

# read in metadata
s9.data<-read.csv(file="A2_Americas_Metadata.csv", row.names=1)

# keep only the tips that are found in the tree
s9.data.subset <- s9.data %>% filter(rownames(s9.data) %in% s9.tree.midpoint_rooted$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s9.tree.midpoint_rooted, s9.data.subset)
summary(chk)

# create a named vector
s9.trait <- setNames(s9.data.subset[,"trait"], rownames(s9.data.subset))
unique(s9.trait)

plotTree(s9.tree.midpoint_rooted,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s9.trait),names(s9.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s9.tree.midpoint_rooted),2,
             dimnames=list(s9.tree.midpoint_rooted$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s9.tree.midpoint_rooted)) 
  pies[s9.tree.midpoint_rooted$tip.label[i],
       X[[s9.tree.midpoint_rooted$tip.label[i]]]]<-
  rep(1/length(X[[s9.tree.midpoint_rooted$tip.label[i]]]),
      length(X[[s9.tree.midpoint_rooted$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s9.tree.midpoint_rooted$tip.label, tip = 1:Ntip(s9.tree.midpoint_rooted), adj = -0.05, 
          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(cols),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)

#######################################################
# fit model
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
s9.ard<-fitMk(s9.tree.midpoint_rooted,s9.trait,model="ARD")
s9.ard

### SYM
s9.sym<-fitMk(s9.tree.midpoint_rooted,s9.trait,model="SYM")
s9.sym

### ER
s9.er<-fitMk(s9.tree.midpoint_rooted,s9.trait,model="ER")
s9.er

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s9.ard),
                    logLik(s9.sym),
                    logLik(s9.er)),
           AIC=c(AIC(s9.ard),
                 AIC(s9.sym),
                 AIC(s9.er)))


# calculate delta AIC
aic_vals <- c(
  ER  = AIC(s9.er),
  SYM = AIC(s9.sym),
  ARD = AIC(s9.ard))
delta_aic <- aic_vals - min(aic_vals)
delta_aic

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit
# SYM/ER are equally supported so plot both

####################################################### 

# Marginal ASR SYM
####################################################### 
# marginal ancestral reconstruction
s9_marginal_sym<-ancr(s9.sym, type="marginal")
s9_marginal_sym
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s9.sym)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s9_QMatrix_americas_SYM_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s9.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s9_marginal_sym,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s9.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s9_marginal_sym,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 


# Marginal ASR ER
####################################################### 
# marginal ancestral reconstruction
s9_marginal_er<-ancr(s9.er, type="marginal")
s9_marginal_er
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s9.er)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s9_QMatrix_americas_ER_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s9.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s9_marginal_er,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s9.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s9_marginal_er,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

########################
#
# END SEGMENT 9
#
########################



















#### Segment 10 ####

########################
# SIMPLE REGIONS 
# AMERICAS
# SEGMENT 10
########################

# plot base tree
#######################################################
# read in tree
s10.tree<-read.tree(file="./trees/BTV_Americas_s10_aln.fasta.treefile")
s10.tree

# does this tree have any zero-length branches?
any(s10.tree$edge.length == 0)
# yes it does

plot_tree_with_colored_zero_length_branches <- function(tree) {
  
  # how many tips and nodes does this tree have?
  number_tips  = Ntip(tree)
  number_nodes = tree$Nnode
  
  # get indexes of all the zero-length edges
  zero_length_edges <- which(tree$edge.length == 0) 
  
  # all the downstream nodes/tips of zero-length branches
  # all of them
  zero_length_children <- tree$edge[tree$edge.length == 0, 2]
  # just the tips
  zero_length_tips  <- zero_length_children[zero_length_children <= number_tips]
  # just the nodes
  zero_length_nodes <- zero_length_children[zero_length_children > number_tips]
  
  # create color schemes
  tip_colors  = rep("black", number_tips)
  edge_colors = rep("black", length(tree$edge.length))
  edge_widths = rep(0.5, length(tree$edge.length))
  
  # color any 0-length tips or branches red, and thick
  tip_colors[zero_length_tips] <- "red"
  edge_colors[zero_length_edges] <- "red"
  edge_widths[zero_length_edges] <- 2
  
  # make the plot
  plot(tree, 
       tip.color = tip_colors, 
       edge.color = edge_colors, 
       edge.width = edge_widths, 
       cex = 0.25)
}

# plot s10.tree with tips downstream of 0-length branches colored red
plot_tree_with_colored_zero_length_branches(s10.tree)

# now convert 0-length branches into branches of length 1e-6
# using base R indexing
s10.tree.no_zero <- s10.tree
s10.tree.no_zero$edge.length[s10.tree.no_zero$edge.length == 0] <- 1e-6

# does this have any 0-length branches?
any(s10.tree.no_zero$edge.length == 0 )

# plot it
plot_tree_with_colored_zero_length_branches(s10.tree.no_zero)

# midpoint root tree
s10.tree.midpoint_rooted <- midpoint_root(s10.tree)

# write out midpoint rooted trees in .nexus format
write.tree(s10.tree.midpoint_rooted, file="./trees/rooted/s10_tree_R_midpoint_rooted.newick")

# make sure that tip label name structure matches data
s10.tree.midpoint_rooted$tip.label <- gsub("/", "_", s10.tree.midpoint_rooted$tip.label)
s10.tree.midpoint_rooted$tip.label <- gsub("\\|", "_", s10.tree.midpoint_rooted$tip.label)

# read in metadata
s10.data<-read.csv(file="A2_Americas_Metadata.csv", row.names=1)

# keep only the tips that are found in the tree
s10.data.subset <- s10.data %>% filter(rownames(s10.data) %in% s10.tree.midpoint_rooted$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s10.tree.midpoint_rooted, s10.data.subset)
summary(chk)

# create a named vector
s10.trait <- setNames(s10.data.subset[,"trait"], rownames(s10.data.subset))
unique(s10.trait)

plotTree(s10.tree.midpoint_rooted,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s10.trait),names(s10.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s10.tree.midpoint_rooted),2,
             dimnames=list(s10.tree.midpoint_rooted$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s10.tree.midpoint_rooted)) 
  pies[s10.tree.midpoint_rooted$tip.label[i],
       X[[s10.tree.midpoint_rooted$tip.label[i]]]]<-
  rep(1/length(X[[s10.tree.midpoint_rooted$tip.label[i]]]),
      length(X[[s10.tree.midpoint_rooted$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s10.tree.midpoint_rooted$tip.label, tip = 1:Ntip(s10.tree.midpoint_rooted), adj = -0.05, 
          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(cols),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)

#######################################################
# fit model
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
s10.ard<-fitMk(s10.tree.midpoint_rooted,s10.trait,model="ARD")
s10.ard

### SYM
s10.sym<-fitMk(s10.tree.midpoint_rooted,s10.trait,model="SYM")
s10.sym

### ER
s10.er<-fitMk(s10.tree.midpoint_rooted,s10.trait,model="ER")
s10.er

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s10.ard),
                    logLik(s10.sym),
                    logLik(s10.er)),
           AIC=c(AIC(s10.ard),
                 AIC(s10.sym),
                 AIC(s10.er)))


# calculate delta AIC
aic_vals <- c(
  ER  = AIC(s10.er),
  SYM = AIC(s10.sym),
  ARD = AIC(s10.ard))
delta_aic <- aic_vals - min(aic_vals)
delta_aic

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit
# ARD fits best but does not beat out ER/SYM. Delta AIC < 2. 
# Plot ER/SYM

####################################################### 

# Marginal ASR SYM
####################################################### 
# marginal ancestral reconstruction
s10_marginal_sym<-ancr(s10.sym, type="marginal")
s10_marginal_sym
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s10.sym)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s10_QMatrix_americas_SYM_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s10.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s10_marginal_sym,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s10.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s10_marginal_sym,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 


# Marginal ASR ER
####################################################### 
# marginal ancestral reconstruction
s10_marginal_er<-ancr(s10.er, type="marginal")
s10_marginal_er
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s10.er)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s10_QMatrix_americas_ER_2Regions_103026.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s10.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s10_marginal_er,
     args.plotTree=list(type="phylogram", fsize=0.3),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

# Plot Marginal ASR without sample names
####################################################### 
# plot marginal only
state_order <- sort(unique(s10.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s10_marginal_er,
     args.plotTree=list(type="phylogram", ftype="off"),
     args.nodelabels=list(piecol=cols,cex=0.15),
     args.tiplabels=list(piecol=cols, cex=0.15),
     legend = FALSE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
add.scale.bar(x = lastPP$x.lim[1] + 0.05 * diff(lastPP$x.lim), 
              y = lastPP$y.lim[1] - 0.02 * diff(lastPP$y.lim),
              length = 0.05, lwd = 2, cex = 0.8)
my_legend_order <- c("Outside USA", "USA")
legend_colors <- cols[my_legend_order]
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(legend_colors),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)
# check to make sure colors are correct
cols
####################################################### 

########################
#
# END SEGMENT 10
#
########################