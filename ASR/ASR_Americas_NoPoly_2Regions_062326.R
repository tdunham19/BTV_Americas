#######################################################
# Aim 2 - Ancestral State Reconstruction (ASR)
# Using correct cluster files
# Americas Region Categories: Inside USA & Outside USA
# n = 163 (removed samples causing polytomies)
# TD 01.08.26
# TD 06.23.26 Repeated with trees made with bootstrapping
#######################################################

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
s1.tree<-read.nexus(file="BTV_Americas_s1_no_poly_aln_midpoint_rooted_nexus_REPEAT.tre")
s1.tree

# read in metadata
s1.data<-read.csv(file="A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s1.data.subset <- s1.data %>%
  filter(rownames(s1.data) %in% s1.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s1.tree, s1.data.subset)
summary(chk)

# create a named vector
s1.trait <- setNames(s1.data.subset[,"trait"], rownames(s1.data.subset))
unique(s1.trait)

plotTree(s1.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s1.trait),names(s1.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s1.tree),2,
             dimnames=list(s1.tree$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s1.tree)) 
  pies[s1.tree$tip.label[i],
       X[[s1.tree$tip.label[i]]]]<-
  rep(1/length(X[[s1.tree$tip.label[i]]]),
      length(X[[s1.tree$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s1.tree$tip.label, tip = 1:Ntip(s1.tree), adj = -0.05, 
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
# Since there no polymorphic states present use polyMk
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
# s1.ard<-fitMk(s1.tree,s1.trait,model="ARD")
# s1.ard
# save the object
# saveRDS(s1.ard, file = "s1_americas_ARD_NoPoly_2Regions_062526.rds")
# restore the object
s1.ard <- readRDS("s1_americas_ARD_NoPoly_2Regions_062526.rds")

# fit other models

### SYM
# s1.sym<-fitMk(s1.tree,s1.trait,model="SYM")
# s1.sym
# save the object
# saveRDS(s1.sym, file = "s1_americas_SYM_NoPoly_2Regions_062526.rds")
# restore the object
s1.sym <- readRDS("s1_americas_SYM_NoPoly_2Regions_062526.rds")

### ER
# s1.er<-fitMk(s1.tree,s1.trait,model="ER")
# s1.er
# save the object
# saveRDS(s1.er, file = "s1_americas_ER_NoPoly_2Regions_062526.rds")
# restore the object
s1.er <- readRDS("s1_americas_ER_NoPoly_2Regions_062526.rds")

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s1.ard),
                    logLik(s1.sym),
                    logLik(s1.er)),
           AIC=c(AIC(s1.ard),
                 AIC(s1.sym),
                 AIC(s1.er)))

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit

# transition_model    logLik      AIC
# 1              ARD -38.68624 81.37248 ****
# 2              SYM -39.97289 81.94577
# 3               ER -39.97289 81.94577

# ARD fits the situation the best and is the best model fit!

# plot the models
# plot(s1.ard)
# plot(s1.sym)
# plot(s1.er)

s1_aov<-anova(s1.ard,s1.sym,s1.er)
summary(s1_aov)

####################################################### 

# Marginal ASR on ARD model
####################################################### 
# marginal ancestral reconstruction
s1_marginal<-ancr(s1.ard, type="marginal")
s1_marginal
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s1.ard)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s1_QMatrix_americas_ARD_NoPoly_2Regions_062526.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s1.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s1_marginal,
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
plot(s1_marginal,
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
s2.tree<-read.nexus(file="BTV_Americas_s2_no_poly_aln_midpoint_rooted_nexus.tre")
s2.tree

# read in metadata
s2.data<-read.csv(file="A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s2.data.subset <- s2.data %>%
  filter(rownames(s2.data) %in% s2.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s2.tree, s2.data.subset)
summary(chk)

# create a named vector
s2.trait <- setNames(s2.data.subset[,"trait"], rownames(s2.data.subset))
unique(s2.trait)

plotTree(s2.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s2.trait),names(s2.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s2.tree),2,
             dimnames=list(s2.tree$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s2.tree)) 
  pies[s2.tree$tip.label[i],
       X[[s2.tree$tip.label[i]]]]<-
  rep(1/length(X[[s2.tree$tip.label[i]]]),
      length(X[[s2.tree$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s2.tree$tip.label, tip = 1:Ntip(s2.tree), adj = -0.05, 
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
# Since there no polymorphic states present use polyMk
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
# s2.ard<-fitMk(s2.tree,s2.trait,model="ARD")
# s2.ard
# save the object
# saveRDS(s2.ard, file = "s2_americas_ARD_NoPoly_2Regions_062326.rds")
# restore the object
s2.ard <- readRDS("s2_americas_ARD_NoPoly_2Regions_062326.rds")

# fit other models

### SYM
# s2.sym<-fitMk(s2.tree,s2.trait,model="SYM")
# s2.sym
# save the object
# saveRDS(s2.sym, file = "s2_americas_SYM_NoPoly_2Regions_062326.rds")
# restore the object
s2.sym <- readRDS("s2_americas_SYM_NoPoly_2Regions_062326.rds")

### ER
# s2.er<-fitMk(s2.tree,s2.trait,model="ER")
# s2.er
# save the object
# saveRDS(s2.er, file = "s2_americas_ER_NoPoly_2Regions_062326.rds")
# restore the object
s2.er <- readRDS("s2_americas_ER_NoPoly_2Regions_062326.rds")

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s2.ard),
                    logLik(s2.sym),
                    logLik(s2.er)),
           AIC=c(AIC(s2.ard),
                 AIC(s2.sym),
                 AIC(s2.er)))

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit

# transition_model    logLik      AIC
# 1              ARD -64.09777 132.1955 ******
# 2              SYM -66.66393 135.3279
# 3               ER -66.66393 135.3279

# ARD fits the situation the best and is the best model fit!

# plot the models
# plot(s2.ard)
# plot(s2.sym)
# plot(s2.er)

####################################################### 

# Marginal ASR on ARD model
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
# write.csv(Q, file = "s2_QMatrix_americas_ARD_NoPoly_2Regions_062326.csv", row.names = TRUE)
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
legend("bottomright",
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
s3.tree<-read.nexus(file="BTV_Americas_s3_no_poly_aln_midpoint_rooted_nexus_REPEAT.tre")
s3.tree

# read in metadata
s3.data<-read.csv(file="A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s3.data.subset <- s3.data %>%
  filter(rownames(s3.data) %in% s3.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s3.tree, s3.data.subset)
summary(chk)

# create a named vector
s3.trait <- setNames(s3.data.subset[,"trait"], rownames(s3.data.subset))
unique(s3.trait)

plotTree(s3.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s3.trait),names(s3.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s3.tree),2,
             dimnames=list(s3.tree$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s3.tree)) 
  pies[s3.tree$tip.label[i],
       X[[s3.tree$tip.label[i]]]]<-
  rep(1/length(X[[s3.tree$tip.label[i]]]),
      length(X[[s3.tree$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s3.tree$tip.label, tip = 1:Ntip(s3.tree), adj = -0.05, 
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
# Since there no polymorphic states present use polyMk
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
# s3.ard<-fitMk(s3.tree,s3.trait,model="ARD")
# s3.ard
# save the object
# saveRDS(s3.ard, file = "s3_americas_ARD_NoPoly_2Regions_062526.rds")
# restore the object
s3.ard <- readRDS("s3_americas_ARD_NoPoly_2Regions_062526.rds")

# fit other models

### SYM
# s3.sym<-fitMk(s3.tree,s3.trait,model="SYM")
# s3.sym
# save the object
# saveRDS(s3.sym, file = "s3_americas_SYM_NoPoly_2Regions_062526.rds")
# restore the object
s3.sym <- readRDS("s3_americas_SYM_NoPoly_2Regions_062526.rds")

### ER
# s3.er<-fitMk(s3.tree,s3.trait,model="ER")
# s3.er
# save the object
# saveRDS(s3.er, file = "s3_americas_ER_NoPoly_2Regions_062526.rds")
# restore the object
s3.er <- readRDS("s3_americas_ER_NoPoly_2Regions_062526.rds")

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s3.ard),
                    logLik(s3.sym),
                    logLik(s3.er)),
           AIC=c(AIC(s3.ard),
                 AIC(s3.sym),
                 AIC(s3.er)))

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit

# transition_model    logLik      AIC
# 1              ARD -44.75846 93.51691
# 2              SYM -44.75846 91.51691 ***
# 3               ER -44.75846 91.51691 ***

# ARD fits the situation the best ????

# plot the models
# plot(s3.ard)
# plot(s3.sym)
# plot(s3.er)

####################################################### 

#### ******* SYM AND ER FIT BETTER ******* ####

# Marginal ASR on ARD model
# Still using ARD model even though it was not the best fit!!!
####################################################### 
# marginal ancestral reconstruction
s3_marginal<-ancr(s3.ard, type="marginal")
s3_marginal
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s3.ard)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s3_QMatrix_americas_ARD_NoPoly_2Regions_062526.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s3.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s3_marginal,
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
legend("bottomright",
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
plot(s3_marginal,
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
legend("bottomright",
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
s4.tree<-read.nexus(file="BTV_Americas_s4_no_poly_aln_midpoint_rooted_nexus.tre")
s4.tree

# read in metadata
s4.data<-read.csv(file="A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s4.data.subset <- s4.data %>%
  filter(rownames(s4.data) %in% s4.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s4.tree, s4.data.subset)
summary(chk)

# create a named vector
s4.trait <- setNames(s4.data.subset[,"trait"], rownames(s4.data.subset))
unique(s4.trait)

plotTree(s4.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s4.trait),names(s4.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s4.tree),2,
             dimnames=list(s4.tree$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s4.tree)) 
  pies[s4.tree$tip.label[i],
       X[[s4.tree$tip.label[i]]]]<-
  rep(1/length(X[[s4.tree$tip.label[i]]]),
      length(X[[s4.tree$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s4.tree$tip.label, tip = 1:Ntip(s4.tree), adj = -0.05, 
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
# Since there no polymorphic states present use polyMk
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
# s4.ard<-fitMk(s4.tree,s4.trait,model="ARD")
# s4.ard
# save the object
# saveRDS(s4.ard, file = "s4_americas_ARD_NoPoly_2Regions_062326.rds")
# restore the object
s4.ard <- readRDS("s4_americas_ARD_NoPoly_2Regions_062326.rds")

# fit other models

### SYM
# s4.sym<-fitMk(s4.tree,s4.trait,model="SYM")
# s4.sym
# save the object
# saveRDS(s4.sym, file = "s4_americas_SYM_NoPoly_2Regions_062326.rds")
# restore the object
s4.sym <- readRDS("s4_americas_SYM_NoPoly_2Regions_062326.rds")

### ER
# s4.er<-fitMk(s4.tree,s4.trait,model="ER")
# s4.er
# save the object
# saveRDS(s4.er, file = "s4_americas_ER_NoPoly_2Regions_062326.rds")
# restore the object
s4.er <- readRDS("s4_americas_ER_NoPoly_2Regions_062326.rds")

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s4.ard),
                    logLik(s4.sym),
                    logLik(s4.er)),
           AIC=c(AIC(s4.ard),
                 AIC(s4.sym),
                 AIC(s4.er)))

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit

# transition_model    logLik      AIC
# 1              ARD -36.42646 76.85292
# 2              SYM -36.52801 75.05602 ***
# 3               ER -36.52801 75.05602 ***

# ARD fits the situation the best ????

# plot the models
# plot(s4.ard)
# plot(s4.sym)
# plot(s4.er)

####################################################### 

#### ******* SYM AND ER FIT BETTER ******* ####

# Marginal ASR on ARD model
# Still using ARD model even though it was not the best fit!!!
####################################################### 
# marginal ancestral reconstruction
s4_marginal<-ancr(s4.ard, type="marginal")
s4_marginal
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s4.ard)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s4_QMatrix_americas_ARD_NoPoly_2Regions_062326.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s4.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s4_marginal,
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
plot(s4_marginal,
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
s5.tree<-read.nexus(file="BTV_Americas_s5_no_poly_aln_midpoint_rooted_nexus_REPEAT.tre")
s5.tree

# read in metadata
s5.data<-read.csv(file="A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s5.data.subset <- s5.data %>%
  filter(rownames(s5.data) %in% s5.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s5.tree, s5.data.subset)
summary(chk)

# create a named vector
s5.trait <- setNames(s5.data.subset[,"trait"], rownames(s5.data.subset))
unique(s5.trait)

plotTree(s5.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s5.trait),names(s5.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s5.tree),2,
             dimnames=list(s5.tree$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s5.tree)) 
  pies[s5.tree$tip.label[i],
       X[[s5.tree$tip.label[i]]]]<-
  rep(1/length(X[[s5.tree$tip.label[i]]]),
      length(X[[s5.tree$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s5.tree$tip.label, tip = 1:Ntip(s5.tree), adj = -0.05, 
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
# Since there no polymorphic states present use polyMk
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
# s5.ard<-fitMk(s5.tree,s5.trait,model="ARD")
# s5.ard
# save the object
# saveRDS(s5.ard, file = "s5_americas_ARD_NoPoly_2Regions_062526.rds")
# restore the object
s5.ard <- readRDS("s5_americas_ARD_NoPoly_2Regions_062526.rds")

# fit other models

### SYM
# s5.sym<-fitMk(s5.tree,s5.trait,model="SYM")
# s5.sym
# save the object
# saveRDS(s5.sym, file = "s5_americas_SYM_NoPoly_2Regions_062526.rds")
# restore the object
s5.sym <- readRDS("s5_americas_SYM_NoPoly_2Regions_062526.rds")

### ER
# s5.er<-fitMk(s5.tree,s5.trait,model="ER")
# s5.er
# save the object
# saveRDS(s5.er, file = "s5_americas_ER_NoPoly_2Regions_062526.rds")
# restore the object
s5.er <- readRDS("s5_americas_ER_NoPoly_2Regions_062526.rds")

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s5.ard),
                    logLik(s5.sym),
                    logLik(s5.er)),
           AIC=c(AIC(s5.ard),
                 AIC(s5.sym),
                 AIC(s5.er)))

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit

# transition_model    logLik      AIC
# 1              ARD -27.83936 59.67872 *********
# 2              SYM -30.17326 62.34652
# 3               ER -30.17326 62.34652

# ARD fits the situation the best and the data the best!!

# plot the models
# plot(s5.ard)
# plot(s5.sym)
# plot(s5.er)

####################################################### 


# Marginal ASR on ARD model
####################################################### 
# marginal ancestral reconstruction
s5_marginal<-ancr(s5.ard, type="marginal")
s5_marginal
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s5.ard)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s5_QMatrix_americas_ARD_NoPoly_2Regions_062526.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s5.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s5_marginal,
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
plot(s5_marginal,
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
s6.tree<-read.nexus(file="BTV_Americas_s6_no_poly_aln_midpoint_rooted_nexus.tre")
s6.tree

# read in metadata
s6.data<-read.csv(file="A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s6.data.subset <- s6.data %>%
  filter(rownames(s6.data) %in% s6.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s6.tree, s6.data.subset)
summary(chk)

# create a named vector
s6.trait <- setNames(s6.data.subset[,"trait"], rownames(s6.data.subset))
unique(s6.trait)

plotTree(s6.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s6.trait),names(s6.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s6.tree),2,
             dimnames=list(s6.tree$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s6.tree)) 
  pies[s6.tree$tip.label[i],
       X[[s6.tree$tip.label[i]]]]<-
  rep(1/length(X[[s6.tree$tip.label[i]]]),
      length(X[[s6.tree$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s6.tree$tip.label, tip = 1:Ntip(s6.tree), adj = -0.05, 
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
# Since there no polymorphic states present use polyMk
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
# s6.ard<-fitMk(s6.tree,s6.trait,model="ARD")
# s6.ard
# save the object
# saveRDS(s6.ard, file = "s6_americas_ARD_NoPoly_2Regions_062326.rds")
# restore the object
s6.ard <- readRDS("s6_americas_ARD_NoPoly_2Regions_062326.rds")

# fit other models

### SYM
# s6.sym<-fitMk(s6.tree,s6.trait,model="SYM")
# s6.sym
# save the object
# saveRDS(s6.sym, file = "s6_americas_SYM_NoPoly_2Regions_062326.rds")
# restore the object
s6.sym <- readRDS("s6_americas_SYM_NoPoly_2Regions_062326.rds")

### ER
# s6.er<-fitMk(s6.tree,s6.trait,model="ER")
# s6.er
# save the object
# saveRDS(s6.er, file = "s6_americas_ER_NoPoly_2Regions_062326.rds")
# restore the object
s6.er <- readRDS("s6_americas_ER_NoPoly_2Regions_062326.rds")

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s6.ard),
                    logLik(s6.sym),
                    logLik(s6.er)),
           AIC=c(AIC(s6.ard),
                 AIC(s6.sym),
                 AIC(s6.er)))

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit

# transition_model    logLik      AIC
# 1              ARD -55.81724 115.6345 ********
# 2              SYM -59.01570 120.0314
# 3               ER -59.01570 120.0314

# ARD fits the situation the best and the data the best!!

# plot the models
# plot(s6.ard)
# plot(s6.sym)
# plot(s6.er)

####################################################### 


# Marginal ASR on ARD model
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
# write.csv(Q, file = "s6_QMatrix_americas_ARD_NoPoly_2Regions_062326.csv", row.names = TRUE)
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
s7.tree<-read.nexus(file="BTV_Americas_s7_no_poly_aln_midpoint_rooted_nexus.tre")
s7.tree

# read in metadata
s7.data<-read.csv(file="A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s7.data.subset <- s7.data %>%
  filter(rownames(s7.data) %in% s7.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s7.tree, s7.data.subset)
summary(chk)

# create a named vector
s7.trait <- setNames(s7.data.subset[,"trait"], rownames(s7.data.subset))
unique(s7.trait)

plotTree(s7.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s7.trait),names(s7.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s7.tree),2,
             dimnames=list(s7.tree$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s7.tree)) 
  pies[s7.tree$tip.label[i],
       X[[s7.tree$tip.label[i]]]]<-
  rep(1/length(X[[s7.tree$tip.label[i]]]),
      length(X[[s7.tree$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s7.tree$tip.label, tip = 1:Ntip(s7.tree), adj = -0.05, 
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
# Since there no polymorphic states present use polyMk
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
# s7.ard<-fitMk(s7.tree,s7.trait,model="ARD")
# s7.ard
# save the object
# saveRDS(s7.ard, file = "s7_americas_ARD_NoPoly_2Regions_062326.rds")
# restore the object
s7.ard <- readRDS("s7_americas_ARD_NoPoly_2Regions_062326.rds")

# fit other models

### SYM
# s7.sym<-fitMk(s7.tree,s7.trait,model="SYM")
# s7.sym
# save the object
# saveRDS(s7.sym, file = "s7_americas_SYM_NoPoly_2Regions_062326.rds")
# restore the object
s7.sym <- readRDS("s7_americas_SYM_NoPoly_2Regions_062326.rds")

### ER
# s7.er<-fitMk(s7.tree,s7.trait,model="ER")
# s7.er
# save the object
# saveRDS(s7.er, file = "s7_americas_ER_NoPoly_2Regions_062326.rds")
# restore the object
s7.er <- readRDS("s7_americas_ER_NoPoly_2Regions_062326.rds")

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s7.ard),
                    logLik(s7.sym),
                    logLik(s7.er)),
           AIC=c(AIC(s7.ard),
                 AIC(s7.sym),
                 AIC(s7.er)))

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit

# transition_model    logLik      AIC
# 1              ARD -46.54003 97.08007 *********
# 2              SYM -48.47054 98.94109
# 3               ER -48.47054 98.94109

# ARD fits the situation the best and the data the best!!

# plot the models
# plot(s7.ard)
# plot(s7.sym)
# plot(s7.er)

####################################################### 

# Marginal ASR on ARD model
####################################################### 
# marginal ancestral reconstruction
s7_marginal<-ancr(s7.ard, type="marginal")
s7_marginal
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s7.ard)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s7_QMatrix_americas_ARD_NoPoly_2Regions_062326.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s7.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s7_marginal,
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
plot(s7_marginal,
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
s8.tree<-read.nexus(file="BTV_Americas_s8_no_poly_aln_midpoint_rooted_nexus.tre")
s8.tree

# read in metadata
s8.data<-read.csv(file="A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s8.data.subset <- s8.data %>%
  filter(rownames(s8.data) %in% s8.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s8.tree, s8.data.subset)
summary(chk)

# create a named vector
s8.trait <- setNames(s8.data.subset[,"trait"], rownames(s8.data.subset))
unique(s8.trait)

plotTree(s8.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s8.trait),names(s8.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s8.tree),2,
             dimnames=list(s8.tree$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s8.tree)) 
  pies[s8.tree$tip.label[i],
       X[[s8.tree$tip.label[i]]]]<-
  rep(1/length(X[[s8.tree$tip.label[i]]]),
      length(X[[s8.tree$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s8.tree$tip.label, 
          tip = 1:Ntip(s8.tree), 
          adj = -0.05, 
          frame = "none", 
          cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(cols),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)

#######################################################

# fit model
# Since there no polymorphic states present use polyMk
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
# s8.ard<-fitMk(s8.tree,s8.trait,model="ARD")
# s8.ard
# save the object
# saveRDS(s8.ard, file = "s8_americas_ARD_NoPoly_2Regions_062326.rds")
# restore the object
s8.ard <- readRDS("s8_americas_ARD_NoPoly_2Regions_062326.rds")

# fit other models

### SYM
# s8.sym<-fitMk(s8.tree,s8.trait,model="SYM")
# s8.sym
# save the object
# saveRDS(s8.sym, file = "s8_americas_SYM_NoPoly_2Regions_062326.rds")
# restore the object
s8.sym <- readRDS("s8_americas_SYM_NoPoly_2Regions_062326.rds")

### ER
# s8.er<-fitMk(s8.tree,s8.trait,model="ER")
# s8.er
# save the object
# saveRDS(s8.er, file = "s8_americas_ER_NoPoly_2Regions_062326.rds")
# restore the object
s8.er <- readRDS("s8_americas_ER_NoPoly_2Regions_062326.rds")

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s8.ard),
                    logLik(s8.sym),
                    logLik(s8.er)),
           AIC=c(AIC(s8.ard),
                 AIC(s8.sym),
                 AIC(s8.er)))

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit

# transition_model    logLik      AIC
# 1              ARD -40.95275 85.90550
# 2              SYM -41.81588 85.63176 ***
# 3               ER -41.81588 85.63176 ***

# ARD fits the situation the best ??? Not much of a difference

# plot the models
# plot(s8.ard)
# plot(s8.sym)
# plot(s8.er)

####################################################### 


#### ******* SYM AND ER FIT BETTER ******* ####

# Marginal ASR on ARD model
# Still using ARD model even though it was not the best fit!!!
# Not much of a difference between ARD/ER/SYM
####################################################### 
# marginal ancestral reconstruction
s8_marginal<-ancr(s8.ard, type="marginal")
s8_marginal
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s8.ard)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s8_QMatrix_americas_ARD_NoPoly_2Regions_062326.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s8.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s8_marginal,
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
plot(s8_marginal,
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
s9.tree<-read.nexus(file="BTV_Americas_s9_no_poly_aln_midpoint_rooted_nexus.tre")
s9.tree

# read in metadata
s9.data<-read.csv(file="A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s9.data.subset <- s9.data %>%
  filter(rownames(s9.data) %in% s9.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s9.tree, s9.data.subset)
summary(chk)

# create a named vector
s9.trait <- setNames(s9.data.subset[,"trait"], rownames(s9.data.subset))
unique(s9.trait)

plotTree(s9.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s9.trait),names(s9.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s9.tree),2,
             dimnames=list(s9.tree$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s9.tree)) 
  pies[s9.tree$tip.label[i],
       X[[s9.tree$tip.label[i]]]]<-
  rep(1/length(X[[s9.tree$tip.label[i]]]),
      length(X[[s9.tree$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s9.tree$tip.label, 
          tip = 1:Ntip(s9.tree), 
          adj = -0.05, 
          frame = "none", 
          cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(cols),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)

#######################################################

# fit model
# Since there no polymorphic states present use polyMk
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
# s9.ard<-fitMk(s9.tree,s9.trait,model="ARD")
# s9.ard
# save the object
# saveRDS(s9.ard, file = "s9_americas_ARD_NoPoly_2Regions_062326.rds")
# restore the object
s9.ard <- readRDS("s9_americas_ARD_NoPoly_2Regions_062326.rds")

# fit other models

### SYM
# s9.sym<-fitMk(s9.tree,s9.trait,model="SYM")
# s9.sym
# save the object
# saveRDS(s9.sym, file = "s9_americas_SYM_NoPoly_2Regions_062326.rds")
# restore the object
s9.sym <- readRDS("s9_americas_SYM_NoPoly_2Regions_062326.rds")

### ER
# s9.er<-fitMk(s9.tree,s9.trait,model="ER")
# s9.er
# save the object
# saveRDS(s9.er, file = "s9_americas_ER_NoPoly_2Regions_062326.rds")
# restore the object
s9.er <- readRDS("s9_americas_ER_NoPoly_2Regions_062326.rds")

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s9.ard),
                    logLik(s9.sym),
                    logLik(s9.er)),
           AIC=c(AIC(s9.ard),
                 AIC(s9.sym),
                 AIC(s9.er)))

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit

# transition_model    logLik      AIC
# 1              ARD -47.99576 99.99153
# 2              SYM -48.07364 98.14729 ****
# 3               ER -48.07364 98.14729 ****

# ARD fits the situation the best ??? Not much of a difference

# plot the models
# plot(s9.ard)
# plot(s9.sym)
# plot(s9.er)

####################################################### 


#### ******* SYM AND ER FIT BETTER ******* ####

# Marginal ASR on ARD model
# Still using ARD model even though it was not the best fit!!!
# Not much of a difference between ARD/ER/SYM
####################################################### 
# marginal ancestral reconstruction
s9_marginal<-ancr(s9.ard, type="marginal")
s9_marginal
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s9.ard)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s9_QMatrix_americas_ARD_NoPoly_2Regions_062326.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s9.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s9_marginal,
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
state_order <- sort(unique(s9.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s9_marginal,
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
s10.tree<-read.nexus(file="BTV_Americas_s10_no_poly_aln_midpoint_rooted_nexus.tre")
s10.tree

# read in metadata
s10.data<-read.csv(file="A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s10.data.subset <- s10.data %>%
  filter(rownames(s10.data) %in% s10.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s10.tree, s10.data.subset)
summary(chk)

# create a named vector
s10.trait <- setNames(s10.data.subset[,"trait"], rownames(s10.data.subset))
unique(s10.trait)

plotTree(s10.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3)
X<-strsplit(setNames(as.character(s10.trait),names(s10.trait)),"+",
            fixed=TRUE)
pies<-matrix(0,Ntip(s10.tree),2,
             dimnames=list(s10.tree$tip.label,
                           c("Outside USA", "USA")))
for(i in 1:Ntip(s10.tree)) 
  pies[s10.tree$tip.label[i],
       X[[s10.tree$tip.label[i]]]]<-
  rep(1/length(X[[s10.tree$tip.label[i]]]),
      length(X[[s10.tree$tip.label[i]]]))
cols<-c("red","blue")
par(fg="transparent")
tiplabels(pie=pies,piecol=cols,cex=0.15)
tiplabels(text = s10.tree$tip.label, 
          tip = 1:Ntip(s10.tree), 
          adj = -0.05, 
          frame = "none", 
          cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
legend("left",
       legend = my_legend_order,
       pch = 16, col = as.vector(cols),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)

#######################################################

# fit model
# Since there no polymorphic states present use polyMk
# Decided to use default pi
# Default pi = equal: Want to give each state an equal chance at the root (no idea what root is)
####################################################### 

### ARD
# s10.ard<-fitMk(s10.tree,s10.trait,model="ARD")
# s10.ard
# save the object
# saveRDS(s10.ard, file = "s10_americas_ARD_NoPoly_2Regions_062326.rds")
# restore the object
s10.ard <- readRDS("s10_americas_ARD_NoPoly_2Regions_062326.rds")

# fit other models

### SYM
# s10.sym<-fitMk(s10.tree,s10.trait,model="SYM")
# s10.sym
# save the object
# saveRDS(s10.sym, file = "s10_americas_SYM_NoPoly_2Regions_062326.rds")
# restore the object
s10.sym <- readRDS("s10_americas_SYM_NoPoly_2Regions_062326.rds")

### ER
# s10.er<-fitMk(s10.tree,s10.trait,model="ER")
# s10.er
# save the object
# saveRDS(s10.er, file = "s10_americas_ER_NoPoly_2Regions_062326.rds")
# restore the object
s10.er <- readRDS("s10_americas_ER_NoPoly_2Regions_062326.rds")

# compare models 
data.frame(transition_model=c("ARD","SYM","ER"),
           logLik=c(logLik(s10.ard),
                    logLik(s10.sym),
                    logLik(s10.er)),
           AIC=c(AIC(s10.ard),
                 AIC(s10.sym),
                 AIC(s10.er)))

# RESULTS
# lower AIC values mean a better model fit
# Use AIC - Penalizes more complicated models and is a better identifier of model fit

# transition_model    logLik      AIC
# 1              ARD -42.22785 88.45569
# 2              SYM -43.16139 88.32278 ***
# 3               ER -43.16139 88.32278 ***

# ARD fits the situation the best ??? Not much of a difference

# plot the models
# plot(s10.ard)
# plot(s10.sym)
# plot(s10.er)

####################################################### 


#### ******* SYM AND ER FIT BETTER ******* ####

# Marginal ASR on ARD model
# Still using ARD model even though it was not the best fit!!!
# Not much of a difference between ARD/ER/SYM
####################################################### 
# marginal ancestral reconstruction
s10_marginal<-ancr(s10.ard, type="marginal")
s10_marginal
####################################################### 

# pull rates from Q matrix
#######################################################
Q <- as.Qmatrix(s10.ard)
rownames(Q)
colnames(Q)
Q
# write.csv(Q, file = "s10_QMatrix_americas_ARD_NoPoly_2Regions_062326.csv", row.names = TRUE)
#######################################################

# Plot Marginal ASR
####################################################### 
# plot marginal only
state_order <- sort(unique(s10.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s10_marginal,
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
state_order <- sort(unique(s10.trait))
my_colors <- c("red", "blue")
cols <- setNames(my_colors[seq_along(state_order)], state_order)
plot(s10_marginal,
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
# END SEGMENT 10
#
########################