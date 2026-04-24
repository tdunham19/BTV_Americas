#######################################################
# Aim 2 - Phylogenetic Tree Visualization
# Plotting all segment trees on same scale
# TD 02.05.26
#######################################################

# load required libraries
library(phytools)
library(tidyverse)

#######################################################


#### Segment 1 ####
# plot tree
#######################################################
# read in tree
s1.tree<-read.nexus(file="../../ASR/Americas/BTV_Americas_s1_no_poly_aln_midpoint_rooted_nexus.tre")
s1.tree

# read in metadata
s1.data<-read.csv(file="../../ASR/Americas/A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s1.data.subset <- s1.data %>%
  filter(rownames(s1.data) %in% s1.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s1.tree, s1.data.subset)
summary(chk)

# create a named vector
s1.trait <- setNames(s1.data.subset[,"trait"], rownames(s1.data.subset))
unique(s1.trait)

# plotTree(s1.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3,x.lim = c(0,20))
plot(s1.tree, x.lim = c(0.1,5.75),show.tip.label = FALSE)
add.scale.bar(x = 0, y = -2, cex = 0.5)
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
# tiplabels(text = s1.tree$tip.label, tip = 1:Ntip(s1.tree), adj = -0.05, 
#          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
# legend("topleft",
#        legend = my_legend_order,
#        pch = 16, col = as.vector(cols),
#        horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
#        y.intersp = 1.2)

# Save as 8 x 8 pdf
#######################################################


#### Segment 2 ####
# plot tree
#######################################################
# read in tree
s2.tree<-read.nexus(file="../../ASR/Americas/BTV_Americas_s2_no_poly_aln_midpoint_rooted_nexus.tre")
s2.tree

# read in metadata
s2.data<-read.csv(file="../../ASR/Americas/A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s2.data.subset <- s2.data %>%
  filter(rownames(s2.data) %in% s2.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s2.tree, s2.data.subset)
summary(chk)

# create a named vector
s2.trait <- setNames(s2.data.subset[,"trait"], rownames(s2.data.subset))
unique(s2.trait)

# plotTree(s2.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3,x.lim = c(0,20))
plot(s2.tree, x.lim = c(0.1,5.75),show.tip.label = FALSE)
add.scale.bar(x = 0, y = -2, cex = 0.5)
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
# tiplabels(text = s2.tree$tip.label, tip = 1:Ntip(s2.tree), adj = -0.05, 
#          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
legend("topleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(cols),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)

# Save as 8 x 8 pdf
#######################################################


#### Segment 3 ####
# plot tree
#######################################################
# read in tree
s3.tree<-read.nexus(file="../../ASR/Americas/BTV_Americas_s3_no_poly_aln_midpoint_rooted_nexus.tre")
s3.tree

# read in metadata
s3.data<-read.csv(file="../../ASR/Americas/A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s3.data.subset <- s3.data %>%
  filter(rownames(s3.data) %in% s3.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s3.tree, s3.data.subset)
summary(chk)

# create a named vector
s3.trait <- setNames(s3.data.subset[,"trait"], rownames(s3.data.subset))
unique(s3.trait)

# plotTree(s3.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3,x.lim = c(0,20))
plot(s3.tree, x.lim = c(0.1,5.75),show.tip.label = FALSE)
add.scale.bar(x = 0, y = -2, cex = 0.5)
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
# tiplabels(text = s3.tree$tip.label, tip = 1:Ntip(s3.tree), adj = -0.05, 
#          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
# legend("topleft",
#        legend = my_legend_order,
#        pch = 16, col = as.vector(cols),
#        horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
#        y.intersp = 1.2)

# Save as 8 x 8 pdf
#######################################################


#### Segment 4 ####
# plot tree
#######################################################
# read in tree
s4.tree<-read.nexus(file="../../ASR/Americas/BTV_Americas_s4_no_poly_aln_midpoint_rooted_nexus.tre")
s4.tree

# read in metadata
s4.data<-read.csv(file="../../ASR/Americas/A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s4.data.subset <- s4.data %>%
  filter(rownames(s4.data) %in% s4.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s4.tree, s4.data.subset)
summary(chk)

# create a named vector
s4.trait <- setNames(s4.data.subset[,"trait"], rownames(s4.data.subset))
unique(s4.trait)

# plotTree(s4.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3,x.lim = c(0,20))
plot(s4.tree, x.lim = c(0.1,5.75),show.tip.label = FALSE)
add.scale.bar(x = 0, y = -2, cex = 0.5)
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
# tiplabels(text = s4.tree$tip.label, tip = 1:Ntip(s4.tree), adj = -0.05, 
#          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
# legend("topleft",
#        legend = my_legend_order,
#        pch = 16, col = as.vector(cols),
#        horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
#        y.intersp = 1.2)

# Save as 8 x 8 pdf
#######################################################


#### Segment 5 ####
# plot tree
#######################################################
# read in tree
s5.tree<-read.nexus(file="../../ASR/Americas/BTV_Americas_s5_no_poly_aln_midpoint_rooted_nexus.tre")
s5.tree

# read in metadata
s5.data<-read.csv(file="../../ASR/Americas/A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s5.data.subset <- s5.data %>%
  filter(rownames(s5.data) %in% s5.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s5.tree, s5.data.subset)
summary(chk)

# create a named vector
s5.trait <- setNames(s5.data.subset[,"trait"], rownames(s5.data.subset))
unique(s5.trait)

# plotTree(s5.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3,x.lim = c(0,20))
plot(s5.tree, x.lim = c(0.1,5.75),show.tip.label = FALSE)
add.scale.bar(x = 0, y = -2, cex = 0.5)
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
# tiplabels(text = s5.tree$tip.label, tip = 1:Ntip(s5.tree), adj = -0.05, 
#          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
# legend("topleft",
#        legend = my_legend_order,
#        pch = 16, col = as.vector(cols),
#        horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
#        y.intersp = 1.2)

# Save as 8 x 8 pdf
#######################################################


#### Segment 6 ####
# plot tree
#######################################################
# read in tree
s6.tree<-read.nexus(file="../../ASR/Americas/BTV_Americas_s6_no_poly_aln_midpoint_rooted_nexus.tre")
s6.tree

# read in metadata
s6.data<-read.csv(file="../../ASR/Americas/A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s6.data.subset <- s6.data %>%
  filter(rownames(s6.data) %in% s6.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s6.tree, s6.data.subset)
summary(chk)

# create a named vector
s6.trait <- setNames(s6.data.subset[,"trait"], rownames(s6.data.subset))
unique(s6.trait)

# plotTree(s6.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3,x.lim = c(0,20))
plot(s6.tree, x.lim = c(0.1,5.75),show.tip.label = FALSE)
add.scale.bar(x = 0, y = -2, cex = 0.5)
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
# tiplabels(text = s6.tree$tip.label, tip = 1:Ntip(s6.tree), adj = -0.05, 
#          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(cols),
       horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
       y.intersp = 1.2)

# Save as 8 x 8 pdf
#######################################################


#### Segment 7 ####
# plot tree
#######################################################
# read in tree
s7.tree<-read.nexus(file="../../ASR/Americas/BTV_Americas_s7_no_poly_aln_midpoint_rooted_nexus.tre")
s7.tree

# read in metadata
s7.data<-read.csv(file="../../ASR/Americas/A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s7.data.subset <- s7.data %>%
  filter(rownames(s7.data) %in% s7.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s7.tree, s7.data.subset)
summary(chk)

# create a named vector
s7.trait <- setNames(s7.data.subset[,"trait"], rownames(s7.data.subset))
unique(s7.trait)

# plotTree(s7.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3,x.lim = c(0,20))
plot(s7.tree, x.lim = c(0.1,5.75),show.tip.label = FALSE)
add.scale.bar(x = 0, y = -2, cex = 0.5)
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
# tiplabels(text = s7.tree$tip.label, tip = 1:Ntip(s7.tree), adj = -0.05, 
#          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
legend("bottomleft",
       legend = my_legend_order,
       pch = 16, col = as.vector(cols),
       horiz = FALSE, cex = 0.6, bty = "n", pt.cex = 2,
       y.intersp = 1.2)

# Save as 8 x 8 pdf
#######################################################


#### Segment 8 ####
# plot tree
#######################################################
# read in tree
s8.tree<-read.nexus(file="../../ASR/Americas/BTV_Americas_s8_no_poly_aln_midpoint_rooted_nexus.tre")
s8.tree

# read in metadata
s8.data<-read.csv(file="../../ASR/Americas/A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s8.data.subset <- s8.data %>%
  filter(rownames(s8.data) %in% s8.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s8.tree, s8.data.subset)
summary(chk)

# create a named vector
s8.trait <- setNames(s8.data.subset[,"trait"], rownames(s8.data.subset))
unique(s8.trait)

# plotTree(s8.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3,x.lim = c(0,20))
plot(s8.tree, x.lim = c(0.1,5.75),show.tip.label = FALSE)
add.scale.bar(x = 0, y = -2, cex = 0.5)
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
# tiplabels(text = s8.tree$tip.label, tip = 1:Ntip(s8.tree), adj = -0.05, 
#          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
# legend("topleft",
#        legend = my_legend_order,
#        pch = 16, col = as.vector(cols),
#        horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
#        y.intersp = 1.2)

# Save as 8 x 8 pdf
#######################################################


#### Segment 9 ####
# plot tree
#######################################################
# read in tree
s9.tree<-read.nexus(file="../../ASR/Americas/BTV_Americas_s9_no_poly_aln_midpoint_rooted_nexus.tre")
s9.tree

# read in metadata
s9.data<-read.csv(file="../../ASR/Americas/A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s9.data.subset <- s9.data %>%
  filter(rownames(s9.data) %in% s9.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s9.tree, s9.data.subset)
summary(chk)

# create a named vector
s9.trait <- setNames(s9.data.subset[,"trait"], rownames(s9.data.subset))
unique(s9.trait)

# plotTree(s9.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3,x.lim = c(0,20))
plot(s9.tree, x.lim = c(0.1,5.75),show.tip.label = FALSE)
add.scale.bar(x = 0, y = -2, cex = 0.5)
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
# tiplabels(text = s9.tree$tip.label, tip = 1:Ntip(s9.tree), adj = -0.05, 
#          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
# legend("topleft",
#        legend = my_legend_order,
#        pch = 16, col = as.vector(cols),
#        horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
#        y.intersp = 1.2)

# Save as 8 x 8 pdf
#######################################################


#### Segment 10 ####
# plot tree
#######################################################
# read in tree
s10.tree<-read.nexus(file="../../ASR/Americas/BTV_Americas_s10_no_poly_aln_midpoint_rooted_nexus.tre")
s10.tree

# read in metadata
s10.data<-read.csv(file="../../ASR/Americas/A2_Americas_Metadata_No_Polytomy.csv", row.names=1)

# keep only the tips that are found in the tree
s10.data.subset <- s10.data %>%
  filter(rownames(s10.data) %in% s10.tree$tip.label)

# checking to see mismatches between data and tree
chk<-geiger::name.check(s10.tree, s10.data.subset)
summary(chk)

# create a named vector
s10.trait <- setNames(s10.data.subset[,"trait"], rownames(s10.data.subset))
unique(s10.trait)

# plotTree(s10.tree,ftype="off",lwd=1,type="phylogram",offset = 0.03,fsize=0.3,x.lim = c(0,20))
plot(s10.tree, x.lim = c(0.1,5.75),show.tip.label = FALSE)
add.scale.bar(x = 0, y = -2, cex = 0.5)
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
# tiplabels(text = s10.tree$tip.label, tip = 1:Ntip(s10.tree), adj = -0.05, 
#          frame = "none", cex = 0.5)
par(fg="black")
my_legend_order <- c("Outside USA", "USA")
# legend("topleft",
#        legend = my_legend_order,
#        pch = 16, col = as.vector(cols),
#        horiz = FALSE, cex = 0.8, bty = "n", pt.cex = 2,
#        y.intersp = 1.2)

# Save as 8 x 8 pdf
#######################################################
