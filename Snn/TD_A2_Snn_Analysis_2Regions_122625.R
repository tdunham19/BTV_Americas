#######################################################
# Nearest Neighbor Statistic Analysis 
# 
# This script reads in sets of sequences and calculates Snn
# Snn = nearest neighbor statistics 
# according to Hudson (2011) 
#
# Hudson RR. A new statistic for detecting genetic differentiation. 
# Genetics. 2000 Aug;155(4):2011-4. 
# doi: 10.1093/genetics/155.4.2011. PMID: 10924493; PMCID: PMC1461195.
#
# Original script written by: Mark Stenglein 8/28/2025
# Script adapted by: Tillie Dunham 09/23/2025 for Aim 2 work
# TD edited on 12.26.25 to perform Snn on 2 regions from ASR analysis
# all samples: USA and Outside USA

#######################################################

#Load Required Libraries
library(tidyverse)
library(ape)
library(profvis)
library(readxl)
library(ggplot2)

#######################################################


# Segment 1 Snn - BTV Americas All samples
################################################################################
# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read.csv("AIM2_Metadata_Nextstrain_renamed_traits_Americas_2_regions.csv")

# pull out just the columns we need
metadata <- metadata %>% select(accession = accession, group)

num_permutations <- 5000

process_alignment <-  function (
    fasta_msa = "BTV_Americas_s1_aln.fasta", 
    prefix="Segment 1") {
  
  # DEBUG
  # fasta_msa = "BTV_Americas_s1_aln.fasta"
  # prefix="Segment 1"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  # accessions <- str_extract(labels(msa), ".*?(?=_-_)")
  accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment1_snn <- process_alignment("BTV_Americas_s1_aln.fasta", "Segment 1")

Segment1_snn$plot


# END Segment 1 Snn
################################################################################




# Segment 2 Snn - BTV Americas All samples
################################################################################
# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read.csv("AIM2_Metadata_Nextstrain_renamed_traits_Americas_2_regions.csv")

# pull out just the columns we need
metadata <- metadata %>% select(accession = accession, group)

num_permutations <- 5000

process_alignment <-  function (
    fasta_msa = "BTV_Americas_s2_aln.fasta", 
    prefix="Segment 2") {
  
  # DEBUG
  # fasta_msa = "BTV_Americas_s2_aln.fasta"
  # prefix="Segment 2"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  # accessions <- str_extract(labels(msa), ".*?(?=_-_)")
  accessions <- str_extract(labels(msa), ".*?(?=\\s)")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment2_snn <- process_alignment("BTV_Americas_s2_aln.fasta", "Segment 2")

Segment2_snn$plot


# END Segment 2 Snn
################################################################################




# Segment 3 Snn - BTV Americas All samples
################################################################################
# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read.csv("AIM2_Metadata_Nextstrain_renamed_traits_Americas_2_regions.csv")

# pull out just the columns we need
metadata <- metadata %>% select(accession = accession, group)

num_permutations <- 5000

process_alignment_s3 <-  function (
    fasta_msa = "BTV_Americas_s3_aln.fasta", 
    prefix="Segment 3") {
  
  # DEBUG
  # fasta_msa = "BTV_Americas_s3_aln.fasta"
  # prefix="Segment 3"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  accessions <- str_extract(labels(msa), ".*?(?=\\s)")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn_s3(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn_s3(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn_s3 <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment3_snn <- process_alignment_s3("BTV_Americas_s3_aln.fasta", "Segment 3")

Segment3_snn$plot


# END Segment 3 Snn
################################################################################




# Segment 4 Snn - BTV Americas All samples
################################################################################
# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read.csv("AIM2_Metadata_Nextstrain_renamed_traits_Americas_2_regions.csv")

# pull out just the columns we need
metadata <- metadata %>% select(accession = accession, group)

num_permutations <- 5000

process_alignment <-  function (
    fasta_msa = "BTV_Americas_s4_aln.fasta", 
    prefix="Segment 4") {
  
  # DEBUG
  # fasta_msa = "BTV_Americas_s4_aln.fasta"
  # prefix="Segment 4"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  # accessions <- str_extract(labels(msa), ".*?(?=_-_)")
  accessions <- str_extract(labels(msa), ".*?(?=\\s)")
  # accessions <- str_extract(labels(msa), ".*")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment4_snn <- process_alignment("BTV_Americas_s4_aln.fasta", "Segment 4")

Segment4_snn$plot


# END Segment 4 Snn
################################################################################




# Segment 5 Snn - BTV Americas All samples
################################################################################
# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read.csv("AIM2_Metadata_Nextstrain_renamed_traits_Americas_2_regions.csv")

# pull out just the columns we need
metadata <- metadata %>% select(accession = accession, group)

num_permutations <- 5000

process_alignment <-  function (
    fasta_msa = "BTV_Americas_s5_aln.fasta", 
    prefix="Segment 5") {
  
  # DEBUG
  # fasta_msa = "BTV_Americas_s5_aln.fasta"
  # prefix="Segment 5"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  # accessions <- str_extract(labels(msa), ".*?(?=_-_)")
  # accessions <- str_extract(labels(msa), ".*")
  accessions <- str_extract(labels(msa), ".*?(?=\\s)")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment5_snn <- process_alignment("BTV_Americas_s5_aln.fasta", "Segment 5")

Segment5_snn$plot


# END Segment 5 Snn
################################################################################




# Segment 6 Snn - BTV Americas All samples
################################################################################
# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read.csv("AIM2_Metadata_Nextstrain_renamed_traits_Americas_2_regions.csv")

# pull out just the columns we need
metadata <- metadata %>% select(accession = accession, group)

num_permutations <- 5000

process_alignment <-  function (
    fasta_msa = "BTV_Americas_s6_aln.fasta", 
    prefix="Segment 6") {
  
  # DEBUG
  # fasta_msa = "BTV_Americas_s6_aln.fasta"
  # prefix="Segment 6"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  # accessions <- str_extract(labels(msa), ".*?(?=_-_)")
  # accessions <- str_extract(labels(msa), ".*")
  accessions <- str_extract(labels(msa), ".*?(?=\\s)")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment6_snn <- process_alignment("BTV_Americas_s6_aln.fasta", "Segment 6")

Segment6_snn$plot


# END Segment 6 Snn
################################################################################




# Segment 7 Snn - BTV Americas All samples
################################################################################
# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read.csv("AIM2_Metadata_Nextstrain_renamed_traits_Americas_2_regions.csv")

# pull out just the columns we need
metadata <- metadata %>% select(accession = accession, group)

num_permutations <- 5000

process_alignment <-  function (
    fasta_msa = "BTV_Americas_s7_aln.fasta", 
    prefix="Segment 7") {
  
  # DEBUG
  # fasta_msa = "BTV_Americas_s7_aln.fasta"
  # prefix="Segment 7"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  # accessions <- str_extract(labels(msa), ".*?(?=_-_)")
  # accessions <- str_extract(labels(msa), ".*")
  accessions <- str_extract(labels(msa), ".*?(?=\\s)")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment7_snn <- process_alignment("BTV_Americas_s7_aln.fasta", "Segment 7")

Segment7_snn$plot


# END Segment 7 Snn
################################################################################




# Segment 8 Snn - BTV Americas All samples
################################################################################
# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read.csv("AIM2_Metadata_Nextstrain_renamed_traits_Americas_2_regions.csv")

# pull out just the columns we need
metadata <- metadata %>% select(accession = accession, group)

num_permutations <- 5000

process_alignment <-  function (
    fasta_msa = "BTV_Americas_s8_aln.fasta", 
    prefix="Segment 8") {
  
  # DEBUG
  # fasta_msa = "BTV_Americas_s8_aln.fasta"
  # prefix="Segment 8"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  # accessions <- str_extract(labels(msa), ".*?(?=_-_)")
  # accessions <- str_extract(labels(msa), ".*")
  accessions <- str_extract(labels(msa), ".*?(?=\\s)")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment8_snn <- process_alignment("BTV_Americas_s8_aln.fasta", "Segment 8")

Segment8_snn$plot


# END Segment 8 Snn
################################################################################




# Segment 9 Snn - BTV Americas All samples
################################################################################
# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read.csv("AIM2_Metadata_Nextstrain_renamed_traits_Americas_2_regions.csv")

# pull out just the columns we need
metadata <- metadata %>% select(accession = accession, group)

num_permutations <- 5000

process_alignment <-  function (
    fasta_msa = "BTV_Americas_s9_aln.fasta", 
    prefix="Segment 9") {
  
  # DEBUG
  # fasta_msa = "BTV_Americas_s9_aln.fasta"
  # prefix="Segment 9"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  # accessions <- str_extract(labels(msa), ".*?(?=_-_)")
  # accessions <- str_extract(labels(msa), ".*")
  accessions <- str_extract(labels(msa), ".*?(?=\\s)")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment9_snn <- process_alignment("BTV_Americas_s9_aln.fasta", "Segment 9")

Segment9_snn$plot


# END Segment 9 Snn
################################################################################




# Segment 10 Snn - BTV Americas All samples
################################################################################
# -------------------------------
# Metadata (sample status)
# -------------------------------

# read in sample metadata including status
metadata <- read.csv("AIM2_Metadata_Nextstrain_renamed_traits_Americas_2_regions.csv")

# pull out just the columns we need
metadata <- metadata %>% select(accession = accession, group)

num_permutations <- 5000

process_alignment <-  function (
    fasta_msa = "BTV_Americas_s10_aln.fasta", 
    prefix="Segment 10") {
  
  # DEBUG
  # fasta_msa = "BTV_Americas_s10_aln.fasta"
  # prefix="Segment 10"
  
  # read in sequences
  msa <- read.dna(fasta_msa, format="fasta")
  
  # relabel sequences with just the first part
  # accessions <- str_extract(labels(msa), ".*?(?=_-_)")
  # accessions <- str_extract(labels(msa), ".*")
  accessions <- str_extract(labels(msa), ".*?(?=\\s)")
  
  # update labels: shorten to just accessions
  rownames(msa) <- accessions
  
  # drop sequences for which we could not identify an accession
  msa <- msa[!(is.na(accessions)), ]
  accessions <- accessions[!(is.na(accessions))]
  
  # calculate distances
  # TN93 distance model
  dist_model = "TN93" 
  msa_dist <- as.matrix(dist.dna(msa, model = dist_model))
  
  acc_loc <- tibble(accession = rownames(msa_dist))
  # this will make groups in order of matrix accession
  acc_loc <- left_join(acc_loc, metadata)  
  # create a vector of the group matching the order of the distance matrix
  groups <- acc_loc %>% pull(group)
  
  # calculate Snn from distance matrix and sample metadata
  this_snn <- calculate_snn(msa_dist, groups)
  
  perm_snns <- c()
  
  # do permutation testing
  for (p in 1:num_permutations) {
    # first, scramble groups using sampling without replacement
    # using sampling without replacement will keep the # of samples from each
    # location (group) the same between permutations, as specified in Hudson 2011
    scrambled_groups  <- sample(groups, replace=F)
    perm_snn             <- calculate_snn(msa_dist, scrambled_groups)
    perm_snns            <- c(perm_snns, perm_snn)
  }
  
  # calculate p-value from Snns from permutation tests
  this_snn_p_val <- sum(perm_snns >= this_snn) / num_permutations
  
  if (this_snn_p_val == 0) {
    this_snn_p_val_sci <- sprintf("%0.1e", 1/num_permutations)
    this_snn_p_val <- paste0("<", this_snn_p_val_sci)
  }
  
  snn_text <- paste0(prefix, 
                     " ", 
                     sprintf("%0.3f", this_snn), 
                     " p = ", 
                     this_snn_p_val)
  print(snn_text)
  
  # plot a histogram of permutation pvalues and this 
  perm_snns_tibble <- tibble(snn = perm_snns)
  plot_snn <- ggplot(perm_snns_tibble) +
    geom_histogram(aes(x=snn), bins=100, fill="grey20", color=NA) +
    geom_vline(xintercept = this_snn, color="red", linewidth = 0.5) +
    annotate("text", x=0.4, y=200, label=snn_text, hjust = 0, size=3) +
    xlab("Snn") +
    ylab("Count") +
    theme_bw(base_size = 14) 
  
  ggsave(paste0(prefix, "_Snn_permutation_testing.pdf"), 
         units="in", width=7.5, height=5)
  
  
  return(list(snn    = this_snn,
              pvalue = this_snn_p_val,
              plot   = plot_snn,
              text   = snn_text))
}

calculate_snn <- function(distance_matrix, groups_vector) {
  
  Snn <- 0
  
  # DEBUG
  # distance_matrix = btv_msa_dist
  # groups_vector = group_loc
  
  # iterate through samples to calculate Snn
  # there is probably a more R way to do this but I am doing it this way
  # because it makes sense to me
  num_samples <- nrow(distance_matrix)
  
  # Loop through each individual
  for (i in 1:num_samples) {
    
    # get the row of distances for this sample
    row = distance_matrix[i,]
    
    # exclude self
    row[i] <- NA  
    
    # for each row, find the minimum distance (ignore NA)
    min_val <- min(row, na.rm = TRUE)
    
    # this is a vector of indices of nearest neighbors, with distance == minimum_distance
    nn_i <- which(row == min_val)
    
    # this is the number of nearest neighbors with the same group as sample i
    num_shared_group <- sum(groups_vector[nn_i] == groups_vector[i])
    
    # this is the number of nearest neighbors (regardless of group)
    num_NN         <- length(nn_i)
    
    # the contribution to SNN for this sample = num_NN_with_shared_loc / num_NN / num_samples
    fraction_shared_group <- (num_shared_group / num_NN) / num_samples
    
    if (is.na(fraction_shared_group)) {
      message(paste0("Warning, sample ", names(row)[i], " produced an NA value"))
    }
    
    Snn <- Snn + fraction_shared_group
  }
  
  # return Snn
  Snn 
}

Segment10_snn <- process_alignment("BTV_Americas_s10_aln.fasta", "Segment 10")

Segment10_snn$plot


# END Segment 10 Snn
################################################################################




# Putting it all together 
################################################################################
# -------------------------------
# Create table
# -------------------------------
snn_supp_table <- 
  tibble( segment = c("BTV Segment 1", "BTV Segment 2", "BTV Segment 3", 
                      "BTV Segment 4", "BTV Segment 5", "BTV Segment 6", 
                      "BTV Segment 7", "BTV Segment 8", "BTV Segment 9", "BTV Segment 10"),
          snn  = c(Segment1_snn$snn, Segment2_snn$snn, Segment3_snn$snn, 
                   Segment4_snn$snn, Segment5_snn$snn, Segment6_snn$snn, 
                   Segment7_snn$snn, Segment8_snn$snn, Segment9_snn$snn, Segment10_snn$snn),
          pval = c(Segment1_snn$pval, Segment2_snn$pval, Segment3_snn$pval, 
                   Segment4_snn$pval, Segment5_snn$pval, Segment6_snn$pval, 
                   Segment7_snn$pval, Segment8_snn$pval, Segment9_snn$pval, Segment10_snn$pval))

write.table(snn_supp_table, 
            file="BTV_Americas_2Regions_all_segments_supplemental_table_snn_122625.txt",
            quote=F, sep="\t", row.names=F)


# -------------------------------
# Adjust p-values
# -------------------------------

# Adjust p-values

# format p-values as numeric
Segment1_snn$pval <- as.numeric(gsub("<", "", Segment1_snn$pval))
Segment2_snn$pval <- as.numeric(gsub("<", "", Segment2_snn$pval))
Segment3_snn$pval <- as.numeric(gsub("<", "", Segment3_snn$pval))
Segment4_snn$pval <- as.numeric(gsub("<", "", Segment4_snn$pval))
Segment5_snn$pval <- as.numeric(gsub("<", "", Segment5_snn$pval))
Segment6_snn$pval <- as.numeric(gsub("<", "", Segment6_snn$pval))
Segment7_snn$pval <- as.numeric(gsub("<", "", Segment7_snn$pval))
Segment8_snn$pval <- as.numeric(gsub("<", "", Segment8_snn$pval))
Segment9_snn$pval <- as.numeric(gsub("<", "", Segment9_snn$pval))
Segment10_snn$pval <- as.numeric(gsub("<", "", Segment10_snn$pval))

# create table and apply bonferroni correction
pval_supp_table <- 
  tibble( segment = c("s1", "s2", "s3", 
                      "s4", "s5", "s6", 
                      "s7", "s8", "s9", "s10"),
          snn  = c(Segment1_snn$snn, Segment2_snn$snn, Segment3_snn$snn, 
                   Segment4_snn$snn, Segment5_snn$snn, Segment6_snn$snn, 
                   Segment7_snn$snn, Segment8_snn$snn, Segment9_snn$snn, Segment10_snn$snn),
          pval = c(Segment1_snn$pval, Segment2_snn$pval, Segment3_snn$pval, 
                   Segment4_snn$pval, Segment5_snn$pval, Segment6_snn$pval, 
                   Segment7_snn$pval, Segment8_snn$pval, Segment9_snn$pval, Segment10_snn$pval)
  ) %>%
  mutate(pval_bonf = p.adjust(pval, method = "bonferroni", n = length(pval)))

# save table
# write.table(pval_supp_table, 
#             file="BTV_Americas_2Regions_all_segments_pval_adjusted_supplemental_table_snn_122625.txt",
#             quote=F, sep="\t", row.names=F)

################################################################################

