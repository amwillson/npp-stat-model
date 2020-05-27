#!/usr/bin/Rscript
# Load packages to library
library(plotrix)
library(dplR)
library(fields)
library(reshape2)
library(plyr)

## Before running this script, download data from the Wiki for your site
## Then, add these data to the "data" folder in npp-stat-model under a folder for your site
## Put all rwl files into one subdirectory "rwl"
## And convert the meta file to CSV if it is not already in that format
## It wouldn't hurt to also check out these files to make sure you understand any aspects that are unique to that dataset

# source("config_HMC")
# Data directory, used to streamline data location between sites
dataDir = 'data/'
# Site name, must be the same as site-specific directories
site = 'SYLVANIA'
# Data and model versions
dvers = "v0.1"
mvers = "v0.1"

# Number of plots, specified in data
nPlots <- 1
# Conversion factor
ftPerMeter <- 3.2808399

# Specify first and last years
lastYear  <- 2015
firstYear <- 1940
# Interpolate years between first and last
years <- firstYear:lastYear
# Number of years
nT <- length(years)

# List of files with tree ring increments
rwFiles <- list.files(paste0(dataDir, site, '/', "rwl"))
# In case there are other files in this directory, take only rwl files
rwFiles <- rwFiles[grep(".rwl$", rwFiles)]
rwData <- list()
for(fn in rwFiles) {
    id <- gsub(".rw", "", fn)
    # Insert the contents of each file into the rwData list
    rwData[[id]] <- t(read.tucson(file.path("data", site, "rwl", fn)))  # rows are tree, cols are times
}

# Load CSV with tree data (DBH, location within plot, species)
treeMeta = read.csv(paste0(dataDir, site, '/alexander_sylvania_June2018.csv'))

# Make dataframe from the list of tree ring increments
incr = ldply(rwData, rbind)
# Insert core ID's into the data frame
incr = incr[,c(".id", sort(colnames(incr)[2:ncol(incr)]))]
rownames(incr) = as.vector(unlist(lapply(rwData, rownames)))
incr[,1] = rownames(incr)

######################################################################################################################################
## make nimble data
######################################################################################################################################
if (!file.exists('data/dump')){
  dir.create('data/dump')
} 

# Melt tree ring data frame  
incr_data = melt(incr)
colnames(incr_data) = c('id', 'year', 'incr')
# Assign plot to each core in data frame (in this case, all 1)
incr_data$plot   = rep(1, nrow(incr_data))#as.numeric(substr(incr_data$id, 3, 3))
# Assign tree ID to each core (= core ID)
incr_data$TreeID = incr_data$id 
# Assign ID to each core (= tree, not core, more than 1 core per tree)
incr_data$id     = as.numeric(substr(incr_data$id, 4, 6))
# Assign year to each observation (= dating from tree rings)
incr_data$year = as.vector(incr_data$year)

# Vector of individual tree IDs (removes directionality from ID)
tree_ids = unique(substr(incr_data$TreeID, 1, 6))
# Number of trees with cores
N_trees  = length(tree_ids)
# Vector sequence from 1 to the total number of trees to ID each core numerically
stat_ids = seq(1, N_trees)

# Add stat_ids to data frame by matching tree ID's (each tree gets a separate stat ID)
incr_data$stat_id = stat_ids[match(substr(incr_data$TreeID, 1, 6), tree_ids)]
for (n in 1:nrow(incr_data)){
  print(n)
  # Add species to tree ring data by matching species in whole-tree dataset with the ID
  incr_data$taxon[n] = as.vector(treeMeta$species[which(as.numeric(substr(treeMeta$ID, 4, 6)) == incr_data$id[n])])
}

# Number of species
N_taxa = length(unique(incr_data$taxon))
# Data frame of 4-letter species codes and corresponding numbers
taxaMatch = data.frame(species=sort(unique(incr_data$taxon)), number=seq(1, N_taxa))

# Vector 4-letter species code corresponding to each tree with tree ring observations
taxon = aggregate(taxon~stat_id, incr_data, unique)[,2]
# Convert to numbers corresponding to each species
taxon = taxaMatch$number[match(taxon, taxaMatch$species)]

# Vector of plot number corresponding to each tree (in this case, all 1's)
plot_id = aggregate(plot~stat_id, incr_data, unique)[,2]

##########################################################################################################################
## STAN DATA
##########################################################################################################################
# incr = incr[,-1]

# Last year of growth observed in tree rings
year_end = max(as.numeric(incr_data$year), na.rm=TRUE)
# First year of growth observed in tree rings
year_start = min(as.numeric(incr_data$year), na.rm=TRUE)
# Number of years in the series
N_years = year_end - year_start + 1
# Interpolate years
years = seq(year_start, year_end)

# order by tree and year (for each tree ring)
incr_data = incr_data[order(incr_data$stat_id, incr_data$year),]
# Remove all NA values (only keep the series for which trees were alive/existed)
incr_data = incr_data[which(!is.na(incr_data$incr)),]
# Make units of time in "years since first recorded tree ring"
incr_data$year = as.numeric(incr_data$year) - year_start + 1
# Number of measurements
N_inc   = nrow(incr_data)
# Vector of years for each increment
m2t     = incr_data$year
# Vector of tree ID's for each increment
m2tree  = incr_data$stat_id
# Vector of plot for each increment (all 1's)
m2plot  = incr_data$plot
# Vector of species numbers for each increment
m2taxon = taxaMatch$number[match(incr_data$taxon, taxaMatch$species)]
# Vector of increments (actual data)
Xobs    = incr_data$incr
# If the increment is 0 (i.e. no growth), make it a really small increment to fit the assuptions of the log transform
Xobs[Xobs==0] = 0.0001
# Log transform
logXobs = log(Xobs)

# Data frame of the first and last year of observation (years corresponding to first and last increment) for each tree
year_idx = data.frame(year_start=as.numeric(aggregate(year~stat_id, data=incr_data, FUN=min, na.rm=TRUE)[,2]), 
                      year_end=as.numeric(aggregate(year~stat_id, incr_data, max)[,2]))
# year_idx[,2] = rep(N_years, nrow(year_idx))
# year_idx = year_idx - year_start + 1

# make pdbh
# Make data frame of tree ID number from whole-tree dataset, sequential tree ID, and last year of observation
pdbh = aggregate(year~stat_id+id, incr_data, max, na.rm=TRUE)
# Order by tree ID number
pdbh = pdbh[order(pdbh$stat_id),]

for (n in 1:nrow(pdbh)){
  # Add DBH from whole-tree dataset using ID to match observations to those in increment dataset
  pdbh$dbh[n] = treeMeta$dbh[which(as.numeric(substr(treeMeta$ID, 4, 6)) == pdbh$id[n])]
  # Same with distance (location within plot)
  pdbh$distance[n] = treeMeta$distance[which(as.numeric(substr(treeMeta$ID, 4, 6)) == pdbh$id[n])]
}
# Number of trees with DBH observations
N_pdbh = nrow(pdbh)
# Log transform of DBH
logPDobs = log(pdbh$dbh)
# Vector of Tree ID's
pdbh_tree_id = pdbh$stat_id
# logPDobs[is.na(logPDobs)] = -999
# Vector of year of observation for each tree
pdbh_year_id = pdbh$year
# Vector of locations within plot
distance = pdbh$distance

idx_stack = data.frame(meas=numeric(0), tree_id=numeric(0), year=numeric(0))
n = 1
for (tree in 1:N_trees){
  # Vector of years of increment observations for a given tree
  year = seq(year_idx[tree,1], year_idx[tree,2])
  # Vector of observation index for each increment of a given tree
  meas = seq(n, n+length(year)-1)
  # Increment counter
  n = n + length(year)
  # Add to a new data frame
  idx_stack = rbind(idx_stack, data.frame(meas=meas, tree_id=rep(tree, length(year)), year=year))
}

# Find location of first observation for each tree
idx_tree = which(!duplicated(idx_stack$tree_id))
# Data frame of first and last observation for each tree
idx_tree = data.frame(idx_tree, c(idx_tree[-1]-1, nrow(idx_stack)))

# Vector of tree IDs from the new data frame
x2tree  = idx_stack$tree_id
# Vector of years from the new data frame
x2year  = idx_stack$year 

# Number of observations (= increments) in new data frame
N_vals   = nrow(idx_stack)

meas2x = vector(length=N_inc)
for (i in 1:N_inc) {
  print(i)
  # ID number for given tree from increment data frame
  id = incr_data$stat_id[i]
  # Year for given tree from increment data frame
  year = incr_data$year[i]

  # Find observations with given ID and year in new data frame
  meas2x[i] = which((idx_stack$tree_id == id) & (idx_stack$year == year))
}

# pdbh$year = rep(N_years, nrow(pdbh))

pdbh2val = vector(length=N_pdbh)
for (i in 1:N_pdbh){
  # ID number for given tree from DBH data frame
  id = pdbh$stat_id[i]
  # Year for given tree from DBH data frame (always the alst year in this case)
  year = pdbh$year[i]
  
  print(i)
  which((idx_stack$tree_id == id) & (idx_stack$year == year))
  
  # Find observations with given ID and year in new data frame
  pdbh2val[i] = which((idx_stack$tree_id == id) & (idx_stack$year == year))
}

# Create site-specific directories and sub-directories
site_dir <- file.path('sites',site)
if (!file.exists(site_dir)){
  dir.create(site_dir)
  dir.create(file.path(site_dir,'data'))
  dir.create(file.path(site_dir,'output'))
  dir.create(file.path(site_dir, 'figures'))
}

# Vector of plot ID for each tree (always 1 in this case)
plot_id = rep(1, N_trees)

# Save data in RDS file
saveRDS(list(N_trees=N_trees, 
           N_years=N_years,
           N_vals=N_vals,
           N_inc = N_inc,
           logXobs=logXobs, 
           logPDobs=logPDobs,
           year_idx=year_idx, 
           N_taxa=N_taxa,
           pdbh_year=pdbh_year_id,
           idx_tree =idx_tree, 
           pdbh2val=pdbh2val,
           x2tree=x2tree,
           x2year=x2year,
           meas2x = meas2x, 
           m2taxon = m2taxon,
           taxon = taxon,
           taxaMatch=taxaMatch,
           plot_id = plot_id,
           years = years,
           m2tree = m2tree,
           m2t = m2t,
           distance=distance),
           file=paste0('sites/', site, '/data/tree_data_', site ,'_STAN_', dvers, '.RDS'))
