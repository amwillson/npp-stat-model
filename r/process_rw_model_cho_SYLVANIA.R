# THIS SCRIPT ASSUMES NO CENSUS DATA AVAILABLE # 

rm(list=ls())
# setwd('~/Desktop/npp-stat-model')

library(ggplot2)
library(reshape2)
library(abind)
library(dplyr)

## TO DO ##
# config file 
# make file names easier to find/read or just consistent across sites
# settle on directory set-up
# determine what year ring widths were collected so we can remove those values 

## QUESTIONS FOR ANDRIA ##
# how many iterations should we keep?

## IMPORTANT NOTES ##
# ring widths collected around 2013; incomplete ring widths and most recent year with incomplete information => thrown out 

#############################################################################################################
### A. Set up working environment and directory for script ##################################################
#############################################################################################################

########################################
# Variables to be adjusted for each site
########################################

# date for save purposes
date = '3June2020' 

# set up output 
site <- 'SYLVANIA'

# data version
dvers = 'v0.1'
#mvers = 'v4.0'

# burn-in range
burn = 1000

# load model output data for site
# data_dir = file.path('data','sites',site)
data_dir = file.path('sites',site, 'data')
dat = readRDS(paste0(data_dir,'/tree_data_', site, '_STAN_v0.1.RDS'))
# dat = readRDS(paste0(data_dir,'/NPP_STAT_MODEL_ROOSTER.RDS'))
fnames = c(paste0('ring_model_t_pdbh_', site, '_v0.1.Rdata'))
models = c('Model RW')

# are dbh and ab files already created? do we want to recreate them? 
run_ringwidth = TRUE

# the year ring widths were collected (this year and next years are removed)
rw_year = 2018

########################################

output_dir <- file.path('sites',site)
if (!file.exists(output_dir)){
  dir.create(output_dir)
  dir.create(file.path(output_dir,'output'))
  dir.create(file.path(output_dir, 'figures'))
}

# get ring width data 
output_dir = file.path('sites',site, 'output')
nmodels = length(fnames)
post = list()
for (i in 1:length(fnames)) {
  fname_model = fnames[i]
  load(file   = paste0(output_dir,'/', fname_model))
  post[[i]]   = post
}  
niter =length(post$lp__)

# how many iterations to keep?
keep = 250 #niter - burn

# extract and rename required variables
N_trees = dat$N_trees
N_years = dat$N_years
logXobs = dat$logXobs
logPDobs = dat$logPDobs
idx_tree = dat$idx_tree
x2tree = dat$x2tree
x2year = dat$x2year
taxaMatch = dat$taxaMatch
years = dat$years
plot_id = dat$plot_id
taxon = dat$taxon
m2tree = dat$m2tree
m2t = dat$m2t
pdbh_year = dat$pdbh_year
trees = seq(1, N_trees)

# match species acronyms to level3a available species/pfts 
acro_level3a = read.csv('data/acronym_to_level3a_v0.2.csv', stringsAsFactors = FALSE)
acro_level3a = left_join(taxaMatch, acro_level3a, by = c('species'='acronym'))
choj = read.csv('data/level3a_to_chojnacky_v0.4.csv', stringsAsFactors = FALSE)
choj = left_join(acro_level3a, choj, id = 'level3a') 

array_or_vector_p <- function(object) {
  ## Test if object is what we consider an array for the purposes of the code below.
  
  ## Lists are not arrays, but vectors are, even if they don't have a dim attribute. The following give an error:
  ## - objects which are neither lists, (non-list) vectors, or arrays
  ## - vectors and arrays with zero length
  ## - arrays with rank lower than 1
  if (is.list(object))
    FALSE
  else {
    stopifnot(is.array(object) || is.vector(object))
    stopifnot(length(object) > 0)
    if (is.array(object))
      stopifnot(length(dim(object)) > 0)
    TRUE
  }
}
subarray <- function(array, index) {
  ## equivalent to array[index, , ...]
  stopifnot(array_or_vector_p(array))
  if (is.vector(array))
    array[index]
  else
    do.call("[",c(list(array,index),rep(TRUE,length(dim(array))-1)))
}
thin_subarrays <- function(subarrays, n=30) {
  ## thin subarrays. for a single n, draw that many, otherwise use that as indexes.
  indexes <- if (length(n)==1) 1:n else n
  lapply(subarrays, function(array) subarray(array, indexes))
}

post[[1]] = thin_subarrays(post[[1]], n=keep)

#############################################################################################################
### B. Process diameter and ab predictions for ring width data ##############################################
#############################################################################################################

##### Step 1: Estimate diameter and aboveground biomass for all trees from ring width data #####

# this function takes the stat model results (dbh, rw) and returns a list of dbh arrays (tree X years X iterations)
# for each pft in the data, a list of the tree indices in the order they were used in the list, and the order of the
# pfts within the list 
build_dbh_p <- function(out, rw2year, rw2tree, Nyears, tree_list, keep, taxaMatch, taxon_list){
  
  # extract colnames and number of iterations (different value types in data, want diameter)
  col_names = names(out)
  niter   = length(out[[1]])
  
  # get species from taxon list 
  pft = as.vector(taxaMatch[match(taxon_list, taxaMatch$number),'species'])
  pfts  = sort(unique(pft))
  N_pft = length(pfts) 
  
  # create list for dbh data for each pft 
  dbh_p = list(length=N_pft)
  tree_order = vector(length=0)
  
  # loop through pfts
  for (p in 1:N_pft){
    
    # identify all trees of this pft and put in list to track 
    pft_trees = which(pft == pfts[p])
    tree_order = c(tree_order, tree_list[pft_trees])
    
    # create save array to put dbh information for pft in 
    # dbh saved for keep iterations for every year and every tree of each pft 
    N_pft_trees = length(pft_trees)
    dbh_p[[p]] = array(NA, c(N_pft_trees, Nyears, keep))
    
    # loop through trees of pft
    for (i in 1:N_pft_trees) {
      print(paste0('Tree ', i))
      tree = tree_list[pft_trees[i]]
      
      # connects individual tree to place in out and which years there is data for the tree 
      tree_idx   = which(rw2tree == tree)
      tree_years = rw2year[tree_idx]
      
      # extracts all diameter measurements for tree in iterations to keep 
      mu_dbh = out$D[(niter-keep+1):(niter),tree_idx]
      
      # adjust negative dbh values
      if (any(mu_dbh < 0)) {
        mu_dbh[mu_dbh < 0] = 0
      }
      
      # saves extracted dbh measurements into correct years for the pft
      dbh_p[[p]][i,tree_years,] = t(mu_dbh)
    }
  }
  return(list(dbh_p=dbh_p, tree_order=tree_order, pft_order = pfts))
}

# this function takes the built dbh array list and returns a biomass array (tree X years X iterations) where
# biomass is calculated using Chojnacky's allometric coefficients in the biomass equation
build_ab_choj <- function(dbh_p, Ntrees, Nyears, keep, choj_pfts){
  
  # build save array 
  ab_p = array(NA, c(Ntrees, Nyears, keep))
  
  # loop through iterations to keep 
  for (k in 1:keep){
    print(paste0("Iteration ", k))
    
    # get dbhs for all trees for all years in iteration (number of trees X years)
    dbh_list = dbh_p[,,k]
    
    # calculate the biomass contribution of each tree for each year in iteration
    # ln(biomass-kg) = beta0 + beta1 * ln(diameter-cm)
    b0 = choj$beta0[choj_pfts]
    b1 = choj$beta1[choj_pfts]
    ab_p[,,k] = exp(b0 + b1 * log(dbh_list))
  }
  
  return(ab_p=ab_p)
}

# if ab and diameter estimates haven't already been found 
if (run_ringwidth){
  
  # for first model with census 
  # obtain diameter estimates for all trees for
  dbh_1 = build_dbh_p(post[[1]], x2year, x2tree, N_years, trees, keep, taxaMatch, taxon)
  dbh_p_1 = dbh_1$dbh_p
  tree_order_1 = dbh_1$tree_order
  pft_order_1  = dbh_1$pft_order
  
  # make list of arrays into one large array across species (all trees X years X iterations)
  dbh_p_1_org = abind(dbh_p_1, along=1)
  dbh_p_1 = dbh_p_1_org[sort(tree_order_1, index.return=TRUE)$ix,,]
  
  # get choj pft index for coefficients in matrix 
  choj_pfts_1 = taxon[tree_order_1]
  
  # find biomass for all trees, years, and iterations
  ab_p_1 = build_ab_choj(dbh_p_1, N_trees, N_years, keep, choj_pfts_1)
  
  # save files 
  #saveRDS(dbh_p_1, paste0('output/',site,'/dbh_p_1_', site, '_', mvers, '.rds'))
  #saveRDS(ab_p_1, paste0('output/',site,'/ab_p_1_', site, '_', mvers, '.rds'))
  
}

# added separate save variable in case we need to go back to this step
dbh_p_1_raw = dbh_p_1
ab_p_1_raw = ab_p_1

##### Step 2: Remove the data points for trees with diameters less than 5 cm #####

# remove trees that have a mean diameter less than 5 centimeters across all iterations for a given year for both models
for (tree in 1:N_trees){
  for (year in 1:N_years){
    dbh_mean = mean(dbh_p_1[tree, year, ], na.rm=TRUE)
    if (is.na(dbh_mean)){next}
    if (dbh_mean < 5){
      #print(tree)
      dbh_p_1[tree, year, ] = rep(NA, keep)
      ab_p_1[tree, year, ] = rep(NA, keep)
    }
  }
}

# added save variable in case we need to back to this step
dbh_p_1_red = dbh_p_1
ab_p_1_red = ab_p_1

##### Step 3: Melt data matrices into data frames for easier use #####

# melt down to data frames
rownames(dbh_p_1) <- seq(1, N_trees)
dbh_melt = melt(dbh_p_1)
colnames(dbh_melt) = c("tree", "year", "iter", "dbh")
dbh_melt$year = years[dbh_melt$year]
dbh_melt$plot = plot_id[dbh_melt$tree]

rownames(ab_p_1) <- seq(1, N_trees)
ab_melt = melt(ab_p_1)
colnames(ab_melt) = c("tree", "year", "iter", "ab")
ab_melt$year = years[ab_melt$year]
ab_melt$plot = plot_id[ab_melt$tree]

# put models together
dbh_melt_all = data.frame(dbh_melt, model=rep("Model RW", nrow(dbh_melt)))
ab_melt_all = data.frame(ab_melt, model=rep("Model RW", nrow(ab_melt)))

dbh_melt_all$taxon = taxon[dbh_melt_all$tree]
dbh_melt_all$taxon   = taxaMatch[match(dbh_melt_all$taxon, taxaMatch$number), 'species']

# XXX: don't think we need these lines
# # aggregate for some plots later
# dbh_mean = aggregate(dbh~tree+year+plot, dbh_melt, mean, na.rm=TRUE) 
# dbh_sum = aggregate(dbh~year+plot, dbh_mean, sum, na.rm=TRUE)
# 
# ab_mean = aggregate(ab~tree+plot+year+model, ab_melt_all, mean) 
# ab_sum = aggregate(ab~year+model+plot, ab_mean, sum)

##### Step 4: Calculate biomass and biomass increment #####

# ab_small: based on trees from 10 to 20 cm DBH from inner most plot of radius 13 m
pdbh = exp(logPDobs)
distance = dat$distance
idx_small  = which(pdbh<20)
idx_medium = which((pdbh>=20)&(pdbh<30))
idx_large  = which(pdbh>=30)

# in Kg/plot, rescale so it is Mg/ha
# but plot sizes differ depending on tree size due to sampling design
ab_pr1 = ab_p_1
ab_pr1[idx_small,,] = ab_p_1[idx_small,,]/(13^2*pi) * 10
ab_pr1[idx_medium,,] = ab_p_1[idx_medium,,]/(20^2*pi) * 10
ab_pr1[idx_large,,] = ab_p_1[idx_large,,]/(30^2*pi) * 10
# ab_pr1 = ab_p_1/(30^2*pi) * 10

# melt predicted biomass for both models
ab_p1_melt = melt(ab_pr1)
colnames(ab_p1_melt) = c('tree_id', 'year_id', 'iter', 'ab')
ab_p1_melt$plot = plot_id[ab_p1_melt$tree_id]

# create total 
ab_p_melt = data.frame(ab_p1_melt, model=rep('Model RW'))
ab_p_melt$taxon   = taxon[ab_p_melt$tree_id]
ab_p_melt$taxon   = taxaMatch[match(ab_p_melt$taxon, taxaMatch$number), 'species']
ab_p_melt$year = years[ab_p_melt$year_id]

msd_p1 = dbh_p_1[distance<=13,,]
msd_p1_melt = melt(msd_p1)
colnames(msd_p1_melt) = c('tree_id', 'year_id', 'iter', 'msd')
msd_p1_melt$plot = plot_id[msd_p1_melt$tree_id]
msd_p1_melt$year = years[msd_p1_melt$year_id]

msd_p = aggregate(msd~year + plot + iter, msd_p1_melt, mean, na.rm=TRUE)

# save aboveground biomass predictions
#saveRDS(ab_p_melt, file=paste0('output/',site,'/ab_p_', site, '_', mvers, '.RDS'))

# determine biomass increment
abi_p1 = apply(ab_pr1, c(1,3), function(x) diff(x))
abi_p1 = aperm(abi_p1, c(2, 1, 3))
abi_p1_melt = melt(abi_p1)
colnames(abi_p1_melt) = c('tree_id', 'year_id', 'iter', 'abi')

# put them together
abi_p_melt = data.frame(abi_p1_melt, model=rep('Model RW'))
abi_p_melt$taxon   = taxon[abi_p_melt$tree_id]
abi_p_melt$taxon   = taxaMatch[match(abi_p_melt$taxon, taxaMatch$number), 'species']
abi_p_melt$plot = plot_id[abi_p_melt$tree_id]
abi_p_melt$year    = years[abi_p_melt$year_id] 

# save aboveground biomass increment predictions
#saveRDS(abi_p_melt, file=paste0('output/',site,'/abi_p_', site, '_', mvers, '.RDS'))

##### Step 5: Prep AGB and ABI for saving #####

abi_p_melt = abi_p_melt[,c('tree_id', 'year', 'plot', 'taxon', 'model', 'iter', 'abi')]
ab_p_melt = ab_p_melt[,c('tree_id', 'year', 'plot', 'taxon', 'model', 'iter', 'ab')]

ab_p_melt = ab_p_melt[ab_p_melt$year<rw_year,]
abi_p_melt = abi_p_melt[abi_p_melt$year<rw_year,]

colnames(ab_p_melt)[which(colnames(ab_p_melt) == "ab")] = 'value'
colnames(abi_p_melt)[which(colnames(abi_p_melt) == "abi")] = 'value'

#############################################################################################################
### C. Process dbh and ab from most recent dbh measurement ##################################################
#############################################################################################################

##### Step 1: Track diameter back from most recent dbh measurement using increment data #####

dbh_m = array(NA, c(N_trees, N_years))

# loop through all trees 
for (i in 1:N_trees) {
  print(i)
  tree = trees[i]
  tree_idx = which(trees == tree)
  pft  = as.vector(taxaMatch[which(taxaMatch$number == taxon[tree]),'species'])
  
  dbh_year = pdbh_year[i]
  dbh_tree = exp(logPDobs[tree])
  
  # for now use average increments
  incr_tree  = exp(logXobs[which(m2tree == tree)])
  incr_years = years[m2t[which(m2tree == tree)]]
  
  incr_mean = aggregate(incr_tree~incr_years, FUN=median)
  
  if (years[dbh_year]==max(incr_mean$incr_years)){
    incr_tree = incr_mean[,2]
    incr_years = incr_mean[,1]
    
    # incr_years = incr_years[-(length(incr_years))]
    # incr_tree = incr_tree[-(length(incr_tree))]
    
    tree_years = seq(min(incr_years), max(incr_years))
    
    incr_cumsum = rev(cumsum(rev(incr_tree[incr_years %in% tree_years])))
    dbh_calc = c(dbh_tree - 2*incr_cumsum/10, dbh_tree)
    
    year_idx = match(tree_years, years)
    
    if (length(dbh_calc) > length(year_idx)){
      dbh_calc = dbh_calc[-1]
    }
    
    dbh_m[tree_idx,year_idx] = dbh_calc
  } else {
    print("Wonky ring width measurement after DBH measurement!")
  }
}

##### Step 2: Get biomass using Chojnacky coefficients for all trees #####

# ln(biomass-kg) = beta0 + beta1 * ln(diameter-cm)
b0 = choj$beta0[taxon]
b1 = choj$beta1[taxon]
ab_m = exp(b0 + b1*log(dbh_m))

# gets rid of data for trees with diameter less than 5 cm
for (n in 1:N_trees){
  for (t in 1:N_years){
    if (is.na(dbh_m[n,t])){
      next
    } else if (dbh_m[n,t] < 5) {
      dbh_m[n,t] = NA
      ab_m[n,t]  = NA
    }  
  }
}

##### Step 3: Melt down and plot #####

# melt to data frames
rownames(dbh_m) <- seq(1, N_trees)
dbh_m_melt = melt(dbh_m)
colnames(dbh_m_melt) = c("tree", "year", "dbh")
dbh_m_melt$year = years[dbh_m_melt$year]
dbh_m_melt$plot = plot_id[dbh_m_melt$tree]

# in Kg/plot, rescale so it is Mg/ha
# ab_mr = ab_m/(20^2*pi) * 10
ab_mr = ab_m
ab_mr[idx_small,] = ab_m[idx_small,]/(13^2*pi) * 10
ab_mr[idx_medium,] = ab_m[idx_medium,]/(20^2*pi) * 10
ab_mr[idx_large,] = ab_m[idx_large,]/(30^2*pi) * 10
# ab_mr = ab_m/(30^2*pi) * 10

# melt measured
ab_m_melt = melt(ab_mr)
colnames(ab_m_melt) = c('tree_id', 'year_id', 'ab')
ab_m_melt$taxon   = taxon[ab_m_melt$tree_id]
ab_m_melt$taxon   = taxaMatch[match(ab_m_melt$taxon, taxaMatch$number), 'species']
ab_m_melt$plot = plot_id[ab_m_melt$tree_id]
ab_m_melt$year    = years[ab_m_melt$year_id]
ab_m_melt = ab_m_melt[ab_m_melt$year<rw_year,]

# measured biomass increment
abi_m = t(apply(ab_mr, 1, function(x) diff(x)))

# melt measured
abi_m_melt = melt(abi_m)
colnames(abi_m_melt) = c('tree_id', 'year_id', 'abi')
abi_m_melt$taxon   = taxon[abi_m_melt$tree_id]
abi_m_melt$taxon   = taxaMatch[match(abi_m_melt$taxon, taxaMatch$number), 'species']
abi_m_melt$plot = plot_id[abi_m_melt$tree_id]
abi_m_melt$year    = years[abi_m_melt$year_id]
abi_m_melt = abi_m_melt[abi_m_melt$year<rw_year,]

# note dbh not reweighted for sampling
dbh_mean = aggregate(dbh~tree+year+plot, dbh_melt, mean, na.rm=TRUE)
dbh_sum = aggregate(dbh~year+plot, dbh_mean, sum, na.rm=TRUE)
dbh_m_sum = aggregate(dbh~year+plot, dbh_m_melt, sum)

ggplot() + 
  geom_line(data=dbh_sum, aes(x=year, y=dbh, colour = 'Model')) +
  geom_line(data=dbh_m_sum, aes(x=year, y=dbh, colour = 'Empirical')) + 
  labs(title ='Sum of DBH Values for DBH and RW Obs', colour = 'Plot ID') +
  scale_color_manual(values = c('Model' = 'darkgreen', 'Empirical' = 'lightgreen'), name = '')
ggsave(paste0('sites/', site, '/figures/sum_dbh_by_plot_', site,'.pdf'))

################ ABCXYZ
#############################################################################################################
### D. Save WIKI RDS and CSV file + clean up working directory ##############################################
#############################################################################################################

# clean up working directory because the next file is huge
variables = ls()
not_rm = variables %in% c('ab_p_melt','abi_p_melt','output_dir','site',
           'ab_m_melt','abi_m_melt',
           'years', 'ab_sum', 'dbh_melt_all', 'msd_p')
rm = variables[!not_rm]
rm(list=rm)

# get rid of NA values and save to one large data frame
ab_p_melt = ab_p_melt %>% filter(!is.na(value))
abi_p_melt = abi_p_melt %>% filter(!is.na(value))
preds = rbind(data.frame(ab_p_melt, type="AB"), data.frame(abi_p_melt, type="ABI"))

# save file for WIKI
write.csv(preds, file=paste0('sites/', site, '/output/', 'NPP_STAT_MODEL_', site, '.csv'), row.names=FALSE)
saveRDS(preds, paste0('sites/', site, '/output/','NPP_STAT_MODEL_',site,'.RDS'))


dbh_preds = dbh_melt_all[,c('tree', 'year', 'plot', 'taxon', 'model', 'iter', 'dbh')]
dbh_preds = data.frame(dbh_preds, type=rep('DBH', nrow(dbh_preds)))
colnames(dbh_preds) = c('tree_id', 'year', 'plot', 'taxon', 'model', 'iter', 'value', 'type')

dbh_preds = dbh_preds %>% filter(!is.na(value))

# save file for WIKI
write.csv(dbh_preds, file=paste0('sites/', site, '/output/', 'DBH_STAT_MODEL_', site, '.csv'), row.names=FALSE)
saveRDS(dbh_preds, paste0('sites/', site, '/output/','DBH_STAT_MODEL_',site,'.RDS'))

write.csv(msd_p, file=paste0('sites/', site, '/output/','MSD_STAT_MODEL_', site, '.csv'), row.names=FALSE)
saveRDS(msd_p, paste0('sites/', site, '/output/','MSD_STAT_MODEL_', site, '.RDS'))
#############################################################################################################
### E. Data Visualizations ##################################################################################
#############################################################################################################

###### Data manipulations for plotting #####

# for aboveground biomass iterations for stat model data
ab_p_sum_by_iter = aggregate(value ~ year+iter+plot+model, ab_p_melt, function(x) sum(x, na.rm=TRUE))
ab_p_quants = aggregate(value~year+plot+model, data=ab_p_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
ab_p_quants = data.frame(ab_p_quants)
ab_p_quants = cbind(ab_p_quants[,c(1:3)], ab_p_quants[,4])
colnames(ab_p_quants) = c('year','plot','model','ab25', 'ab50', 'ab975') 

# for aboveground biomass measurements 
ab_m_sum_by_iter = aggregate(ab ~ year+plot, ab_m_melt, function(x) sum(x, na.rm=TRUE))
ab_m_quants = aggregate(ab~year+plot, data=ab_m_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
ab_m_quants = data.frame(ab_m_quants)
ab_m_quants = cbind(ab_m_quants[,c(1:2)], ab_m_quants[,3])
colnames(ab_m_quants) = c('year','plot','ab25', 'ab50', 'ab975') 

# for aboveground biomass increment estimates for stat model results
abi_p_sum_by_iter = aggregate(value ~ year+plot+iter+model, abi_p_melt, function(x) sum(x, na.rm=TRUE))
abi_p_quants = aggregate(value~year+plot+model, data=abi_p_sum_by_iter, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
abi_p_quants = data.frame(abi_p_quants)
abi_p_quants = cbind(abi_p_quants[,1:3], abi_p_quants[,4])
colnames(abi_p_quants) = c('year','plot','model','ab25', 'ab50', 'ab975')  

# sums for measured
abi_m_sum = aggregate(abi ~ year+plot, abi_m_melt, function(x) sum(x, na.rm=TRUE))

# MK :: super hack. apologies...
#abi_m_sum$year = abi_m_sum$year-1

##### Plot 1: AGB plots for different model and empirical ######

cols = c('Model RW'='#8c2d04', 'Empirical RW'='#238b45')
cols_fill = c('Model RW'="#fdd0a2", 'Empirical RW'="white")

# same y scales
ggplot() +  
  geom_ribbon(data=ab_p_quants, aes(x=year, ymin=ab25, ymax=ab975, fill=model, colour=model), alpha=0.4) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab50, colour=model), size=1) + 
  geom_ribbon(data=ab_m_quants, aes(x=year, ymin=ab25, ymax=ab975, colour='Empirical RW', fill='Empirical RW'), alpha=0.4) +
  geom_line(data=ab_m_quants, aes(x=year, y=ab50, colour='Empirical RW'), size=1) + 
  facet_grid(plot~.) + 
  scale_color_manual(values=cols, name='Method') +
  scale_fill_manual(values=cols_fill, name='Method') +
  scale_x_continuous(breaks=seq(min(years), max(years), by=5)) +
  theme_bw() + 
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14)) +
  ylab("Biomass (Mg/ha)") + 
  xlab('Year')
ggsave(file=paste0('sites/', site, '/figures/AGB_by_plot_',site,'.png'), type="cairo")

# free y scales 
ggplot() +  
  geom_ribbon(data=ab_p_quants, aes(x=year, ymin=ab25, ymax=ab975, color = model, fill = model), alpha=0.4) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab50, colour=model), size=1) + 
  geom_ribbon(data=ab_m_quants, aes(x=year, ymin=ab25, ymax=ab975, colour='Empirical RW', fill='Empirical RW'), alpha=0.4) +
  geom_line(data=ab_m_quants, aes(x=year, y=ab50, colour='Empirical RW'), size=1) + 
  geom_line(data=ab_p_quants, aes(x=year, y=ab25, colour=model), linetype=1, size=0.5) +
  geom_line(data=ab_p_quants, aes(x=year, y=ab975, colour=model), linetype=1, size=0.5) +
  facet_grid(plot~., scales="free_y") + scale_color_manual(values=cols, name='Method') +
  scale_fill_manual(values=cols_fill, name='Method') +
  theme_bw() + 
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14)) +
  ylab("Biomass (Mg/ha)") + xlab('Year') +
  scale_x_continuous(breaks=seq(min(years), max(years), by=5))

ggsave(file=paste0('sites/', site,'/figures/AGB_by_plot_',site,'_freey.png'), type="cairo")

##### Plot 2: AGB increment plots for model and empirical #####
ggplot() +  
  geom_ribbon(data=abi_p_quants, aes(x=year, ymin=ab25, ymax=ab975, fill=model, colour=model), alpha=0.4) +
  geom_line(data=abi_p_quants, aes(x=year, y=ab50, colour=model), size=1) + 
  geom_line(data=abi_m_sum, aes(x=year, y=abi, colour='Empirical RW', fill='Empirical RW'), size=1) + 
  facet_grid(plot~.) +
  scale_color_manual(values=cols, name='Method') +
  scale_fill_manual(values=cols_fill, name='Method') + 
  theme_bw() +
  ylab("Biomass Increment (Mg/ha/year)") + xlab('Year') +
  scale_x_continuous(breaks=seq(min(years), max(years), by=5), limits=c(min(years), max(years)))
ggsave(file=paste0('sites/', site,'/figures/AGBI_by_plot_',site,'.png'), type="cairo")

