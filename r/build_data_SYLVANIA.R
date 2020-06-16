#!/usr/bin/Rscript
library(plotrix)
library(dplR)
library(fields)
library(reshape2)
library(plyr)

# source("config_HMC")
dataDir = '~/Documents/projects/npp-stat-model/data'
site = 'SYLVANIA'
dvers = "v0.1"
mvers = "v0.1"

treeMeta = read.csv("data/SYLVANIA/alexander_sylvania_June2018.csv")

rwFiles <- list.files(paste0("data/", site, '/', "rwl"))
rwFiles <- rwFiles[grep(".rwl$", rwFiles)]
rwData <- list()
for(fn in rwFiles) {
  id <- gsub(".rw", "", fn)
  rwData[[id]] <- t(read.tucson(file.path("data", site, "rwl", fn)))  # rows are tree, cols are times
}

incr = ldply(rwData, rbind)
incr = incr[,c(".id", sort(colnames(incr)[2:ncol(incr)]))]
rownames(incr) = as.vector(unlist(lapply(rwData, rownames)))
incr[,1] = rownames(incr)


incr = incr[, which((as.numeric(colnames(incr))>1950)|(colnames(incr)==".id"))]

######################################################################################################################################
## make STAN data
######################################################################################################################################

incr_data = melt(incr)
colnames(incr_data) = c('id', 'year', 'incr')
incr_data$plot   = rep(1, nrow(incr_data))#as.numeric(substr(incr_data$id, 3, 3))
incr_data$TreeID = incr_data$id 
incr_data$id     = as.numeric(substr(incr_data$id, 4, 6))
incr_data$year   = as.vector(incr_data$year)

tree_ids = unique(substr(incr_data$TreeID, 1, 6))
N_trees  = length(tree_ids)
stat_ids = seq(1, N_trees)

incr_data$stat_id = stat_ids[match(substr(incr_data$TreeID, 1, 6), tree_ids)]
for (n in 1:nrow(incr_data)){
  print(n)
  incr_data$taxon[n] = as.vector(treeMeta$species[which(as.numeric(substr(treeMeta$ID, 4, 6)) == incr_data$id[n])])
}

N_taxa = length(unique(incr_data$taxon))
taxaMatch = data.frame(species=sort(unique(incr_data$taxon)), number=seq(1, N_taxa))

taxon = aggregate(taxon~stat_id, incr_data, unique)[,2]
taxon = taxaMatch$number[match(taxon, taxaMatch$species)]

plot_id = aggregate(plot~stat_id, incr_data, unique)[,2]

year_end = max(as.numeric(incr_data$year), na.rm=TRUE)
year_start = min(as.numeric(incr_data$year), na.rm=TRUE)
N_years = year_end - year_start + 1
years = seq(year_start, year_end)

# order by tree and year
incr_data = incr_data[order(incr_data$stat_id, incr_data$year),]
incr_data = incr_data[which(!is.na(incr_data$incr)),]
incr_data$year = as.numeric(incr_data$year) - year_start + 1
N_inc   = nrow(incr_data) # number of measurements
m2t     = incr_data$year
m2tree  = incr_data$stat_id
m2plot  = incr_data$plot
m2taxon = taxaMatch$number[match(incr_data$taxon, taxaMatch$species)]
Xobs    = incr_data$incr
Xobs[Xobs==0] = 0.0001
logXobs = log(Xobs)

year_idx = data.frame(year_start=as.numeric(aggregate(year~stat_id, data=incr_data, FUN=min, na.rm=TRUE)[,2]), 
                      year_end=as.numeric(aggregate(year~stat_id, incr_data, max)[,2]))

# make pdbh
pdbh = aggregate(year~stat_id+plot+id, incr_data, max, na.rm=TRUE)
pdbh = pdbh[order(pdbh$stat_id),]

for (n in 1:nrow(pdbh)){
  pdbh$dbh[n] = treeMeta$dbh[which(as.numeric(substr(treeMeta$ID, 4, 6)) == pdbh$id[n])]
  pdbh$distance[n] = treeMeta$distance[which(as.numeric(substr(treeMeta$ID, 4, 6)) == pdbh$id[n])]
}
N_pdbh = nrow(pdbh)
logPDobs = log(pdbh$dbh)
pdbh_tree_id = pdbh$stat_id
pdbh_year_id = pdbh$year
distance = pdbh$distance

idx_stack = data.frame(meas=numeric(0), tree_id=numeric(0), year=numeric(0))
n = 1
for (tree in 1:N_trees){
  year = seq(year_idx[tree,1], year_idx[tree,2])
  meas = seq(n, n+length(year)-1)
  n = n + length(year)
  idx_stack = rbind(idx_stack, data.frame(meas=meas, tree_id=rep(tree, length(year)), year=year))
}

idx_tree = which(!duplicated(idx_stack$tree_id))
idx_tree = data.frame(idx_tree, c(idx_tree[-1]-1, nrow(idx_stack)))

x2tree  = idx_stack$tree_id
x2year  = idx_stack$year 

N_vals   = nrow(idx_stack)

meas2x = vector(length=N_inc)
for (i in 1:N_inc) {
  print(i)
  id = incr_data$stat_id[i]
  year = incr_data$year[i]
  
  meas2x[i] = which((idx_stack$tree_id == id) & (idx_stack$year == year))
}

pdbh2val = vector(length=N_pdbh)
for (i in 1:N_pdbh){
  id = pdbh$stat_id[i]
  year = pdbh$year[i]
  
  print(i)
  which((idx_stack$tree_id == id) & (idx_stack$year == year))
  
  pdbh2val[i] = which((idx_stack$tree_id == id) & (idx_stack$year == year))
}

site_dir <- file.path('sites',site)
if (!file.exists(site_dir)){
  dir.create(site_dir)
  dir.create(file.path(site_dir,'data'))
  dir.create(file.path(site_dir,'output'))
  dir.create(file.path(site_dir, 'figures'))
}

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