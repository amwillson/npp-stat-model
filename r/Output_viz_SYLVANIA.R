## This script provides a rough visualization of the STAN model output

## Prepare workspace
rm(list = ls())
library(ggplot2)
library(reshape2)
library(dplyr)
setwd('/Users/Aly/Google Drive 2/McLachlan/SU20/npp-stat-model/')

## Load model output
out = readRDS('sites/SYLVANIA/output/ring_model_t_pdbh_SYLVANIA_v0.1.RDS')
data = readRDS('sites/SYLVANIA/data/tree_data_SYLVANIA_STAN_v0.1.RDS')

# Unlist
for(i in 1:length(out)){
  temp = out[[i]]
  eval(parse(text=paste(names(out)[[i]],"= temp")))
}
rm(out)
rm(temp)

# Reformat
beta_sd = as.vector(beta_sd)
beta_t_sd = as.vector(beta_t_sd)
beta0 = as.vector(beta0)
sig_x = as.vector(sig_x)
sig_x_obs = as.vector(sig_x_obs)

#################
## Trace plots ##
#################
# Trace plots of vectors
ggplot() +
  geom_line(aes(x = 1:2500, y = beta_sd))

ggplot() +
  geom_line(aes(x = 1:2500, y = beta_t_sd))

ggplot() +
  geom_line(aes(x = 1:2500, y = beta0))

ggplot() +
  geom_line(aes(x = 1:2500, y = sig_x))

ggplot() +
  geom_line(aes(x = 1:2500, y = sig_x_obs))

# Trace plots of matrices
colnames(beta) = 1:63
rownames(beta) = 1:2500

beta_melt = melt(beta)
colnames(beta_melt) = c('Iteration', 'Tree', 'Val')

# Plot trace plot by tree
ggplot(beta_melt, aes(x = Iteration, y = Val, group = Tree, color = as.factor(Tree))) +
  geom_line(show.legend = F)

# Extract species from data to include in facet wrapping
spec = matrix(, nrow = data$N_trees, ncol = data$N_years)
for(i in 1:data$N_trees){
  for(j in 1:data$N_years){
    spec[i,j] = data$taxon[i]
  }
}

colnames(spec) = data$years
rownames(spec) = 1:data$N_trees

spec_melt = melt(spec)
colnames(spec_melt) = c('Tree', 'Year', 'Species')

# Add species to beta data frame
beta_melt = beta_melt %>%
  full_join(spec_melt, by = 'Tree')

labs = c('1' = 'ACSA', '2' = 'BEAL', '3' = 'THOC', '4' = 'TSCA')

# These are commented out because they are very slow
# They do work and can  be viewed

# Plot by species
#ggplot(beta_melt, aes(x = Iteration, y = Val, group = Tree, color = as.factor(Tree))) +
#  geom_line(show.legend = F) +
#  facet_wrap(~Species, labeller = labeller(Species = labs))

# Plot by species using color instead of facets
#ggplot(beta_melt, aes(x = Iteration, y = Val, group = Tree, color = as.factor(Species))) +
#  geom_line() +
#  scale_color_discrete(name = 'Species', breaks = c('1', '2', '3', '4'), labels = c('ACSA', 'BEAL', 'THOC', 'TSCA'))

# Repeat for beta t
colnames(beta_t) = data$years
rownames(beta_t) = 1:2500

beta_t_melt = melt(beta_t)
colnames(beta_t_melt) = c('Iteration', 'Year', 'Val')

ggplot(beta_t_melt, aes(x = Iteration, y = Val, group = Year, color = as.factor(Year))) +
  geom_line(show.legend = F)

# Repeat for D0
colnames(D0) = 1:data$N_trees
rownames(D0) = 1:2500

D0_melt = melt(D0)
colnames(D0_melt) = c('Iteration', 'Tree', 'Val')

D0_melt = D0_melt %>%
  full_join(spec_melt, by = 'Tree')

ggplot(D0_melt, aes(x = Iteration, y = Val, group = Tree, color = as.factor(Tree))) +
  geom_line(show.legend = F)

###################
## Density plots ##
###################

ggplot() +
  stat_density(aes(x = beta_sd))

ggplot() +
  stat_density(aes(x = beta_t_sd))

ggplot() +
  stat_density(aes(x = beta0))

ggplot() +
  stat_density(aes(x = sig_x))

ggplot() +
  stat_density(aes(x = sig_x_obs))

ggplot(beta_melt, aes(x = Val, fill = as.factor(Species))) +
  stat_density() +
  facet_wrap(~Tree)

ggplot(beta_t_melt, aes(x = Val)) +
  stat_density() +
  facet_wrap(~Year)

ggplot(D0_melt, aes(x = Val, fill = as.factor(Species))) +
  stat_density() +
  facet_wrap(~Tree)
