## This script provides visualization for the data products processed
## in the build_data script

## These should be compared with assumptions made in the STAN model
## to ensure that the assumptions made are appropriate given the data

## Prepare workspace

rm(list = ls())
library(ggplot2)
library(reshape2)
library(dplyr)
setwd('/Users/Aly/Google Drive 2/McLachlan/SU20/npp-stat-model/')

## Load data
data = readRDS('sites/SYLVANIA/data/tree_data_SYLVANIA_STAN_v0.1.RDS')

#########################
## Prior visualization ##
#########################

## This part of the script just sets up some vectors of values to visualize 
## the distributions from which priors are drawn in the model

beta.mu = rnorm(1000, 0, 1.0/0.0001)
sigma.x.obs = runif(1000, 1e-6, 2.0)
sigma.d.obs = runif(1000, 1e-6, 1000)
sigma.x = runif(1000, 1e-6, 1000)
beta.sigma = runif(1000, 1e-6, 1000)
beta.t.sigma = runif(1000, 1e-6, 1000)
d.i.0 = runif(1000, 0, 80)

priors = c(beta.mu, sigma.x.obs, sigma.d.obs, sigma.x, beta.sigma, beta.t.sigma, d.i.0)
priors = as.data.frame(priors)
priors$var = c(rep('beta.mu', times = 1000), rep('sigma.x.obs', times = 1000), rep('sigma.d.obs', times = 1000), rep('sigma.x', times = 1000), rep('beta.sigma', times = 1000), rep('beta.t.sigma', times = 1000), rep('d.i.0', times = 1000))

ggplot(priors, aes(x = priors, group = var, color = var)) +
  geom_density(show.legend = F) +
  facet_wrap(~var, scales = 'free')

##################
## Data visuals ##
##################

## This section just plots the data products from the build_data script

## Variables of interest:
  ## logXobs = observed increments
  ## logPDobs = observed DBH

## Most of the rest of the stuff in data helps with indexing in the STAN model

# Increment data
xobs = data$logXobs
xobs = exp(xobs)
Xobs = matrix(, nrow = data$N_trees, ncol = data$N_years)
spec = matrix(, nrow = data$N_trees, ncol = data$N_years)

# Organize increment data by tree and by year
for(i in 1:data$N_trees){
  for(j in 1:data$N_years){
    species = data$taxon[i]
    spec[i,j] = species
    x = xobs[which(data$m2tree == i & data$m2t == data$years[j] - 1900)]
    if(length(x) >= 1){
      x = mean(x)
    }
    if(length(x) == 0){
      x = 0
    }
    Xobs[i,j] = x
  }
}

colnames(Xobs) = data$years
rownames(Xobs) = 1:data$N_trees

Xobs_melt = melt(Xobs)

colnames(Xobs_melt) = c('Tree', 'Year', 'X')

# Plot all increments for each tree together
ggplot(Xobs_melt, aes(x = Year, y = X, group = Tree, color = as.factor(Tree))) +
  geom_line(show.legend = F)

# Same, but with individual tree facets
ggplot(Xobs_melt, aes(x = Year, y = X, group = Tree, color = as.factor(Tree))) +
  geom_line(show.legend = F) +
  facet_wrap(~Tree)

# Use average increments and current DBH to get initial DBH (D0)

# Set initial DBH (= current DBH because we are working backwards)
D = matrix(, nrow = data$N_trees, ncol = data$N_years)
D[,data$N_years] = exp(data$logPDobs)

for(i in 1:data$N_trees){
  for(j in (data$N_years - 1):1){
    inc = Xobs_melt$X[which(Xobs_melt$Tree == i & Xobs_melt$Year == j + 1900)]
    D[i,j] = D[i,j+1] - (2 * inc / 10)
  }
}

D0 = D[,1]

# Plot histogram of estimated initial DBH
ggplot() +
  geom_histogram(aes(x = D0), bins = 40)

# Plot density with prior distribution
ggplot() +
  stat_density(aes(x = D0, fill = 'Calc'), alpha = 0.5) +
  stat_density(aes(x = priors$priors[which(priors$var == 'd.i.0')], fill = 'Prior'), alpha = 0.5) +
  scale_fill_manual(values = c('Calc' = 'lightblue', 'Prior' = 'seagreen'), name = '')

max(D0)

# Organize DBH estimates
colnames(D) = c(data$years)
rownames(D) = c(1:data$N_trees)

D_melt = melt(D)
colnames(D_melt) = c('Tree', 'Year', 'DBH')

ggplot(D_melt, aes(x = Year, y = DBH, group = Tree, color = as.factor(Tree))) +
  geom_line(show.legend = F)

ggplot(D_melt, aes(x = Year, y = DBH, group = Tree, color = as.factor(Tree))) +
  geom_line(show.legend = F) +
  facet_wrap(~Tree)

##########################
## Visualize by species ##
##########################

## This considers the behavior of each species over time
## In case there are species-specific idiosyncracies to the data/site

# Organize species from increment processing above
colnames(spec) = data$years
rownames(spec) = 1:data$N_trees

spec_melt = melt(spec)
colnames(spec_melt) = c('Tree', 'Year', 'Species')

# Add to Xobs_melt
X_spec = Xobs_melt %>%
  full_join(spec_melt, by = c('Tree', 'Year'))

# Plot increments with species as colors
ggplot(X_spec, aes(x = Year, y = X, group = Tree, color = as.factor(Species))) +
  geom_line() +
  scale_color_discrete(name = 'Species', breaks = c('1', '2', '3', '4'), labels = c('ACSA', 'BEAL', 'THOC', 'TSCA'))

labs = c('1' = 'ACSA', '2' = 'BEAL', '3' = 'THOC', '4' = 'TSCA')

# Plot increments with species as facets
ggplot(X_spec, aes(x = Year, y = X, group = Tree, color = as.factor(Tree))) +
  geom_line(show.legend = F) +
  facet_wrap(~Species, labeller = labeller(Species = labs))

# Add DBH to data frame with species
X_D_spec = X_spec %>%
  full_join(D_melt, by = c('Tree', 'Year'))

ggplot(X_D_spec, aes(x = Year, y = DBH, group = Tree, color = as.factor(Species))) +
  geom_line() +
  scale_color_discrete(name = 'Species', breaks = c('1', '2', '3', '4'), labels = c('ACSA', 'BEAL', 'THOC', 'TSCA'))

ggplot(X_D_spec, aes(x = Year, y = DBH, group = Tree, color = as.factor(Tree))) +
  geom_line(show.legend = F) +
  facet_wrap(~Species, labeller = labeller(Species = labs))