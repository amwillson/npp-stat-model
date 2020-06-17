## Prepare workspace
rm(list = ls())
library(rstan)
setwd('/Users/Aly/Google Drive 2/McLachlan/SU20/npp-stat-model/')

# source('config_HMC')
# Data and model versions (should match those in build_data)
dvers="v0.1"
mvers="v0.1"
# Specify site (should correspond with site-specific directories)
site = "SYLVANIA"

## Before running the rest of this script, you need to consider what kind of data you have for your site
## At Sylvania, there is only tree ring data and no census data
## Therefore, we take the measurement error for DBH from a site that has both of these datasets
## This procedure is in the process of being revised
## If your site has both, you should NOT do this

## Get required data

# File name for data coming from the last step, build_data
fname_data = paste0('tree_data_', site, '_STAN_', dvers)
# File name for model (this is a STAN model, only one is included in the repo right now)
fname_model = "ring_model_t_pdbh"

# Load in RDS from build_data
dat = readRDS(paste0('sites/', site, '/data/', fname_data, '.RDS'))

# Load in data file from another site with tree ring and census data
# This will be in your environment as "out"
# This file is a matrix of parameter estimates from another site in the PalEON network
# We are only interested in extracting the DBH measurement error (sigma d obs) parameter estimates from the matrix
load('/Users/Aly/Google Drive 2/McLachlan/SU20/ring_model_t_pdbh_HMC_NOCOVAR_v3.0.Rdata')

# Make a list of column names from the model output loaded above 
# (these correspond to the paramter estimated in that column)
col_names = sapply(strsplit(colnames(out), '\\['), function(x) x[[1]])
# Visualize: notice that only a few estimates are much higher than the mode
hist(out[,which(col_names=="sig_d_obs")])
# Take mean of measurement error columns
sig_d_obs = mean(out[,which(col_names=="sig_d_obs")])

# Add to data (from RDS)
dat$sig_d_obs = sig_d_obs

#######################################################################################################################################
# full model but with zero covariance; not efficient
#######################################################################################################################################
## This part runs the STAN model found in the ring_model_t_pdbh_sigd_STAN file
## using the processed data that is the output from build_data

compiled <- stan_model(file = 'models/stan/ring_model_t_pdbh_sigd_STAN.stan')

fit <- sampling(compiled, 
                data = dat, 
                iter = 5000, 
                chains = 1,
                verbose=TRUE)

rm(compiled)

post=rstan::extract(fit)
rm(fit)

saveRDS(post, file = paste0('sites/', site, '/output/', fname_model, '_', site, '_', mvers, '.RDS'))
