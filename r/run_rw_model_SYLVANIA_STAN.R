library(rstan)

# source('config_HMC')
dvers="v0.1"
mvers="v0.1"
site = "SYLVANIA"

fname_data = paste0('tree_data_', site, '_STAN_', dvers)
fname_model = "ring_model_t_pdbh"

dat = readRDS(paste0('sites/', site, '/data/', fname_data, '.RDS'))

load(paste0('output/ring_model_t_pdbh_HMC_NOCOVAR_v3.0.Rdata'))

col_names = sapply(strsplit(colnames(out), '\\['), function(x) x[[1]])
hist(out[,which(col_names=="sig_d_obs")])
sig_d_obs = mean(out[,which(col_names=="sig_d_obs")])

dat$sig_d_obs = sig_d_obs

#######################################################################################################################################
# full model but with zero covariance; not efficient
#######################################################################################################################################

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