dvers="v0.1"
mvers="v0.1"
site = "SYLVANIA"

fname_data = paste0('tree_data_', site, '_STAN_', dvers)
fname_model = "ring_model_t_pdbh"

dat = readRDS(paste0('sites/', site, '/data/', fname_data, '.RDS'))
load(paste0('output/ring_model_t_pdbh_HMC_NOCOVAR_v3.0.Rdata'))
#col_names = sapply(strsplit(colnames(out), '\\['), function(x) x[[1]])

PDobs = exp(dat$logPDobs)
Xobs = exp(dat$logXobs)

D0 = PDobs
for (tree in 1:dat$N_trees) {
  print(paste0("Tree: ", tree))
  for (y in 1:length(dat$years)) {
    
    year = dat$years[y]
    xs = Xobs[which((dat$m2tree==tree) & (dat$m2t==year-1950))]
    if (length(xs) >= 1) {
      xs_mean = mean(xs)
    } else if (length(xs)==0) {
      xs_mean = 0
    }
    print(xs)
    D0[tree] = D0[tree] - 2*xs_mean / 10.0;
  }
}
print(D0)

hist(D0)

#hist(exp(dat$logPDobs))
