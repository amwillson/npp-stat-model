data {
  int<lower=0> N_trees;    // number of trees 
  int<lower=0> N_years;    // number of years
  int<lower=0> N_vals;    // number of values to estimate
  int<lower=0> N_inc;     // total number of increments
  int<lower=0> pdbh2val[N_trees];    // number of trees
  int<lower=0> x2tree[N_vals];    // number of trees
  int<lower=0> x2year[N_vals];    // number of trees 
  real logPDobs[N_trees];  // DBH when increment was taken
  real logXobs[N_inc]; // increments
  real<lower=0> sig_d_obs; 
  int<lower=0> idx_tree[N_trees, 2]; // increments
  int<lower=0> year_idx[N_trees, 2]; 
  int<lower=0> meas2x[N_inc]; // index to relate measurements to estimated X vals
}
transformed data {
}
parameters {
  real beta0;
  real beta[N_trees];
  real<lower=0> beta_sd;
  real beta_t[N_years];
  real<lower=0> beta_t_sd;
  real<lower=0> sig_x;
  real<lower=0> sig_x_obs;
  //real<lower=0> sig_d_obs;
  real<lower=0> D0[N_trees];
  vector<lower=0>[N_vals] X;
  //vector<lower=0, upper=400>[N_species] sigma;
  //real<lower=0, upper=200> a;
  
}
transformed parameters {
    // process evolution
  vector<lower=0>[N_vals] D;

  for (tree in 1:N_trees){
    //print(D0[tree]);
    //print( D0[tree] + 2.0 * X[tree, year_idx[tree,1]] / 10.0);
    D[idx_tree[tree,1]] = D0[tree] + 2.0 * X[idx_tree[tree,1]] / 10.0;
    for (val in (idx_tree[tree,1]+1):(idx_tree[tree,2])){
      D[val] = D[val-1] + 2.0 * X[val] / 10.0;
    }
  }
}

model{
  
  beta0     ~ normal(0, 1.0/0.00001);
  sig_x_obs ~ uniform(0, 2.0);
  sig_d_obs ~ uniform(1e-6, 1000);
  
  sig_x     ~ uniform(1e-6, 1000);
  beta_sd   ~ uniform(1e-6, 1000);
  beta_t_sd ~ uniform(1e-6, 1000);
  
    
  for(tree in 1:N_trees) {
    D0[tree] ~ uniform(0, 50);
    beta[tree] ~ normal(beta0, beta_sd);
  }
  
  for(year in 1:N_years) {
    beta_t[year] ~ normal(0, beta_t_sd);
  }
  
  // increment likelihood
  for (val in 1:N_vals){
    //for (year in year_idx[tree, 1]:year_idx[tree,2]){
     X[val] ~ lognormal(beta[x2tree[val]] + beta_t[x2year[val]], sig_x);
   }
   
  for (inc in 1:N_inc){
   logXobs[inc] ~ normal(log(X[meas2x[inc]]), sig_x_obs);
  }
  
  for (tree in 1:N_trees){ 
    if (logPDobs[tree] == -999){
      //} else if () {
    } else {
      //logPDobs[tree] ~ student_t(3, log(D[pdbh2val[tree]-1] + 2*X[pdbh2val[tree]]/10.0), sig_d_obs);
      logPDobs[tree] ~ student_t(3, log(D[pdbh2val[tree]]), sig_d_obs);
    }
  }
}
