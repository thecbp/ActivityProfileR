# Supplement to "Nonparametric Comparisons of Activity Profiles From Wearable Device Data" by Chang and McKeague
# The following code illustrates how we obtain results of the proposed EL test,
# how the empirical rejection rates are obtained, and how we generate data based on the Ornstein--Uhlenbeck process
############# required packages
# For cluster computing:
library(parallel)
library(foreach)
library(doParallel)
# For root solving
library(rootSolve)
library(plyr)
# For Ornstein--Uhlenbeck process generation
library(ESGtoolkit)
############# END required packages

############# User input
# Where R code containing the functions is
dir_path = paste("C:\\Users\\", Sys.getenv("USERNAME"), "\\Dropbox\\papers_mine\\paper_fnl_data", sep = "") 
# Where the results will be saved
dir_path2 = paste("C:\\Users\\", Sys.getenv("USERNAME"), "\\Dropbox\\computing\\R_program\\fnl_result", sep = "") 
# Source R code containing the functions
setwd(dir_path)
source("functions_to_be_sourced.R")
# Set parameters for generating data:
n = 300
p1 = 70 / 300
p2 = 100 / 300
parameters = list()  # a list for storing parameter values below
# Calculate the number of cores for parallel computing
no_cores = detectCores()/2 - 1
# Number of parallel tasks for computing nrep (see below) results
parameters$nsplit = no_cores
# Number of datasets to be generated; should be a factor of nsplit above
parameters$nrep = 2 * parameters$nsplit
# Number of datasets to be generated per parallel task
parameters$nrep_sub = ceiling(parameters$nrep / parameters$nsplit)
# Number of bootstrap samples
parameters$nboot = 1000  
# First entry is main alpha level of interest; second entry is the secondary alpha level
parameters$alpha_vec = c(0.05, 0.01)
# [0, tau]: the time domain of Xt process
parameters$tau = 1
# Grid solution for time domain of Xt
parameters$gridsol = 0.001
# Specify Ornstein--Uhlenbeck process for Xt generation using the simdiff function of package ESGtoolkit
parameters$SDEmodel = "OU"
# Rescale the Ornstein--Uhlenbeck process so that it has similar range of activity levels as NHANES data
parameters$Xtscale = 300
# Parameters for Ornstein--Uhlenbeck process in the simdiff function of package ESGtoolkit: (CIR1, CIR2, CIR3)
# We specify (CIR3 ^ ) / CIR2, CIR2, CIR1 / CIR2 instead, based on the mean and variance
# formulas of Ornstein--Uhlenbeck process 
parameters$CIR3sq_div_CIR2 =  c(64/24, 64/24 * 1.015, 64/24 * 1.05)
parameters$CIR2 = c(2, 1, 0.5)
parameters$CIR1_div_CIR2 = c(-0.4, -0.4022, -0.4022)
# Parameters for beta random variable generation: (beta1, beta2)
# We specify beta1 / beta2, beta2 instead, which allows flexible control of the variance while preserving the mean
parameters$beta1_div_beta2 = c(10, 10, 10)  # same, then mean same
parameters$beta2 = c(10, 10, 10)
# The activity profiles are represented at a resolution of n_agrid equidistant points
parameters$n_agrid = 100
# Parameters related to initial value generation in computing -2logR(a)
parameters$n_mu_grid = 20  # number of grid points for mu initialization
parameters$n_lamb_grid = 20  # number of grid points for lambdas initialization
parameters$n_grid_incbd = 6  # max number of times the number of grid points can be updated
parameters$mu_tol_grid = 20
parameters$lamb_tol_grid = 20
# [m_1, m_2] / tau
parameters$mu_hat_bounds = c(0.05, 0.95)
# starting seed for each replication of data generation
parameters$seedstart = 0
############## END User input

############# Other parameter calculations based on user input
parameters$CIR1 = parameters$CIR1_div_CIR2 * parameters$CIR2  
parameters$CIR3 = sqrt(parameters$CIR3sq_div_CIR2 * parameters$CIR2)
parameters$beta1 = parameters$beta1_div_beta2 * parameters$beta2
############# END Other parameter calculations based on user input

############# Producing p-value of the EL test, given a specific dataset
### making an artificial dataset Ta, by generating Xt from the positive part of an Ornstein--Uhlenbeck process, 
### and multiplied the resulting activity profile by an independent beta random variable:
n_subject = c(ceiling(n * p1), ceiling(n * p2), n - ceiling(n * p1) - ceiling(n * p2))
tau = parameters$tau
n_boot = parameters$nboot
alpha_vec = parameters$alpha_vec
gridsol = parameters$gridsol
t_ncol = tau / gridsol + 1
Xt = list()
beta_to_time = list()

as_all_ranges = sapply(1:length(n_subject), FUN = function (j) {
  sim.CIR <- simdiff(n = n_subject[j], horizon = (t_ncol - 1), frequency = "annual", model = parameters$SDEmodel, x0 = 0, theta1 = parameters$CIR1[j], theta2 = parameters$CIR2[j], theta3 = parameters$CIR3[j])
  Xt[[j]] <<- parameters$Xtscale * pmax(t(sim.CIR), 0)
  

  as_all = sort(c((Xt[[j]])[, -1]))  # -1 bec X0 same for all
  as_quant = quantile(as_all, probs = c(1-parameters$mu_hat_bounds[1], 1 - parameters$mu_hat_bounds[2]))
  beta_to_time[[j]] <<- rbeta(n_subject[j], parameters$beta1[j], parameters$beta2[j])
  ind_keep = (as_all < as_quant[1] & as_all >= as_quant[2])
  range(as_all[ind_keep])
})

as = seq(max(as_all_ranges[1, ]), min(as_all_ranges[2, ]), length = parameters$n_agrid)

Ta = lapply(1:length(n_subject), FUN = function (j) {
  t(sapply(1:n_subject[j], FUN = function (i) {
    observed_Xt = (Xt[[j]])[i, -1][!is.na((Xt[[j]])[i, -1])]
    observed_Xt_len = length(observed_Xt)
    grid_width_i = rep(1 / observed_Xt_len, observed_Xt_len)  # 1 / observed_Xt_len: adjusting for missing data
    mu_a_hat(observed_Xt, grid_widths = grid_width_i, as) * (beta_to_time[[j]])[i]
  }))
})

### implement the function computing p-values and related results
test_out = procedure(Ta = Ta, as = as, n_boot = n_boot, alpha_vec = alpha_vec)
### outputs from the function:
test_out$suptest  # the EL test statistic
test_out$sup_EL_crit  # critical value of the EL test
test_out$out_sup_pval  # p-value of the EL test
############# END Producing p-value of the EL test,  given a specific dataset

############# Power function
powerfn = function(n, parameters, p1, p2, split) {
# Computes the empirical rejection rate of our composite procedure
# Args:
#   n: total sample size 
#   parameters: the parameters input/calculated in the codes above (up to # END Other parameter calculations based on user input)
#   p1: proportion of data for the first group
#   p2: proportion of data for the second group
#   split: number of parallel tasks for computing nrep (see parameters$nrep above) results
# Returns:
#   the empirical rejection rate of the EL test
 ### Parameter organization
  n_subject = c(ceiling(n * p1), ceiling(n * p2), n - ceiling(n * p1) - ceiling(n * p2))
  nrep_sub = parameters$nrep_sub
  tau = parameters$tau
  n_boot = parameters$nboot
  alpha_vec = parameters$alpha_vec
  gridsol = parameters$gridsol
  t_ncol = tau / gridsol + 1
 ### Define output
  suptest_pval_all = 1:nrep_sub * 0
 ### For loop generating nrep_sub number of datasets and computing statistics accordingly
  for (i in 1:nrep_sub) {
    set.seed(parameters$seedstart + nrep_sub * (split - 1) + i) 
    print(parameters$seedstart + nrep_sub * (split - 1) + i)
     ## data generation
    Xt = list()
    beta_to_time = list()
    as_all_ranges = sapply(1:length(n_subject), FUN = function (j) {
    sim.CIR <- simdiff(n = n_subject[j], horizon = (t_ncol - 1), frequency = "annual", model = parameters$SDEmodel, x0 = 0, theta1 = parameters$CIR1[j], theta2 = parameters$CIR2[j], theta3 = parameters$CIR3[j])
    Xt[[j]] <<- parameters$Xtscale * pmax(t(sim.CIR), 0)
    as_all = sort(c((Xt[[j]])[, -1]))  # -1 bec X0 same for all
    as_quant = quantile(as_all, probs = c(1-parameters$mu_hat_bounds[1], 1 - parameters$mu_hat_bounds[2]))
    beta_to_time[[j]] <<- rbeta(n_subject[j], parameters$beta1[j], parameters$beta2[j])
    ind_keep = (as_all < as_quant[1] & as_all >= as_quant[2])
    range(as_all[ind_keep])
    })
    as = seq(max(as_all_ranges[1, ]), min(as_all_ranges[2, ]), length = parameters$n_agrid)
    Ta = lapply(1:length(n_subject), FUN = function (j) {
      t(sapply(1:n_subject[j], FUN = function (i) {
      observed_Xt = (Xt[[j]])[i, -1][!is.na((Xt[[j]])[i, -1])]
      observed_Xt_len = length(observed_Xt)
      grid_width_i = rep(1 / observed_Xt_len, observed_Xt_len)  # 1 / observed_Xt_len: adjusting for missing data
      mu_a_hat(observed_Xt, grid_widths = grid_width_i, as) * (beta_to_time[[j]])[i]
      }))
    })     
 ## Computing results of the EL test
    print(Sys.time())
    test_out = procedure(Ta = Ta, as = as, n_boot = n_boot, alpha_vec = alpha_vec)
    print(Sys.time())
 ## P-value of the proposed EL test
    suptest_pval_all[i] = test_out$out_sup_pval
  }  # END for
# Finally get rej rate
  rp = sum(suptest_pval_all < 0.05) / nrep_sub
  return(rp) 
}  # END powerfn
############# END Power function

############# Calculating power given n,  p1,  p2 (utilizing R parallel computing)
# Initiate cluster
cl = makeCluster(no_cores)
registerDoParallel(cl)
### Compute the results:
time1 = Sys.time()
foreach(split = 1:parameters$nsplit,  
  .combine = c, 
  .packages = c("parallel", "foreach", "iterators", "doParallel", "rootSolve", "Rcpp", "plyr", "coda", "numDeriv", "bbmle", "ESGtoolkit")
) %dopar% {
  out = powerfn(n = n, parameters, p1 = p1, p2 = p2, split)
  setwd(dir_path2)
  save(out,  file = paste("n_", n, "_split_", split, '_nrep_', parameters$nrep, "_p1_", p1, "_p2_", p2, ".Rdata", sep = ""))
}  # END foreach
time2 = Sys.time()
time2-time1 
stopCluster(cl)
### Read the results:
rp_vec = 1:parameters$nsplit * 0
setwd(dir_path2)
for (split in 1:parameters$nsplit) {
  load(paste("n_", n, "_split_", split, '_nrep_', parameters$nrep, "_p1_", p1, "_p2_", p2, ".Rdata", sep=""))
  rp_vec[split] = out
}  # END for
mean(rp_vec)  # the desired power
#############END Calculating power given n,  p1,  p2

