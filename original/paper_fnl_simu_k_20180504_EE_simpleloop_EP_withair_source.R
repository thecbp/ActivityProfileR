######################## Required Packages ########################
#Rlib_loc=paste("C:\\Users\\",Sys.getenv("USERNAME"),"\\Dropbox\\computing\\R_package\\Rlib",sep="")
#for cluster computing:
#install.packages('parallel')
library(parallel)  # inside R core since 2.14.0
#install.packages('foreach')
library(foreach)
#library(iterators)
#install.packages('doParallel')
library(doParallel)
#for root solving
#install.packages('rootSolve')
library(rootSolve)
#library(Rcpp)
#install.packages('plyr')
#library(plyr) # 20181201 try removing it as it causes error
#library(coda)
#library(numDeriv)
#library(bbmle)
#install.packages('emdbook')
#library(emdbook) # 20181201 try removing it as it causes error
#install.packages('fda.usc')
library(fda.usc)
######################## User Inputs ########################
parameters = list()
# number of bootstrap samples
parameters$nboot = 1000  
# first entry is main alpha level of interest; second entry is the secondary alpha level
parameters$alpha_vec = c(0.05, 0.01)
# [0, tau]: the time domain of Xt process
parameters$tau = 1
# n_agrid = 100 gives 1~99 percentiles, = 10 gives 10~90 quantiles etc
#parameters$n_agrid = 100
# number of grid points for mu initialization
parameters$n_mu_grid = 20
# number of grid points for lambdas initialization
parameters$n_lamb_grid = 20
# max number of times the number of grid points can be updated
parameters$n_grid_incbd = 6
# grid to be away from max of mininums in each group, and min of maximums in each group
parameters$mu_tol_grid = 20
# grid to be away from lower and upper bounds of lambdas
parameters$lamb_tol_grid = 20
# lower and upper bounds for mu_hat_all to form [alpha_1, alpha_2] 
parameters$mu_hat_bounds = c(0.05, 0.95)
######################## Functions ########################
# to be applied to vector of accelerometer readings, pick > 0 columns, take log, times c, and sum it up...
g0tcsum = function(vec, effc) {
  sum(log(vec[vec > 0]) * effc[vec > 0])
}  # END g0tcsum
# activity profile generation, given functional data, for discrete t***
mu_a_hat = function(Xt, grid_widths, as) {
# Computes the activity profile of a stochastic process of accelerometer readings
# Args:
#   Xt: a vector of accelerometer readings on a grid of time points
#   grid_widths: the number of time units each grid point contains
#   Xt, grid_widths should be vectors of the same length
#   as: the grid of a's
# Returns:
#   activity profile of Xt
# check how many y - values there is == the number of a values
  n_a = length(as)
  act_profile = 1:n_a * 0
  for (a in 1:n_a) {
    act_profile[a] = sum(grid_widths[Xt > as[a]])
  }  # END for
  return(act_profile)
}  # END mu_a
# testing with as defined by Xt itself:
# mu_a(Xt, rep(1, length(Xt)), sort(unique(Xt))[-length(unique(Xt))])

#, mu_hat_vec
neg2logR_test = function(ai, Ta, n_subject, EL_testcrit = 0, error_out = 0) {  
# Computes -2logR(a)
# Args:
#   ai: index for the given activity level a
#   Ta: a list with Ta[[j]] being the matrix of activity profile data, i-th row representing the i-th subject's activity profile
#   (dim = n_subject[j] x length(as[[j]]))
#   n_subject: a vector, n_subject[j] being the sample size for the j-th group
#   mu_hat_vec: a length(n_subject) x length(as) matrix of estimated mu(a) with mu_hat_vec[j, ] being estimated
#   mean activity levels across a for the j-th group
#   EL_testcrit: critical value for testing H0
# Returns:
#   -2logR(a)
  n = sum(n_subject)
  gamma_js = n_subject / sum(n_subject)
  n_group = length(n_subject)
  min_ups_max_lows = matrix(sapply(1:n_group, FUN = function (j) {
    range((Ta[[j]])[, ai])
  }), nrow = 2)  # 1st row: min of each group; 2nd row: max of each group
##############dealing with non-overlapping regions (include the case when the region is degenerate, i.e., the = case)
  if (max(min_ups_max_lows[1, ]) >= min(min_ups_max_lows[2, ])) {
    if (error_out == 0) {
      return(Inf)
    } else {
      return(list(neg2logR_crit = Inf, 
        still_error = 0   # no status and roots output
      ))
    }
  } 
##############END dealing with non-overlapping regions
# 0726 air: to-do: test while code below
  # related to while loop update  # 0726 air
  n_grid_inc = 0  # 0726 air
  still_error = 1
  test1t = Inf  # 0726 air: if no solution below, then means sup{empty set} = 0 in numerator of -2logR
  while (still_error == 1 & n_grid_inc <= parameters$n_grid_incbd) {  # 0726 air: max of 7 whiles available
    n_grid_inc = n_grid_inc + 1  # 0726 air
    parameters$n_mu_grid = parameters$n_mu_grid * n_grid_inc + 1  # 0726 air
    parameters$n_lamb_grid = parameters$n_lamb_grid * n_grid_inc + 1  # 0726 air
  # END related to while loop udpate  # 0726 air
##############grid construction
    mu_tol = (min(min_ups_max_lows[2, ]) - max(min_ups_max_lows[1, ])) / parameters$mu_tol_grid
    init_mu_a_grid = seq(max(min_ups_max_lows[1, ]) + mu_tol, min(min_ups_max_lows[2, ]) - mu_tol, length = parameters$n_mu_grid)
    ll_l = matrix(0, nrow = parameters$n_mu_grid, ncol = n_group)
    ul_l = matrix(0, nrow = parameters$n_mu_grid, ncol = n_group)
    T_js = list()
    init_Delta_grid = array(unlist(sapply(1:parameters$n_mu_grid, FUN = function (init_mu_ind) {
      init_mu_a = init_mu_a_grid[init_mu_ind]
      sapply(1:n_group, FUN = function (j) {  
      # only (n_group - 1) initialization of Delta needed; 
      # but want ll_l[n_group] and ul_l[n_group] to be calculated
        T_js[[j]] <<- (Ta[[j]])[, ai]
        g_j = ((Ta[[j]])[, ai] - init_mu_a) / gamma_js[j]  # a vector of length n_j
        ll_l[init_mu_ind, j] <<- max((1/n_subject[j] - 1) / g_j[g_j > 0])
        ul_l[init_mu_ind, j] <<- min((1/n_subject[j] - 1) / g_j[g_j < 0])
        lamb_tol = (ul_l[init_mu_ind, j] - ll_l[init_mu_ind, j]) / parameters$lamb_tol_grid
        seq(ll_l[init_mu_ind, j] + lamb_tol, ul_l[init_mu_ind, j] - lamb_tol, length = parameters$n_lamb_grid) 
      })
    })), c(parameters$n_lamb_grid, n_group, parameters$n_mu_grid))
###############END grid construction

##############initialization
# 1. in R, better to just do the matrix (with proper replications)
# then sweep through the rows
# instead of multiple for loops

# all possible initial values in mat1; 
# then use init_check for sum(Delta_j) = 0 check); 
# first column is init_mu and the rest is init Delta
    init_check = matrix(0, nrow = parameters$n_mu_grid, ncol = parameters$n_lamb_grid * (parameters$n_lamb_grid ^ (n_group - 2)))
    init_Delta_grid_expand = matrix(unlist(lapply(1:parameters$n_mu_grid, FUN = function (init_mu_ind) {
      init_mu_a = init_mu_a_grid[init_mu_ind]
      mat1 = cbind(init_mu_a, sapply(1:(n_group - 1), FUN = function (j) {
        rep(rep(init_Delta_grid[, j, init_mu_ind],each = parameters$n_lamb_grid ^ (j - 1)), length = parameters$n_lamb_grid * (parameters$n_lamb_grid ^ (n_group - 2)))
      }))
# ** as scratch paper 0619 p.1, check if \lambda_{n_group - 1} 
# = -sum(init_Delta_grid_expand[, -1, init_mu_ind])
# btw ll_l[init_mu_ind, n_group] and ul_l[init_mu_ind, n_group]
      init_lamb_km1 = -apply(as.matrix(mat1[, -1]), 1, sum)
      init_check[init_mu_ind, ] <<- (init_lamb_km1 > ll_l[init_mu_ind, n_group] & init_lamb_km1 < ul_l[init_mu_ind, n_group])
      return(t(mat1[(init_check[init_mu_ind, ] ==1), ]))
    })), byrow = T, ncol = n_group)
#, c(parameters$n_lamb_grid * (parameters$n_lamb_grid ^ (n_group - 2)), n_group, parameters$n_mu_grid))
# 2nd dimension: (n_group - 1) + 1st mu dimension

    n_init = dim(init_Delta_grid_expand)[1]
    status_grid = matrix(NA, nrow = n_init, ncol = 2 + 3 * n_group - 1)  # use rbind later?
    roots = matrix(NA, nrow = n_init, ncol = 2 * n_group - 1) #details
    multir = list()  # 0726 air
    multir$estim.precis = 1  # 0726 air
    min_ps = rep(-1, times = n_group)  # 0726 air (for all NA output)
#  time1 = Sys.time()
    for (init_ind in 1:n_init) {
      init_mu_a = init_Delta_grid_expand[init_ind, 1]
      init_lambda_js = cumsum(-init_Delta_grid_expand[init_ind, -1])
      init_gam = gamma_js - diff(c(0, init_lambda_js, 0)) * (-init_mu_a)
      lamb_gam = c(init_lambda_js, init_gam)
      options(warn = 2)	#warning=error! will stop the whole function; so use tryCatch-->ok to stop in that
      possibleError <- tryCatch(
        multir<-multiroot(f = model, start = lamb_gam, lbs = T_js, n_subject = n_subject, maxeval = 10000), #can input multir here!
        error=function(e) e
      )
# !inherits(possibleError, "error") # T if there's no error
      options(warn = 0)
      if (!inherits(possibleError, "error")) {  # then multir has meaningful values
        roots[init_ind, ] = multir$root
        lamb = roots[init_ind, ][1:(n_group - 1)]
        gam = roots[init_ind, ][-(1:(n_group - 1))]
        lamb_kp1 = c(0, lamb, 0)  # j-th element is (j-1)-th element of lambdas
        p = list()
        min_ps = sapply(1:n_group, FUN = function (j) {
          p[[j]] <<- 1 / (n * (gam[j] + (lamb_kp1[j] - lamb_kp1[j + 1]) * T_js[[j]]))
          min(p[[j]])
        })
        status_grid[init_ind, ] = c(!inherits(possibleError, "error"), multir$estim.precis, multir$f.root, min_ps)
        if (multir$estim.precis < 10 ^ (-8) & !inherits(possibleError, "error") & sum(min_ps < 0) == 0) break  # decision_stop
      }  # END if
#    print(init_ind)
#    print(Sys.time())
    }  # END for
#  time2 = Sys.time()
#  print(time2 - time1)
    still_error = !(multir$estim.precis < 10 ^ (-8) & !inherits(possibleError, "error") & sum(min_ps < 0) == 0)
  }  # END while  # 0726 air
  if (exists("p") == 1) {  # 0726 air
    test1t = -2 * sum(unlist(lapply(1:n_group, FUN = function (j) {
      return(log(p[[j]]))
    }))) - 2 * sum(n_subject * log(n_subject))
  }
  if (error_out == 0) {
    return(test1t - EL_testcrit)
  } else {
    return(list(neg2logR_crit = test1t - EL_testcrit, 
      still_error = still_error, 
      status = status_grid[init_ind, ], 
      roots = roots[init_ind, ])
    )
  }
# can't use the following for uniroot below:
# return(list(out = test1t - EL_CBcrit, 
#   still_error = still_error,
#   num_lambda = num_neg_loglik$lambda))
}  # END neg2logR_CB


#######################
procedure = function(Ta, as, n_boot = 1000, alpha_vec) { 
# Computes confidence bands for EL, EP, and HW procedures and whether they capture the truth mu_a for all as
# Args:
#   Ta: a list with Ta[[j]] being the matrix of activity profile data, i-th row representing the i-th subject's activity profile
#   (dim = n_subject[j] x length(as[[j]]))
#   as: a vector of activity levels a occurred in the Ta dataset (from smallest to largest)
#   n_boot: number of bootstrap samples
#   mu_a: a length(n_subject) x length(as) matrix of truth (if known) with mu_a[j, ] being 
#   mean activity levels across a for the j-th group
#   alpha_vec: first entry is main alpha level of interest; second entry is the secondary alpha level
# Returns: (to-do)
#   EL confidence band construction; 0 otherwise
#   EL_CB: EL confidence band, with first row being the lower bound and second row being the upper bound
#   EL_CBcrit: critical value for computing EL confidence band
#   EL_CBrej: 1 if the truth mu_a is not included in the EL confidence band; 0 otherwise
#   EP_CB: EP confidence band, with first row being the lower bound and second row being the upper bound
#   EP_CBcrit: critical value for computing EP confidence band
#   EP_CBrej: 1 if the truth mu_a is not included in the EP confidence band; 0 otherwise
#   HW_CB: HW confidence band, with first row being the lower bound and second row being the upper bound
#   HW_CBcrit: critical value for computing HW confidence band
#   HW_CBrej: 1 if the truth mu_a is not included in the HW confidence band; 0 otherwise
#   mu_hat_vec: sample mean of mu_a
#   EL_CB_aicheck: 1 if the EL band reaches the maximum of Tai - a tolerance level 10 ^ (-tolfac_ai)
#   at a; 0 otherwise
#   tolfac_us: a vector of the maximal tolerance 10 ^ (-tolfac_u) (initialize tolfac_u as 3) needed for 
#   -2logR(maximum of Tai - 10 ^ (-tolfac_u), a) to have different sign 
#   from -2logR(sample mean of mu_a, a), for each a
#   tolfac_ls: a vector of the maximal tolerance 10 ^ (-tolfac_l) (initialize tolfac_l as 3) needed for 
#   -2logR(minimum of Tai + 10 ^ (-tolfac_l), a) to have different sign 
#   from -2logR(sample mean of mu_a, a), for each a
#   neg2logR_CB_decisionstop_error: 1 if there is error in computing -2logR(\tilde{\mu}(a), a) during 
  n_Ta = length(as)  # the length of as
  n_subject = sapply(Ta, FUN = function (x) {
    dim(x)[1]
  })  # a vector, n_subject[j] being the sample size for the j-th group
  gamma_js = n_subject / sum(n_subject)

##******S2_hat can't be 0 in EL! So theta_hat_js can be in the demoninator! (new agrid specification
## might solve the problem

  U_hat_star = array(0, c(n_boot, n_Ta, length(n_subject)))
  theta_hat_js = array(0, c(n_boot, n_Ta, length(n_subject)))
  wjs = array(0, c(n_boot, n_Ta, length(n_subject)))
  mu_hat_vec = matrix(0, nrow = length(n_subject), ncol = n_Ta)
  S2_hat_vec = matrix(0, nrow = length(n_subject), ncol = n_Ta)
  boot_indx = list()
  U_hat_star_run = sapply(1:length(n_subject), FUN = function (j) {
    boot_indx[[j]] <<- sample(1:n_subject[j], n_subject[j] * n_boot, replace = T)  # suppose n_subject[j] subjects selected, then repeated n_boot times
    boot_indx_mat = matrix(boot_indx[[j]], byrow = T, nrow = n_boot, ncol = n_subject[j])
    mu_hat = matrix(rep(apply(Ta[[j]], 2, mean), each = n_boot), nrow = n_boot, ncol = n_Ta)
    mu_hat_vec[j, ] <<- mu_hat[1, ]
    Ta_boot = array((Ta[[j]])[boot_indx[[j]], ], c(n_subject[j], n_boot, n_Ta))
    mu_hat_star = apply(Ta_boot, c(2, 3) , mean)
    S2_hat = matrix(apply(Ta[[j]], 2, var) * (n_subject[j] - 1) / n_subject[j], byrow = T, nrow = n_boot, ncol = n_Ta)
    S2_hat_vec[j, ] <<- S2_hat[1, ]
  # a can't be max(Ta) so that 
  # Ta[, ai] == 0 vec and mu_hat == 0 so that S2_hat == 0
  # later (not proved yet): can try apply(Ta[boot_indx[[j]], ], 2, var) * ...
# 0418  S2_hat_star = apply(Ta_boot, c(2, 3), var) * (n_subject - 1) / n_subject
    U_hat_star[, , j] <<- sqrt(n_subject[j]) * (mu_hat_star -  mu_hat) / sqrt(S2_hat)  # n_boot x n_Ta
# 0418   U_hat_star = sqrt(n_subject) * (mu_hat_star -  mu_hat) / sqrt(S2_hat_star)  # n_boot x n_Ta
    theta_hat_js[, , j] <<- S2_hat / gamma_js[j]
    wjs[, , j] <<- 1 / theta_hat_js[, , j]  # need to be normalized to 1 later; but just once for all bootstrap samples
    return(0)  # like Ujps in ELSOc_k
  })
# normalizing wjs to 1: since it's the same for all bootstrap samples, just do it once
  wjs = wjs / array(rep(apply(wjs[1, , ], 1, sum), each = n_boot), c(n_boot, n_Ta, length(n_subject)))
  wjs_vec = t(wjs[1, , ])
# check: apply(wjs[1, , ], 1, sum) == a vector of 1's, ok
  avg_U_hat_star = apply(sqrt(wjs) * U_hat_star, c(1, 2), sum)
  avg_U_hat_star_karray = array(rep(avg_U_hat_star, times = length(n_subject)), c(n_boot, n_Ta, length(n_subject)))
  U2_boot = apply(wjs * (U_hat_star / sqrt(wjs) - avg_U_hat_star_karray) ^ 2, c(1, 2), sum)  # bootstrap SSB(a) in Theorem 1, n_boot x n_Ta
### bootstrap sup_a SSB(a) and int_a SSB(a) in Theorem 1:
  sup_boot = apply(as.matrix(U2_boot), 1, max)  # bootstrap sup_{a \in [alpha_1, alpha_2]} SSB(a) in Theorem 1
  sup_EL_crit = as.vector(quantile(sup_boot, probs = 1 - alpha_vec)) 
# finished small eg checking
  a_big = matrix(rep(as, times = n_boot), byrow = TRUE, nrow = n_boot)  # n_boot x n_Ta
  U2_boot_times_da = as.matrix(U2_boot[, -n_Ta]) * as.matrix(t(apply(a_big, 1, diff)))
  int_da_boot = apply(as.matrix(U2_boot_times_da), 1, sum)
  int_da_EL_crit = as.vector(quantile(int_da_boot, probs = 1 - alpha_vec))
### EL test:
  teststat_pre = 1:n_Ta * 0 
  error_vec = 1:n_Ta * 0
  for (ai in 1:n_Ta) {
    neg2logRa = neg2logR_test(ai, Ta, n_subject, EL_testcrit = 0, error_out = 1) 
    error_vec[ai] = neg2logRa$still_error
    if (error_vec[ai] == 1) next 
    teststat_pre[ai]=neg2logRa$neg2logR_crit
    print(ai)
    print(Sys.time())
  }  # END for
  suptest = max(teststat_pre)
  inttest_pre_da = teststat_pre[-n_Ta] * diff(as)
  inttest_da = sum(inttest_pre_da) 
# EP:
  mu_bar = apply(matrix(gamma_js, nrow = length(n_subject), ncol = n_Ta) * mu_hat_vec, 2, sum)
  mu_bar_vec = matrix(mu_bar, byrow = T, nrow = length(n_subject), ncol = n_Ta)
  n_subject_vec = matrix(n_subject, nrow = length(n_subject), ncol = n_Ta)
  U_hat = sqrt(n_subject_vec) * (mu_hat_vec -  mu_bar_vec) / sqrt(S2_hat_vec)
  avg_U_hat = apply(sqrt(wjs_vec) * U_hat, 2, sum)
  avg_U_hat_vec = matrix(avg_U_hat, byrow = T, nrow = length(n_subject), ncol = n_Ta)
  U2_vec = apply(wjs_vec * (U_hat / sqrt(wjs_vec) - avg_U_hat_vec) ^ 2, 2, sum)  # estimate of SSB(a) in Theorem 1, n_boot x n_Ta
  suptest_EP = max(U2_vec)
  inttest_EP_pre_da = U2_vec[-n_Ta] * diff(as)
  inttest_EP_da = sum(inttest_EP_pre_da) 
# 20180906 anova_onefactor_2004
  Ta_mat = matrix(unlist(lapply(1:length(n_subject), FUN = function (j) {  # for fda.usc's fdata object construction
    c(t(Ta[[j]]))
  })), byrow = T, ncol = parameters$n_agrid)
  anova_onefactor_group = as.factor(rep(1:length(parameters$n_subject), times = parameters$n_subject))
  anova_onefactor_fdata = fdata(Ta_mat, argvals = as, names = list(main = "activity profiles", xlab = "activity level (a)", ylab = "Ta"))  # this object is a list
  anova_onefactor = anova.onefactor(anova_onefactor_fdata, anova_onefactor_group, nboot = parameters$nboot, plot = FALSE)
  return(list(
    out_anova_onefactor_pval = anova_onefactor$pvalue,
    suptest = suptest,
    inttest_da = inttest_da,
    suptest_EP = suptest_EP,
    inttest_EP_da = inttest_EP_da,
    sup_EL_crit = sup_EL_crit,
    int_da_EL_crit = int_da_EL_crit,
    out_sup_pval = mean(sup_boot >= suptest),  
    out_dF_pval = mean(int_da_boot >= inttest_da),
    out_sup_EP_pval = mean(sup_boot >= suptest_EP),  
    out_dF_EP_pval = mean(int_da_boot >= inttest_EP_da),
    as = as,
    mu_hat_vec = mu_hat_vec,
    sup_boot = sup_boot,
    int_da_boot = int_da_boot
  ))
}  # END procedure

