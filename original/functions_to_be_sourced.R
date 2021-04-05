# Supplement to "Nonparametric Comparisons of Activity Profiles From Wearable Device Data" by Chang and McKeague
# The following code provides functions needed to run the file power_calculations.R
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
  act_profile = sapply(1:n_a, FUN = function (a) {
    sum(grid_widths[Xt > as[a]])
  })  # END sapply
  return(act_profile)
}  # END mu_a


model = function(lamb_gam, lbs, n_subject) {
  # Computes estimating equations in -2logR(a) calculation
  # Args:
  #   lamb_gam: a vector of 2 * n_group - 1 length, containing c(lambdas, gammas) mentioned in Appendix B of the paper
  #   lbs: Ta[[j]][, ai] for j = 1, ..., n_group
  #   n_subject: a vector of n_group length, with j-th entry being sample size for the j-th group
  # Returns:
  #   out: a vector of of 2 * n_group - 1 length, each entry being each value of the 2 * n_group - 1 
  #     estimating equations in -2logR(a) calculation: 
  #     c(mu_2 - mu_1, ..., mu_k - mu_{k-1}, sum(p_1) - 1, ..., sum(p_k) - 1)
  n = sum(n_subject)
  gamma_js = n_subject / n
  n_group = length(n_subject)
  out = rep(0, times = 2 * n_group - 1)
  lamb = lamb_gam[1:(n_group - 1)]
  gam = lamb_gam[-(1:(n_group - 1))]
  lamb_kp1 = c(0, lamb, 0)  # j-th element is (j-1)-th element of lambdas
  p = list()
  p[[1]] = 1/(n * (gam[1] + (lamb_kp1[1] - lamb_kp1[2]) * lbs[[1]]))
  EEs = sapply(1:n_group, FUN = function (j) {
    out[n_group + j - 1] <<- sum(p[[j]]) - 1
    if (j != n_group){
      p[[j + 1]] <<- 1/(n * (gam[j + 1] + (lamb_kp1[j + 1] - lamb_kp1[j + 2]) * lbs[[j + 1]]))
      out[j] <<- sum(p[[j + 1]] * lbs[[j + 1]]) - sum(p[[j]] * lbs[[j]])
    }
  })
  
  return(out)
}


neg2logR_test = function(ai, Ta, n_subject, error_out = 0) {  
  # Computes -2logR(a)
  # Args:
  #   ai: index for the given activity level a
  #   Ta: a list with Ta[[j]] being the matrix of activity profile data, i-th row representing the i-th subject's activity profile
  #   n_subject: a vector, n_subject[j] being the sample size for the j-th group
  #   error_out: whether more detailed information is to be returned, including whether there is an error,
  # Returns:
  #   -2logR(a) when error_out = 0; when error_out = 1, an addition output
  #   still_error giving 1 if there is an error in computing -2logR(a), 0 otherwise
  n = sum(n_subject)
  gamma_js = n_subject / sum(n_subject)
  n_group = length(n_subject)
  min_ups_max_lows = matrix(sapply(1:n_group, FUN = function (j) {
    range((Ta[[j]])[, ai])
  }), nrow = 2)  # 1st row: min of each group; 2nd row: max of each group
  ### Dealing with non-overlapping regions
  if (max(min_ups_max_lows[1, ]) >= min(min_ups_max_lows[2, ])) {
    if (error_out == 0) {
      return(Inf)
    } else {
      return(list(neg2logR_crit = Inf, 
                  still_error = 0   
      ))
    }
  } 
  ### Get T_{aj} from the Ta list for given a and j, and scale it so that max range differs by 1
  T_js = lapply(1:n_group, FUN = function (j) {   
    (Ta[[j]])[, ai] / (min(min_ups_max_lows[2, ]) - max(min_ups_max_lows[1, ]))  
  })
  min_ups_max_lows = matrix(sapply(1:n_group, FUN = function (j) {
    range(T_js[[j]])
  }), nrow = 2)  # 1st row: min of each group; 2nd row: max of each group
  
  n_grid_inc = 0  # the number of times we have updated the number of grid points for mu and lambda initialization
  still_error = 1  # 1 if there is an error in computing -2logR(a), 0 otherwise
  test1t = Inf  # if no solution below, then means sup{empty set} = 0 in numerator of -2logR
  n_mu_grid = parameters$n_mu_grid  # number of grid points for mu initialization
  n_lamb_grid = parameters$n_lamb_grid  # number of grid points for lambda initialization
  ### A while loop for computing -2logR(a), allowing n_grid_inc to be updated if still_error = 1
  while (still_error == 1 & n_grid_inc <= parameters$n_grid_incbd) {  
    n_grid_inc = n_grid_inc + 1  
    n_mu_grid = n_mu_grid * n_grid_inc + 1  
    n_lamb_grid = n_lamb_grid * n_grid_inc + 1  
    ### Grid construction for mu and lambda initialization
    mu_tol = (min(min_ups_max_lows[2, ]) - max(min_ups_max_lows[1, ])) / parameters$mu_tol_grid
    init_mu_a_grid = seq(max(min_ups_max_lows[1, ]) + mu_tol, min(min_ups_max_lows[2, ]) - mu_tol, length = n_mu_grid)
    ll_l = matrix(0, nrow = n_mu_grid, ncol = n_group)
    ul_l = matrix(0, nrow = n_mu_grid, ncol = n_group)
    init_Delta_grid = array(unlist(sapply(1:n_mu_grid, FUN = function (init_mu_ind) {
      init_mu_a = init_mu_a_grid[init_mu_ind]
      sapply(1:n_group, FUN = function (j) {  
        g_j = (T_js[[j]] - init_mu_a) / gamma_js[j]  # a vector of length n_j
        ll_l[init_mu_ind, j] <<- max((1/n_subject[j] - 1) / g_j[g_j > 0])
        ul_l[init_mu_ind, j] <<- min((1/n_subject[j] - 1) / g_j[g_j < 0])
        lamb_tol = (ul_l[init_mu_ind, j] - ll_l[init_mu_ind, j]) / parameters$lamb_tol_grid
        seq(ll_l[init_mu_ind, j] + lamb_tol, ul_l[init_mu_ind, j] - lamb_tol, length = n_lamb_grid) 
      })
    })), c(n_lamb_grid, n_group, n_mu_grid))
    ### Mu and lambda initialization 
    init_check = matrix(0, nrow = n_mu_grid, ncol = n_lamb_grid * (n_lamb_grid ^ (n_group - 2)))
    init_Delta_grid_expand = matrix(unlist(lapply(1:n_mu_grid, FUN = function (init_mu_ind) {
      init_mu_a = init_mu_a_grid[init_mu_ind]
      mat1 = cbind(init_mu_a, sapply(1:(n_group - 1), FUN = function (j) {
        rep(rep(init_Delta_grid[, j, init_mu_ind],each = n_lamb_grid ^ (j - 1)), length = n_lamb_grid * (n_lamb_grid ^ (n_group - 2)))
      }))
      init_lamb_km1 = -apply(as.matrix(mat1[, -1]), 1, sum)
      init_check[init_mu_ind, ] <<- (init_lamb_km1 > ll_l[init_mu_ind, n_group] & init_lamb_km1 < ul_l[init_mu_ind, n_group])
      return(t(mat1[(init_check[init_mu_ind, ] ==1), ]))
    })), byrow = T, ncol = n_group)
    ### Find gamma and lambda satisfying the estimating equations related to optimization in -2logR(a) 
    n_init = dim(init_Delta_grid_expand)[1]
    status_grid = matrix(NA, nrow = n_init, ncol = 2 + 3 * n_group - 1)  
    roots = matrix(NA, nrow = n_init, ncol = 2 * n_group - 1) 
    multir = list()  
    multir$estim.precis = 1  
    min_ps = rep(-1, times = n_group)  
    for (init_ind in 1:n_init) {
      init_mu_a = init_Delta_grid_expand[init_ind, 1]
      init_lambda_js = cumsum(-init_Delta_grid_expand[init_ind, -1])
      init_gam = gamma_js - diff(c(0, init_lambda_js, 0)) * (-init_mu_a)
      lamb_gam = c(init_lambda_js, init_gam)
      options(warn = 2)	 # warning = error! will stop the whole function
      possibleError <- tryCatch(
        multir<-multiroot(f = model, start = lamb_gam, lbs = T_js, n_subject = n_subject), 
        error=function(e) e
      )
      options(warn = 0)
      if (!inherits(possibleError, "error")) {  # !inherits(possibleError, "error") == T if there's no error
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
        if (multir$estim.precis < 10 ^ (-8) & !inherits(possibleError, "error") & sum(min_ps < 0) == 0) break  
      }  # END if
    }  # END for
    still_error = !(multir$estim.precis < 10 ^ (-8) & !inherits(possibleError, "error") & sum(min_ps < 0) == 0)
  }  # END while  
  ### Compute -2logR(a) given the p output above
  if (exists("p") == 1) {  
    test1t = -2 * sum(unlist(lapply(1:n_group, FUN = function (j) {
      return(log(p[[j]]))
    }))) - 2 * sum(n_subject * log(n_subject))
  }  # END if
  if (error_out == 0) {
    return(test1t)
  } else {
    return(list(neg2logR_crit = test1t, 
                still_error = still_error
    ))
  }  # END if/else
}  # END neg2logR_CB


procedure = function(Ta, as, n_boot = 1000, alpha_vec) { 
  # Performs the EL test
  # Args:
  #   Ta: a list with Ta[[j]] being the matrix of activity profile data, i-th row representing the i-th subject's activity profile
  #   as: a vector of activity levels a occurred in the Ta dataset (from smallest to largest)
  #   n_boot: number of bootstrap samples
  #   alpha_vec: first entry is main alpha level of interest; second entry is the secondary alpha level
  # Returns:
  #    suptest: the EL test statistic
  #    sup_EL_crit: critical value of the EL test
  #    out_sup_pval: p-value of the EL test   
  n_Ta = length(as)  # the length of as
  n_subject = sapply(Ta, FUN = function (x) {
    dim(x)[1]
  })  # a vector, n_subject[j] being the sample size for the j-th group
  gamma_js = n_subject / sum(n_subject)
  U_hat_star = array(0, c(n_boot, n_Ta, length(n_subject)))
  theta_hat_js = array(0, c(n_boot, n_Ta, length(n_subject)))
  wjs = array(0, c(n_boot, n_Ta, length(n_subject)))
  mu_hat_vec = matrix(0, nrow = length(n_subject), ncol = n_Ta)
  S2_hat_vec = matrix(0, nrow = length(n_subject), ncol = n_Ta)
  boot_indx = list()
  
  U_hat_star_run = sapply(1:length(n_subject), FUN = function (j) {
    boot_indx[[j]] <<- sample(1:n_subject[j], n_subject[j] * n_boot, replace = T)  
    boot_indx_mat = matrix(boot_indx[[j]], byrow = T, nrow = n_boot, ncol = n_subject[j])
    mu_hat = matrix(rep(apply(Ta[[j]], 2, mean), each = n_boot), nrow = n_boot, ncol = n_Ta)
    mu_hat_vec[j, ] <<- mu_hat[1, ]
    Ta_boot = array((Ta[[j]])[boot_indx[[j]], ], c(n_subject[j], n_boot, n_Ta))
    mu_hat_star = apply(Ta_boot, c(2, 3) , mean)
    S2_hat = matrix(apply(Ta[[j]], 2, var) * (n_subject[j] - 1) / n_subject[j], byrow = T, nrow = n_boot, ncol = n_Ta)
    S2_hat_vec[j, ] <<- S2_hat[1, ]
    U_hat_star[, , j] <<- sqrt(n_subject[j]) * (mu_hat_star -  mu_hat) / sqrt(S2_hat)  # n_boot x n_Ta
    theta_hat_js[, , j] <<- S2_hat / gamma_js[j]
    wjs[, , j] <<- 1 / theta_hat_js[, , j]  # need to be normalized to 1 later; but just once for all bootstrap samples
    return(0)  
  })
  # Normalizing wjs to 1: since it's the same for all bootstrap samples, just do it once
  wjs = wjs / array(rep(apply(wjs[1, , ], 1, sum), each = n_boot), c(n_boot, n_Ta, length(n_subject)))
  wjs_vec = t(wjs[1, , ])
  avg_U_hat_star = apply(sqrt(wjs) * U_hat_star, c(1, 2), sum)
  avg_U_hat_star_karray = array(rep(avg_U_hat_star, times = length(n_subject)), c(n_boot, n_Ta, length(n_subject)))
  U2_boot = apply(wjs * (U_hat_star / sqrt(wjs) - avg_U_hat_star_karray) ^ 2, c(1, 2), sum)  # bootstrap SSB(a) in Theorem 1, n_boot x n_Ta
  ### Bootstrap sup_{a \in [alpha_1, alpha_2]} SSB(a) in Theorem 3:
  sup_boot = apply(as.matrix(U2_boot), 1, max)  
  sup_EL_crit = as.vector(quantile(sup_boot, probs = 1 - alpha_vec)) 
  ### EL test:
  teststat_pre = 1:n_Ta * 0 
  error_vec = 1:n_Ta * 0
  for (ai in 1:n_Ta) {
    neg2logRa = neg2logR_test(ai, Ta, n_subject, error_out = 1) 
    error_vec[ai] = neg2logRa$still_error
    if (error_vec[ai] == 1) next 
    teststat_pre[ai]=neg2logRa$neg2logR_crit
  }  # END for
  suptest = max(teststat_pre)
  inttest_pre_da = teststat_pre[-n_Ta] * diff(as)
  inttest_da = sum(inttest_pre_da) 
  return(list(
    suptest = suptest,
    sup_EL_crit = sup_EL_crit,
    out_sup_pval = mean(sup_boot >= suptest)
  ))
}  # END procedure


neg2logR_test_logger = function(ai, Ta, n_subject, error_out = 0) {  
  # Computes -2logR(a)
  # Args:
  #   ai: index for the given activity level a
  #   Ta: a list with Ta[[j]] being the matrix of activity profile data, i-th row representing the i-th subject's activity profile
  #   n_subject: a vector, n_subject[j] being the sample size for the j-th group
  #   error_out: whether more detailed information is to be returned, including whether there is an error,
  # Returns:
  #   -2logR(a) when error_out = 0; when error_out = 1, an addition output
  #   still_error giving 1 if there is an error in computing -2logR(a), 0 otherwise
  n = sum(n_subject)
  gamma_js = n_subject / sum(n_subject)
  n_group = length(n_subject)
  min_ups_max_lows = matrix(sapply(1:n_group, FUN = function (j) {
    range((Ta[[j]])[, ai])
  }), nrow = 2)  # 1st row: min of each group; 2nd row: max of each group
  ### Dealing with non-overlapping regions
  
  min_ups_max_lows_logs[[ai]] <<- min_ups_max_lows
  
  if (max(min_ups_max_lows[1, ]) >= min(min_ups_max_lows[2, ])) {
    if (error_out == 0) {
      return(Inf)
    } else {
      return(list(neg2logR_crit = Inf, 
                  still_error = 0   
      ))
    }
  } 
  ### Get T_{aj} from the Ta list for given a and j, and scale it so that max range differs by 1
  T_js = lapply(1:n_group, FUN = function (j) {   
    (Ta[[j]])[, ai] / (min(min_ups_max_lows[2, ]) - max(min_ups_max_lows[1, ]))  
  })
  
  T_js_logs[[ai]] <<- T_js
  
  min_ups_max_lows = matrix(sapply(1:n_group, FUN = function (j) {
    range(T_js[[j]])
  }), nrow = 2)  # 1st row: min of each group; 2nd row: max of each group
  
  min_ups_max_lows_logs2[[ai]] <<- min_ups_max_lows
  
  n_grid_inc = 0  # the number of times we have updated the number of grid points for mu and lambda initialization
  still_error = 1  # 1 if there is an error in computing -2logR(a), 0 otherwise
  test1t = Inf  # if no solution below, then means sup{empty set} = 0 in numerator of -2logR
  n_mu_grid = parameters$n_mu_grid  # number of grid points for mu initialization
  n_lamb_grid = parameters$n_lamb_grid  # number of grid points for lambda initialization
  ### A while loop for computing -2logR(a), allowing n_grid_inc to be updated if still_error = 1
  while (still_error == 1 & n_grid_inc <= parameters$n_grid_incbd) {  
    n_grid_inc = n_grid_inc + 1  
    n_mu_grid = n_mu_grid * n_grid_inc + 1  
    n_lamb_grid = n_lamb_grid * n_grid_inc + 1  
    ### Grid construction for mu and lambda initialization
    mu_tol = (min(min_ups_max_lows[2, ]) - max(min_ups_max_lows[1, ])) / parameters$mu_tol_grid
    init_mu_a_grid = seq(max(min_ups_max_lows[1, ]) + mu_tol, min(min_ups_max_lows[2, ]) - mu_tol, length = n_mu_grid)
    
    init_mu_a_grid_logs[[ai]] <<- list()
    init_mu_a_grid_logs[[ai]][[n_grid_inc]] <<- init_mu_a_grid
    
    ll_l = matrix(0, nrow = n_mu_grid, ncol = n_group)
    ul_l = matrix(0, nrow = n_mu_grid, ncol = n_group)
    init_Delta_grid = array(unlist(sapply(1:n_mu_grid, FUN = function (init_mu_ind) {
      init_mu_a = init_mu_a_grid[init_mu_ind]
      sapply(1:n_group, FUN = function (j) {  
        g_j = (T_js[[j]] - init_mu_a) / gamma_js[j]  # a vector of length n_j
        ll_l[init_mu_ind, j] <<- max((1/n_subject[j] - 1) / g_j[g_j > 0])
        ul_l[init_mu_ind, j] <<- min((1/n_subject[j] - 1) / g_j[g_j < 0])
        lamb_tol = (ul_l[init_mu_ind, j] - ll_l[init_mu_ind, j]) / parameters$lamb_tol_grid
        seq(ll_l[init_mu_ind, j] + lamb_tol, ul_l[init_mu_ind, j] - lamb_tol, length = n_lamb_grid) 
      })
    })), c(n_lamb_grid, n_group, n_mu_grid))
    
    ll_logs[[ai]] <<- list()
    ll_logs[[ai]][[n_grid_inc]] <<- ll_l
    
    ul_logs[[ai]] <<- list()
    ul_logs[[ai]][[n_grid_inc]] <<- ul_l
    
    init_Delta_grid_logs[[ai]] <<- list()
    init_Delta_grid_logs[[ai]][[n_grid_inc]] <<- init_Delta_grid
    
    ### Mu and lambda initialization 
    init_check = matrix(0, nrow = n_mu_grid, ncol = n_lamb_grid * (n_lamb_grid ^ (n_group - 2)))
    init_Delta_grid_expand = matrix(unlist(lapply(1:n_mu_grid, FUN = function (init_mu_ind) {
      init_mu_a = init_mu_a_grid[init_mu_ind]
      
      mat1 = cbind(init_mu_a, 
                   sapply(1:(n_group - 1), FUN = function (j) {
                     rep(
                       rep(init_Delta_grid[, j, init_mu_ind], each = n_lamb_grid ^ (j - 1)), 
                       length = n_lamb_grid * (n_lamb_grid ^ (n_group - 2)))
                     })
                   )
      
      init_lamb_km1 = -apply(as.matrix(mat1[, -1]), 1, sum)
      
      init_check[init_mu_ind, ] <<- (init_lamb_km1 > ll_l[init_mu_ind, n_group] & init_lamb_km1 < ul_l[init_mu_ind, n_group])
      
      return(t(mat1[(init_check[init_mu_ind, ] ==1), ]))
    })), byrow = T, ncol = n_group)
    
    init_Delta_grid_expand_logs[[ai]] <<- list()
    init_Delta_grid_expand_logs[[ai]][[n_grid_inc]] <<- init_Delta_grid_expand
    
    
    init_Delta_grid_expand_no_filter = matrix(unlist(lapply(1:n_mu_grid, FUN = function (init_mu_ind) {
      init_mu_a = init_mu_a_grid[init_mu_ind]
      
      mat1 = cbind(init_mu_a, sapply(1:(n_group - 1), FUN = function (j) {
        rep(rep(init_Delta_grid[, j, init_mu_ind],each = n_lamb_grid ^ (j - 1)), length = n_lamb_grid * (n_lamb_grid ^ (n_group - 2)))
      }))
      
      mat1_check <<- mat1
      
      init_lamb_km1 = -apply(as.matrix(mat1[, -1]), 1, sum)
      
      init_lamb_km1_check <<- init_lamb_km1
      
      init_check[init_mu_ind, ] <<- (init_lamb_km1 > ll_l[init_mu_ind, n_group] & init_lamb_km1 < ul_l[init_mu_ind, n_group])
      
      init_check_check <<- init_check
      
      ul_check[[ai]] <<- list()
      ll_check[[ai]] <<- list()
      ul_check[[ai]][[n_grid_inc]] <<- list()
      ll_check[[ai]][[n_grid_inc]] <<- list()
      ul_check[[ai]][[n_grid_inc]][[init_mu_ind]] <<- ul_l[init_mu_ind, n_group]
      ll_check[[ai]][[n_grid_inc]][[init_mu_ind]] <<- ll_l[init_mu_ind, n_group]
      
      return(t(mat1))
    })), byrow = T, ncol = n_group)
    
    init_Delta_grid_expand_no_filter_logs[[ai]] <<- list()
    init_Delta_grid_expand_no_filter_logs[[ai]][[n_grid_inc]] <<- init_Delta_grid_expand_no_filter
    
    ### Find gamma and lambda satisfying the estimating equations related to optimization in -2logR(a) 
    n_init = dim(init_Delta_grid_expand)[1]
    status_grid = matrix(NA, nrow = n_init, ncol = 2 + 3 * n_group - 1)  
    roots = matrix(NA, nrow = n_init, ncol = 2 * n_group - 1) 
    multir = list()  
    multir$estim.precis = 1  
    min_ps = rep(-1, times = n_group)  
    for (init_ind in 1:n_init) {
      init_mu_a = init_Delta_grid_expand[init_ind, 1]

      init_lambda_js = cumsum(-init_Delta_grid_expand[init_ind, -1])
      init_gam = gamma_js - diff(c(0, init_lambda_js, 0)) * (-init_mu_a)
      lamb_gam = c(init_lambda_js, init_gam)
      
      init_lambda_kp1 = c(0, init_lambda_js, 0)  # j-th element is (j-1)-th element of lambdas
      p = list()
      init_ps = sapply(1:n_group, FUN = function (j) {
        p[[j]] <<- 1 / (n * (init_gam[j] + (init_lambda_kp1[j] - init_lambda_kp1[j + 1]) * T_js[[j]]))
      })
      
      init_model = model(lamb_gam, T_js, n_subject)
      
      init_model_logs[[ai]] <<- list()
      init_model_logs[[ai]][[n_grid_inc]] <<- list()
      init_model_logs[[ai]][[n_grid_inc]][[init_ind]] <<- init_model
      
      init_p_logs[[ai]] <<- list()
      init_p_logs[[ai]][[n_grid_inc]] <<- list()
      init_p_logs[[ai]][[n_grid_inc]][[init_ind]] <<- init_ps
      
      init_param_logs[[ai]] <<- list()
      init_param_logs[[ai]][[n_grid_inc]] <<- list()
      init_param_logs[[ai]][[n_grid_inc]][[init_ind]] <<- c(init_mu_a, init_lambda_js, init_gam)
      
      # options(warn = 2)	 # warning = error! will stop the whole function
      possibleError <- tryCatch(
        multir<-multiroot(f = model, start = lamb_gam, lbs = T_js, n_subject = n_subject), 
        error=function(e) e
      )
      # options(warn = 0)
      
      root_logs[[ai]] <<- list()
      root_logs[[ai]][[n_grid_inc]] <<- list()
      root_logs[[ai]][[n_grid_inc]][[init_ind]] <<- possibleError
      
      was_error[[ai]] <<- list()
      was_error[[ai]][[n_grid_inc]] <<- list()
      was_error[[ai]][[n_grid_inc]][[init_ind]] <<- !inherits(possibleError, "error")
      
      if (!inherits(possibleError, "error")) {  # !inherits(possibleError, "error") == T if there's no error
        roots[init_ind, ] = multir$root
        lamb = roots[init_ind, ][1:(n_group - 1)]
        gam = roots[init_ind, ][-(1:(n_group - 1))]
        
        gam_logs[[ai]] <<- list()
        gam_logs[[ai]][[n_grid_inc]] <<- list()
        gam_logs[[ai]][[n_grid_inc]][[init_ind]] <<- gam
        
        lamb_logs[[ai]] <<- list()
        lamb_logs[[ai]][[n_grid_inc]] <<- list()
        lamb_logs[[ai]][[n_grid_inc]][[init_ind]] <<- lamb
        
        lamb_kp1 = c(0, lamb, 0)  # j-th element is (j-1)-th element of lambdas
        p = list()
        min_ps = sapply(1:n_group, FUN = function (j) {
          p[[j]] <<- 1 / (n * (gam[j] + (lamb_kp1[j] - lamb_kp1[j + 1]) * T_js[[j]]))
          min(p[[j]])
        })
        status_grid[init_ind, ] = c(!inherits(possibleError, "error"), multir$estim.precis, multir$f.root, min_ps)
        
        p_logs[[ai]] <<- list()
        p_logs[[ai]][[n_grid_inc]] <<- list()
        p_logs[[ai]][[n_grid_inc]][[paste0(n_grid_inc, "_")]] <<- p
        
        min_p_logs[[ai]] <<- list()
        min_p_logs[[ai]][[n_grid_inc]] <<- list()
        min_p_logs[[ai]][[n_grid_inc]][[paste0(n_grid_inc, "_")]] <<- min_ps
        
        print(paste0(init_ind, " ", init_mu_a))
        print(paste0("No Error in Root?:  ", !inherits(possibleError, "error")))
        print(paste0("Precision Met: ", multir$estim.precis < 10 ^ (-8)))
        print(paste0("Number negative ps: ", sum(min_ps < 0)))
        print("")
        
        if (multir$estim.precis < 10 ^ (-8) & !inherits(possibleError, "error") & sum(min_ps < 0) == 0) break  
      }  else {
        print("possibleError had an error!")
      }
      
      
      
    }  # END for
    
    status_grid_logs[[ai]] <<- list()
    status_grid_logs[[ai]][[n_grid_inc]] <<- status_grid
    
    still_error = !(multir$estim.precis < 10 ^ (-8) & !inherits(possibleError, "error") & sum(min_ps < 0) == 0)
  }  # END while  
  ### Compute -2logR(a) given the p output above
  if (exists("p") == 1) {  
    test1t = -2 * sum(unlist(lapply(1:n_group, FUN = function (j) {
      return(log(p[[j]]))
    }))) - 2 * sum(n_subject * log(n_subject))
  }  # END if
  if (error_out == 0) {
    return(test1t)
  } else {
    return(list(neg2logR_crit = test1t, 
                still_error = still_error
    ))
  }  # END if/else
}  # END neg2logR_CB
