################ 20181124 age group analysis: can we make other tests' results less significant by
################ considering less balanced groups??
dir_path_source = "E:\\R_program"
setwd(dir_path_source)
# source model function for root finding
source("mean_EL_k_20180623.R")
# source k-sample test procedure
# (if reinstall R, need to install the packages all over again)
source("paper_fnl_simu_k_20180504_EE_simpleloop_EP_withair_source.R")  # parameters in the following is in there!!
parameters$n_agrid = 436  # updated 20190110 # 100  # updated 20181001
dir_path2 = "E:\\R_program\\fnl_mil_20190103"

# 20181128: load new_demog_data directly from paper_fnl_NHANES_data_analysis_activity_20180619_subgroups.R results
load("new_mydata.Rdata")
#range(ridageyr[!is.na(dmqmilit) & dmqmilit == 1])
#[1]  17 85
#range(ridageyr[!is.na(dmqmilit) & dmqmilit == 2])
#[1] 17 85
agej_ub = 10
n_agemilgroups = matrix(0, nrow = agej_ub + 1, ncol = 4)
n_miss_agemilgroups = matrix(0, nrow = agej_ub + 1, ncol = 4)
n_group = 4

agej = 0
time1 = Sys.time()
#for (agej in 0:agej_ub) {
setwd(dir_path_source)
load("new_demog_data.Rdata")
# before gender restriction: too significant... SH suggested restricting to males only (I: decreasing sample size might help)
attach(new_demog_data)  # seqn at this point is valid, mydata and demog_data available


band_select_groups = list()
# don't restrict to military, and first look at age & gender groups, using some age cutoff to restrict sample size?
attach(new_demog_data)  # seqn at this point is valid, mydata and demog_data available
band_select_groups = list()
#band_select_groups[[1]] = seqn[(!is.na(ridageyr) & ridageyr > 75 - 5* agej & ridageyr <= 85 - 5* agej & !is.na(dmqmilit) & dmqmilit == 1)] # & !is.na(riagendr) & riagendr == 1
#band_select_groups[[2]] = seqn[(!is.na(ridageyr) & ridageyr > 75 - 5* agej & ridageyr <= 85 - 5* agej & !is.na(dmqmilit) & dmqmilit == 2)] # & !is.na(riagendr) & riagendr == 1
#band_select_groups[[3]] = seqn[(!is.na(ridageyr) & ridageyr > 65 - 5* agej & ridageyr <= 75 - 5* agej & !is.na(dmqmilit) & dmqmilit == 1)] # & !is.na(riagendr) & riagendr == 1
#band_select_groups[[4]] = seqn[(!is.na(ridageyr) & ridageyr > 65 - 5* agej & ridageyr <= 75 - 5* agej & !is.na(dmqmilit) & dmqmilit == 2)] # & !is.na(riagendr) & riagendr == 1
# successful group=c(2,1) eg:
#band_select_groups[[1]] = seqn[(!is.na(ridageyr) & ridageyr > 75 - 5* agej & ridageyr <= 85 - 5* agej & !is.na(dmqmilit) & dmqmilit == 1)] # & !is.na(riagendr) & riagendr == 1
#band_select_groups[[2]] = seqn[(!is.na(ridageyr) & ridageyr > 75 - 5* agej & ridageyr <= 85 - 5* agej & !is.na(dmqmilit) & dmqmilit == 2)] # & !is.na(riagendr) & riagendr == 1
#band_select_groups[[3]] = seqn[(!is.na(ridageyr) & ridageyr >= 65 - 5* agej & ridageyr <= 75 - 5* agej & !is.na(dmqmilit) & dmqmilit == 1)] # & !is.na(riagendr) & riagendr == 1
#band_select_groups[[4]] = seqn[(!is.na(ridageyr) & ridageyr >= 65 - 5* agej & ridageyr <= 75 - 5* agej & !is.na(dmqmilit) & dmqmilit == 2)] # & !is.na(riagendr) & riagendr == 1
# recall previous results
band_select_groups[[1]] = seqn[(!is.na(ridageyr) & ridageyr >= 75 - 5* agej & ridageyr <= 85 - 5* agej & !is.na(dmqmilit) & dmqmilit == 1)] # & !is.na(riagendr) & riagendr == 1
band_select_groups[[2]] = seqn[(!is.na(ridageyr) & ridageyr >= 75 - 5* agej & ridageyr <= 85 - 5* agej & !is.na(dmqmilit) & dmqmilit == 2)] # & !is.na(riagendr) & riagendr == 1
band_select_groups[[3]] = seqn[(!is.na(ridageyr) & ridageyr >= 65 - 5* agej & ridageyr < 75 - 5* agej & !is.na(dmqmilit) & dmqmilit == 1)] # & !is.na(riagendr) & riagendr == 1
band_select_groups[[4]] = seqn[(!is.na(ridageyr) & ridageyr >= 65 - 5* agej & ridageyr < 75 - 5* agej & !is.na(dmqmilit) & dmqmilit == 2)] # & !is.na(riagendr) & riagendr == 1
#

# note: can't run the following after doing 1-sample band analysis, bec the source file there will
# replace some of the parameters$ settings such as mu_hat_bounds and n_agrid
    Xt = list()
    sample_indx = list()
    n_group = 4
    n_subject = 1:n_group * 0
    outside_endpts = matrix(0, nrow = 2, ncol = n_group)
    new_mydata_band_select = list()  # 20181118
 
    as_all_ranges = sapply(1:n_group, FUN = function (j) {
      new_mydata_band_select[[j]] <<- new_mydata[new_mydata[, 1] %in% band_select_groups[[j]], ]
      attach(new_mydata_band_select[[j]])
      subject_ids = unique(seqn)
      n_subject[j] <<- length(unique(seqn))
#      n_subject[j] = length(unique(seqn))
      t_ncol = length(unique(paxn))
      Xt[[j]] <<- lapply(1:n_subject[j], FUN = function (i) {
#     Xt[[j]] = lapply(1:n_subject[j],  FUN = function (i) {
        subject_id = subject_ids[i]
#        paxinten[seqn == subject_id]
        paxinten[seqn == subject_id & paxstat==1 & paxcal==1]  # added on 20181118: post-valid trt
      })
      # assume don't know population... so don't build as from original Xt_list population!
      as_all = sort(unique(unlist(Xt[[j]])))  
      #nNA_i = 1:n_subject[j] * 0
      #time1 = Sys.time()
      Ta_all = t(sapply(1:n_subject[j], FUN = function (i) {
        observed_Xt = (Xt[[j]])[[i]][!is.na((Xt[[j]])[[i]])]
        observed_Xt_len = length(observed_Xt)
        # 1 / observed_Xt_len: adjusting for missing data (scratch paper 20180718 p.1L; from original gridsol = 1 / 10080)
        grid_width_i = rep(1 / observed_Xt_len, observed_Xt_len)  # 20180813 for k-sample comparison, use gridsol = 1 / 10080
        # otherwise will have to do t_ncol / observed_Xt_len
        Ta_i = mu_a_hat(observed_Xt, grid_widths = grid_width_i, as_all)  #to-do: correct grid_width_i in data anlaysis><
        #nNA_i[i] <<- sum(is.na(Ta_i))
        #if (nNA_i[i] == 0) return(Ta_i)
      }))  
      #time2 = Sys.time()
      #time2 - time1
      mu_hat_all = apply(Ta_all, 2, mean)
      ind_keep = (mu_hat_all >= parameters$mu_hat_bounds[1] & mu_hat_all <= parameters$mu_hat_bounds[2])
      outside_endpts[, j] <<- c(mean(mu_hat_all < parameters$mu_hat_bounds[1]), mean(mu_hat_all > parameters$mu_hat_bounds[2]))  # proportions of data falling outside
#     outside_endpts[, j] = c(mean(mu_hat_all < parameters$mu_hat_bounds[1]), mean(mu_hat_all > parameters$mu_hat_bounds[2]))  # proportions of data falling outside
      range(as_all[ind_keep])
    })
    # the following numbers - 2 and + 7 are by looking at resulting as and adjust...
    as = seq(max(as_all_ranges[1, ]), min(as_all_ranges[2, ]), length = parameters$n_agrid)  # look at as first, then see what grid
n_subject = 1:n_group * 0
    observed_Xt_lens = list()

    Ta = lapply(1:n_group, FUN = function (j) {
      observed_Xt_lens_tmp = 1:n_subject[j] * 0
      new_mydata_band_select = new_mydata[new_mydata[, 1] %in% band_select_groups[[j]], ]
      attach(new_mydata_band_select)  # due to day first, redundant for eah day
      subject_ids = unique(seqn)  # due to day first, redundant for eah day
      n_subject[j] <<- length(unique(seqn))  # due to day first, redundant for eah day
      
      Ta_j = t(sapply(1:n_subject[j], FUN = function (i) {
        subject_id = subject_ids[i]
        observed_Xt = paxinten[seqn == subject_id & paxstat==1 & paxcal==1]
        observed_Xt_len = length(observed_Xt)
        observed_Xt_lens_tmp[i] <<- observed_Xt_len
        # 1 / observed_Xt_len: adjusting for missing data
        grid_width_i = rep(1 / observed_Xt_len, observed_Xt_len)  # each observed data (number = length((Xt[[j]])[[i]]) 
        mu_a_hat(observed_Xt, grid_widths = grid_width_i, as)
      }))
      observed_Xt_lens[[j]] <<- observed_Xt_lens_tmp
      return(Ta_j)
    })


  for (missj in 1:n_group){
    n_miss_agemilgroups[agej + 1, missj] = sum(observed_Xt_lens[[missj]] != 10080)
  }

  n_subject = sapply(1:n_group, FUN = function (j) {
    dim(Ta[[j]])[1]
  })
  n_agemilgroups[agej + 1, ] = n_subject

  Ta_mat = matrix(unlist(lapply(1:length(n_subject), FUN = function (j) {  # for fda.usc's fdata object construction
    c(t(Ta[[j]]))
  })), byrow = T, ncol = parameters$n_agrid)
anova_onefactor_group = as.factor(rep(1:length(n_subject), times = n_subject))

setwd(dir_path2)

png(paste("mean", agej, sep=''),width = 7, height = 6, units='in', res = 600)
require(graphics)

# 20181128 mean and variance plot!
mean_n_1 = apply(Ta_mat[1:n_subject[1],], 2, mean)
mean_n_2 = apply(Ta_mat[(n_subject[1] + 1):(n_subject[1] + n_subject[2]),], 2, mean)
mean_n_3 = apply(Ta_mat[(n_subject[1] + n_subject[2] + 1):(n_subject[1] + n_subject[2] + n_subject[3]),], 2, mean)
mean_n_4 = apply(Ta_mat[(n_subject[1] + n_subject[2] + n_subject[3] + 1):(n_subject[1] + n_subject[2] + n_subject[3] + n_subject[4]),], 2, mean)
par(mfrow=c(1, 1))
plot(as, mean_n_1, type = "l", ylim = c(min(mean_n_1, mean_n_2, mean_n_3, mean_n_4), max(mean_n_1, mean_n_2, mean_n_3, mean_n_4)), ylab = 'mean')
lines(as, mean_n_2, type = "l", col = 'red')
lines(as, mean_n_3, type = "l", col = 'blue')
lines(as, mean_n_4, type = "l", col = 'green')

dev.off()

png(paste("var", agej, sep=''),width = 7, height = 6, units='in', res = 600)
require(graphics)

# 20181125: investigating variance / sample proportion i.e. weight
var_n_1 = apply(Ta_mat[1:n_subject[1],], 2, var) 
var_n_2 = apply(Ta_mat[(n_subject[1] + 1):(n_subject[1] + n_subject[2]),], 2, var) 
var_n_3 = apply(Ta_mat[(n_subject[1] + n_subject[2] + 1):(n_subject[1] + n_subject[2] + n_subject[3]),], 2, var) 
var_n_4 = apply(Ta_mat[(n_subject[1] + n_subject[2] + n_subject[3] + 1):(n_subject[1] + n_subject[2] + n_subject[3] + n_subject[4]),], 2, var) #par(mfrow=c(1, 1))
plot(as, var_n_1, type = "l", ylim = c(min(var_n_1, var_n_2, var_n_3, var_n_4), max(var_n_1, var_n_2, var_n_3, var_n_4)), ylab = 'var')
lines(as, var_n_2, type = "l", col = 'red')
lines(as, var_n_3, type = "l", col = 'blue')
lines(as, var_n_4, type = "l", col = 'green')

dev.off()

print(agej)

#}

time2 = Sys.time()
time2 - time1
n_agemilgroups
      [,1] [,2] [,3] [,4]
 [1,]  160  279  139  348

n_miss_agemilgroups
      [,1] [,2] [,3] [,4]
 [1,]    0    5    0    2

# missing rate calculation:
j=1
sum(observed_Xt_lens[[j]] != 10080) # 0
j=2
sum(observed_Xt_lens[[j]] != 10080) # 5
j=3
sum(observed_Xt_lens[[j]] != 10080) # 0
j=4
sum(observed_Xt_lens[[j]] != 10080) # 2
j=2
1-(observed_Xt_lens[[j]])[observed_Xt_lens[[j]] != 10080]/10080
[1] 0.07787698 0.35138889 0.05327381 0.05228175 0.47857143
j=4
1-(observed_Xt_lens[[j]])[observed_Xt_lens[[j]] != 10080]/10080
[1] 0.0546627 0.6229167
# so, ranges from 5--62%!


####*** traditional analysis using Xt? see military20 file for conclusion...
# later if asked to do Xt analysis in this dataset, can just compare the two 
# age groups among veterans!!

if (0) {
## analysis

####done
# to check our run time
setwd(dir_path_source)
source("paper_fnl_simu_k_20180504_EE_simpleloop_EP_withair_source_oursonly.R")  # parameters in the following is in there!!

parameters$n_subject = n_subject
set.seed(916)
time1 = Sys.time()
    test_out = procedure(Ta = Ta, as = as, n_boot = parameters$nboot, alpha_vec = parameters$alpha_vec)
time2 = Sys.time()
time2 - time1 #Time difference of 3.992503 mins


test_out$out_sup_pval  
#[1] 0
test_out$out_dF_pval
#[1] 0
test_out$out_sup_EP_pval
#[1] 0
test_out$out_dF_EP_pval
#[1] 0


####done
library(fdANOVA)
setwd(dir_path_source)
source("fanova_tests_20181204.R") # this magically removes the pkg warnings!! 
set.seed(1012)
# check time: run tests individually, found that in smaller number of time points
# lit tests work fine!
time1 = Sys.time()
(fanova1 <- fanova.tests(x = t(Ta_mat), group.label = anova_onefactor_group, test = c("CS")))
#     Analysis of Variance for Functional Data  
#CS test - L2-norm-based parametric bootstrap test for heteroscedastic samples 
#Test statistic = 444.8666  p-value = 0 
time2 = Sys.time()
time2 - time1
#Time difference of 12.78221 secs
time1 = Sys.time()
(fanova1 <- fanova.tests(x = t(Ta_mat), group.label = anova_onefactor_group, test = c("GPF")))
#     Analysis of Variance for Functional Data 
#GPF test - globalizing the pointwise F-test 
#Test statistic = 22.94267  p-value = 5.551115e-16 
time2 = Sys.time()
time2 - time1
#Time difference of 0.500036 secs
time1 = Sys.time()
(fanova1 <- fanova.tests(x = t(Ta_mat), group.label = anova_onefactor_group, test = c("L2B")))
#     Analysis of Variance for Functional Data 
#L2B test - L2-norm-based test with bias-reduced method of estimation 
#Test statistic = 126.0921  p-value = 1.326717e-13 
time2 = Sys.time()
time2 - time1
#Time difference of 0.3125229 secs
time1 = Sys.time()
(fanova1 <- fanova.tests(x = t(Ta_mat), group.label = anova_onefactor_group, test = c("FB")))
#     Analysis of Variance for Functional Data 
#FB test - F-type test with bias-reduced method of estimation 
#Test statistic = 19.24144  p-value = 3.550493e-13 
time2 = Sys.time()
time2 - time1
#Time difference of 0.281271 secs
time1 = Sys.time()
(fanova1 <- fanova.tests(x = t(Ta_mat), group.label = anova_onefactor_group, test = c("FP")))
#     Analysis of Variance for Functional Data 
#FP test - permutation test based on a basis function representation 
#Test statistic = 19.24147  p-value = 0 
time2 = Sys.time()
time2 - time1
#Time difference of 2.404868 mins
time1 = Sys.time()
(fanova1 <- fanova.tests(x = t(Ta_mat), group.label = anova_onefactor_group, test = c("Fmaxb")))
#     Analysis of Variance for Functional Data 
#Fmaxb test - Fmax bootstrap test 
#Test statistic = 32.98689  p-value = 0 
time2 = Sys.time()
time2 - time1
#Time difference of 17.88734 mins
time1 = Sys.time()
(fanova1 <- fanova.tests(x = t(Ta_mat), group.label = anova_onefactor_group, test = c("TRP")))
#     Analysis of Variance for Functional Data 
#TRP - tests based on k = 30 random projections 
#p-value ANOVA = 0 
#p-value ATS = 0 
#p-value WTPS = 0 
time2 = Sys.time()
time2 - time1
Time difference of 1.168577 mins


### done
# same as above
# other existing tests! see 'paper_fnl_data_addothers_20181102.R'
set.seed(1012)
(fanova1 <- fanova.tests(x = t(Ta_mat), group.label = anova_onefactor_group, test = c("CS", "GPF", "L2B", "FB", "FP", "Fmaxb", "TRP")))


### done
#### going to paper_fnl_simu_k_20180504_EE_simpleloop_EP_withair_source_oursonly.R line by line
#### to see where the diff in K_n is biggest:
#### before that:
Ta_orig = Ta
#### then run procedure's argument in the following group_sub = c(1, 2)
which.max(teststat_pre)
round(teststat_pre, 2)
sup_ptEL_crit = apply(as.matrix(U2_boot), 2, quantile, probs = 1 - alpha_vec)
round(sup_ptEL_crit, 2)
# can see sup_ptEL_crit roughly similar across as, whereas teststat_pre increasing as a increases
#### after this
Ta = Ta_orig


#### done
group_sub = c(1, 2)
parameters$n_subject = n_subject[group_sub]
set.seed(916)
time1 = Sys.time()
    test_out = procedure(Ta = Ta[group_sub], as = as, n_boot = parameters$nboot, alpha_vec = parameters$alpha_vec)
time2 = Sys.time()
time2 - time1 #Time difference of 27.72084 secs for c(3, 4)


#group_sub = c(1, 2), c(3, 4), c(1, 3), c(2, 4): 
# >= 75 for groups 1 & 2
test_out$out_sup_pval  
#[1] 0.01, 0.345, 0, 0
test_out$out_dF_pval
#[1] 0.056, 0.385, 0, 0
test_out$out_sup_EP_pval
#[1] 0.013, 0.343, 0, 0
test_out$out_dF_EP_pval
#[1] 0.062, 0.385, 0, 0
test_out$out_anova_onefactor_pval
#[1] 0.119,

#### done
#group_sub = c(1, 2), c(3, 4), c(1, 3), c(2, 4):
#>= 75 for groups 1 & 2
set.seed(1012)
Ta_mat_sub = rbind(Ta_mat[anova_onefactor_group == group_sub[1], ], Ta_mat[anova_onefactor_group == group_sub[2], ])
anova_onefactor_group_sub = c(anova_onefactor_group[anova_onefactor_group %in% group_sub[1]], anova_onefactor_group[anova_onefactor_group %in% group_sub[2]])
(fanova1 <- fanova.tests(x = t(Ta_mat_sub), group.label = anova_onefactor_group_sub, test = c("CS", "GPF", "L2B", "FB", "FP", "Fmaxb", "TRP")))
 

#####detailed results
#>= 75 for groups 2 & 4 
     Analysis of Variance for Functional Data 
 
FP test - permutation test based on a basis function representation 
Test statistic = 47.20625  p-value = 0 
 
CS test - L2-norm-based parametric bootstrap test for heteroscedastic samples 
Test statistic = 192.0041  p-value = 0 
 
L2B test - L2-norm-based test with bias-reduced method of estimation 
Test statistic = 106.5668  p-value = 5.487832e-13 
 
FB test - F-type test with bias-reduced method of estimation 
Test statistic = 47.20614  p-value = 1.570744e-12 
 
GPF test - globalizing the pointwise F-test 
Test statistic = 57.09654  p-value = 3.219647e-15 
 
Fmaxb test - Fmax bootstrap test 
Test statistic = 82.12102  p-value = 0 
 
TRP - tests based on k = 30 random projections 
p-value ANOVA = 0 
p-value ATS = 0 
p-value WTPS = 0 

#>= 75 for groups 1 & 3 
     Analysis of Variance for Functional Data 
 
FP test - permutation test based on a basis function representation 
Test statistic = 9.197982  p-value = 0.002 
 
CS test - L2-norm-based parametric bootstrap test for heteroscedastic samples 
Test statistic = 40.17594  p-value = 0.0018 
 
L2B test - L2-norm-based test with bias-reduced method of estimation 
Test statistic = 18.67711  p-value = 0.001681491 
 
FB test - F-type test with bias-reduced method of estimation 
Test statistic = 9.197967  p-value = 0.001902662 
 
GPF test - globalizing the pointwise F-test 
Test statistic = 10.91561  p-value = 0.0006692706 
 
Fmaxb test - Fmax bootstrap test 
Test statistic = 16.63647  p-value = 1e-04 
 
TRP - tests based on k = 30 random projections 
p-value ANOVA = 4.574612e-05 
p-value ATS = 4.744049e-05 
p-value WTPS = 0 

#>= 75 for groups 3 & 4 
     Analysis of Variance for Functional Data 
 
FP test - permutation test based on a basis function representation 
Test statistic = 0.896886  p-value = 0.349 
 
CS test - L2-norm-based parametric bootstrap test for heteroscedastic samples 
Test statistic = 3.064773  p-value = 0.3254 
 
L2B test - L2-norm-based test with bias-reduced method of estimation 
Test statistic = 2.190022  p-value = 0.3508002 
 
FB test - F-type test with bias-reduced method of estimation 
Test statistic = 0.8968883  p-value = 0.3522181 
 
GPF test - globalizing the pointwise F-test 
Test statistic = 0.6946989  p-value = 0.4164821 
 
Fmaxb test - Fmax bootstrap test 
Test statistic = 1.498893  p-value = 0.3646 
 
TRP - tests based on k = 30 random projections 
p-value ANOVA = 0.5872636 
p-value ATS = 0.5766456 
p-value WTPS = 0.5785 

#>= 75 for groups 1 & 2 
     Analysis of Variance for Functional Data 
 
FP test - permutation test based on a basis function representation 
Test statistic = 2.31906  p-value = 0.104 
 
CS test - L2-norm-based parametric bootstrap test for heteroscedastic samples 
Test statistic = 6.928297  p-value = 0.1224 
 
L2B test - L2-norm-based test with bias-reduced method of estimation 
Test statistic = 4.403178  p-value = 0.1232201 
 
FB test - F-type test with bias-reduced method of estimation 
Test statistic = 2.319055  p-value = 0.1247247 
 
GPF test - globalizing the pointwise F-test 
Test statistic = 3.40863  p-value = 0.06048778 
 
Fmaxb test - Fmax bootstrap test 
Test statistic = 7.4881  p-value = 0.0164 
 
TRP - tests based on k = 30 random projections 
p-value ANOVA = 0.03537065 
p-value ATS = 0.0387567 
p-value WTPS = 0.033

}

###done
# source 1-sample band procedure (have to run k-sample above first for parameters input)
#source("paper_fnl_data_20180602_cts_fixagrid_source.R")
setwd(dir_path_source)
source("paper_fnl_data_20180602_cts_fixagrid_source_addothers.R")  # add EL into Functions in 'paper_fnl_data_simu_forNHANES_20180602_Xt2_cts_fixagrid_addothers_20180903.R'

time1 = Sys.time()
# 20181002 set.seed(916) 
CB_out = lapply(1:n_group, FUN = function (j) {
set.seed(916 + j)  # because there's bootstrap
  procedure(Ta = Ta[[j]], as = as, n_boot = 1000, mu_a = NA, alpha_vec = c(0.05, 0.01))
})
time2 = Sys.time()
time2 - time1 #Time difference of 57.95883 secs
# 20181001 after adding MFD band Time difference of 12.81953 mins

#summary = matrix(0, nrow = 3, ncol = 3)
summary = matrix(0, nrow = 4, ncol = 4)
for (j in 1:4){
  HW_widths = ((CB_out[[j]])$HW_CB)[2, ]-((CB_out[[j]])$HW_CB)[1, ]
  EP_widths = ((CB_out[[j]])$EP_CB)[2, ]-((CB_out[[j]])$EP_CB)[1, ]
  EL_widths = ((CB_out[[j]])$EL_CB)[2, ]-((CB_out[[j]])$EL_CB)[1, ]
  EP_lit_widths = ((CB_out[[j]])$EP_CB_lit)[2, ]-((CB_out[[j]])$EP_CB_lit)[1, ]
  #summary[j, ] = apply(cbind(EL_widths, EP_widths, HW_widths), 2, mean)
  summary[j, ] = apply(cbind(EL_widths, EP_widths, HW_widths, EP_lit_widths), 2, mean) 
}
summary
#EL_widths EP_widths HW_widths EP_lit_widths
           [,1]       [,2]       [,3]       [,4]
[1,] 0.02282743 0.02254470 0.03232025 0.02256009
[2,] 0.01694894 0.01688620 0.02857629 0.01661452
[3,] 0.02729052 0.02696995 0.03919893 0.02519303
[4,] 0.01823673 0.01815081 0.02568813 0.01741610


  j=1
  HW_widths = ((CB_out[[j]])$HW_CB)[2, ]-((CB_out[[j]])$HW_CB)[1, ]
  EP_widths = ((CB_out[[j]])$EP_CB)[2, ]-((CB_out[[j]])$EP_CB)[1, ]
  EL_widths = ((CB_out[[j]])$EL_CB)[2, ]-((CB_out[[j]])$EL_CB)[1, ]
  EP_lit_widths = ((CB_out[[j]])$EP_CB_lit)[2, ]-((CB_out[[j]])$EP_CB_lit)[1, ]

  sum(HW_widths < EL_widths)
#15
(HW_widths - EL_widths)[HW_widths < EL_widths]
range((HW_widths - EL_widths)[HW_widths >= EL_widths])
#[1] 0.0001673772 0.0163048551
  sum(EP_lit_widths < EL_widths)
#436
(EP_lit_widths - EL_widths)[EP_lit_widths < EL_widths]
range((EP_lit_widths - EL_widths)[EP_lit_widths < EL_widths])
#[1] -0.0003297889 -0.0001349805  # diff small compared w HW vs EL


# 20190110: looking at a=99
(a_indx = which.min(abs(as-99))), (a_indx = which.min(abs(as-400)))
#[1] 100, 401
as[100]

  j = 1
# sendentary is time spent <= 99; the rest of time is spent > 99
# we look at the latter, so to get the former we can do 1 - below
# then to get the hours/day, we do (1 - below) * 24
#  HW_widths = 
((CB_out[[j]])$HW_CB)[2, a_indx]
#[1] 0.1788697, [1] 0.08299792
((CB_out[[j]])$HW_CB)[1, a_indx]
#[1] 0.1465494, [1] 0.05067767
#  EL_widths = 
((CB_out[[j]])$EL_CB)[2, a_indx]
#[1] 0.1764038, [1] 0.07586545
((CB_out[[j]])$EL_CB)[1, a_indx]
#[1] 0.1500423, [1] 0.05894621
#  EP_lit_widths = 
((CB_out[[j]])$EP_CB_lit)[2, a_indx]
#  ub.norm 
#0.1757301 , 0.07521322 
((CB_out[[j]])$EP_CB_lit)[1, a_indx]
#  lb.norm 
#0.1496854 , 0.05846131 

#  HW
((CB_out[[j]])$HW_CB)[2, a_indx] * 24
#[1] 1.99195
((CB_out[[j]])$HW_CB)[1, a_indx] * 24
#[1] 1.216264
#  EL 
((CB_out[[j]])$EL_CB)[2, a_indx] * 24
#[1] 1.820771
((CB_out[[j]])$EL_CB)[1, a_indx] * 24
#[1] 1.414709
#  EP 
((CB_out[[j]])$EP_CB_lit)[2, a_indx] * 24
# ub.norm 
#1.805117 
((CB_out[[j]])$EP_CB_lit)[1, a_indx] * 24
# lb.norm 
#1.403072 


################20190109 plotting 2005--2006 NHANES activity profiles & mean est + mean est add CB (from 20181114)
# for the first group: male with military service
dir_path2 = paste("C:\\Users\\", Sys.getenv("USERNAME"), "\\Dropbox\\papers_mine\\paper_fnl_data\\figures", sep = "")
setwd(dir_path2)
left_mar=6
bot_mar=4
n_quant = 5

#0823 to-do: thicken the lines, and do n_quant = 10 or 5, then check in latex

png("paper_fnl_fig_NHANES_CBs_20190109.png",width = 16, height = 6, units='in', res = 600)
require(graphics)

par(mfrow=c(1,3),mar=c(bot_mar,left_mar,1,1.1),oma=c(0.3,0.5,0,0),las=1,mgp=c(3,1,0), cex.axis = 1.5, font.lab = 50, tck = -0.015) #oma: outer margin; mar: inner margin

# 1. activity profiles & mean est


plot(as, (Ta[[1]])[1, ], type = "n", col = 'gray', ylim = c(min(unlist(Ta)), max(unlist(Ta))), xlab = '', ylab = '')
#mtext(expression(a),cex=1.4,line=2.5, side=1 )
mtext(expression(Ta),cex=1.3,line=3.5, side=2 )
mtext('activity level (a)',cex=1.3,line=3, side=1 )

j = 1
rep_subject = 1:(n_quant + 1) * 0
for (quanti in 0:n_quant) {
 quant = quantile((Ta[[j]])[, 1], probs=quanti * (1 / n_quant))
 rep_subject[quanti + 1] = which.min(abs((Ta[[j]])[, 1] - quant))
 points(as, (Ta[[j]])[rep_subject[quanti + 1], ], type="l", lwd = 1.2) #, col = 'gray')
}

# As commented by the Diggle et al book (see Diggle et al pp33-39), plotting all the raw data would be too busy
#  for (i in 1:n_subject[j]) {
#    points(as, (Ta[[j]])[i, ], type="l", col = 'gray')
#  }

j = 2
rep_subject = 1:(n_quant + 1) * 0
for (quanti in 0:n_quant) {
 quant = quantile((Ta[[j]])[, 1], probs=quanti * (1 / n_quant))
 rep_subject[quanti + 1] = which.min(abs((Ta[[j]])[, 1] - quant))
 points(as, (Ta[[j]])[rep_subject[quanti + 1], ], type="l", col = "red", lwd = 1.2) #, lty = 2) #col = 'lightblue') # 
}

j = 3
rep_subject = 1:(n_quant + 1) * 0
for (quanti in 0:n_quant) {
 quant = quantile((Ta[[j]])[, 1], probs=quanti * (1 / n_quant))
 rep_subject[quanti + 1] = which.min(abs((Ta[[j]])[, 1] - quant))
 points(as, (Ta[[j]])[rep_subject[quanti + 1], ], type="l", col = "blue", lwd = 1.2) #, col = 'gray') #lty=6#
}

j = 4
rep_subject = 1:(n_quant + 1) * 0
for (quanti in 0:n_quant) {
 quant = quantile((Ta[[j]])[, 1], probs=quanti * (1 / n_quant))
 rep_subject[quanti + 1] = which.min(abs((Ta[[j]])[, 1] - quant))
 points(as, (Ta[[j]])[rep_subject[quanti + 1], ], type="l", col = "green", lwd = 1.2) #lty = 5 longdash #, col = 'gray') #col = 'pink') #
}

#points(as, (CB_out[[1]])$mu_hat_vec, type="l")
#points(as, (CB_out[[2]])$mu_hat_vec, type="l", lty = 2) #, col = 'blue') #
#points(as, (CB_out[[3]])$mu_hat_vec, type="l", lty = 6) #, col = 'red') # 

# 2. mean est add CB for each group

plot(as, ((CB_out[[1]])$EL_CB)[1, ], type = "l", ylim = c(min(c((CB_out[[1]])$EL_CB, (CB_out[[2]])$EL_CB, (CB_out[[3]])$EL_CB)), max(c((CB_out[[1]])$EL_CB, (CB_out[[2]])$EL_CB, (CB_out[[3]])$EL_CB))), xlab = '', ylab = '')
#mtext(expression(a),cex=1.4,line=2.5, side=1 )
mtext(expression(mu(a)),cex=1.3,line=3.6, side=2 )
mtext('activity level (a)',cex=1.3,line=3, side=1 )

j = 1
points(as, (CB_out[[j]])$mu_hat_vec, type="l", lty = 2)
points(as, ((CB_out[[j]])$EL_CB)[1, ], type="l")
points(as, ((CB_out[[j]])$EL_CB)[2, ], type="l")

j = 2
points(as, (CB_out[[j]])$mu_hat_vec, type="l", col = "red", lty = 2)
points(as, ((CB_out[[j]])$EL_CB)[1, ], type="l", col = "red") #lty = 2
points(as, ((CB_out[[j]])$EL_CB)[2, ], type="l", col = "red")

j = 3
points(as, (CB_out[[j]])$mu_hat_vec, type="l", col = "blue", lty = 2)
points(as, ((CB_out[[j]])$EL_CB)[1, ], type="l", col = "blue") #lty=6
points(as, ((CB_out[[j]])$EL_CB)[2, ], type="l", col = "blue")


j = 4
points(as, (CB_out[[j]])$mu_hat_vec, type="l", col = "green", lty = 2)
points(as, ((CB_out[[j]])$EL_CB)[1, ], type="l", col="green") # lty = 5
points(as, ((CB_out[[j]])$EL_CB)[2, ], type="l", col="green")


# 3. mean est add CB

j = 1
plot(as, ((CB_out[[1]])$EL_CB)[1, ], type = "n", ylim = c(min(c((CB_out[[1]])$EL_CB, (CB_out[[2]])$EL_CB, (CB_out[[3]])$EL_CB)), max(c((CB_out[[1]])$EL_CB, (CB_out[[2]])$EL_CB, (CB_out[[3]])$EL_CB))), xlab = '', ylab = '')
#plot(as, ((CB_out[[j]])$EL_CB)[1, ], type = "l", lty = 3, ylim = c(min(c((CB_out[[j]])$EL_CB, (CB_out[[j]])$EP_CB, (CB_out[[j]])$HW_CB)), max(c((CB_out[[j]])$EL_CB, (CB_out[[j]])$EP_CB, (CB_out[[j]])$HW_CB))), xlab = '', ylab = '')
mtext('activity level (a)',cex=1.3,line=3, side=1 )
mtext(expression(mu(a)),cex=1.3,line=3.6, side=2 )

points(as, ((CB_out[[j]])$EL_CB)[1, ], type="l", lty=1)
points(as, ((CB_out[[j]])$EL_CB)[2, ], type="l", lty=1)
points(as, ((CB_out[[j]])$EP_CB_lit)[1, ], type="l", col="blue")  #lty=6
points(as, ((CB_out[[j]])$EP_CB_lit)[2, ], type="l", col="blue")
points(as, ((CB_out[[j]])$HW_CB)[1, ], type="l", col="red") #lty = 2
points(as, ((CB_out[[j]])$HW_CB)[2, ], type="l", col="red")
points(as, (CB_out[[j]])$mu_hat_vec, type="l", col = 'gray', lwd = 1.7)


dev.off()
