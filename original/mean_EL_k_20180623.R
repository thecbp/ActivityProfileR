# Computes log EL for k means, based on max_usual_EL_20160510.R
model = function(lamb_gam, lbs, n_subject, maxeval = 10000) {
# x: a vector of 2 * n_group - 1 length: with c(lambdas, gammas)
# lbs: Ta[[j]][, ai] for j = 1, ..., n_group
# n_subject replacing gr
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