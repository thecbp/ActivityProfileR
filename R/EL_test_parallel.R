#' Paralllelized version of the EL_test function
#'
#' @param data dataset containing the activity profiles
#' @param grid_length number describing how long the activity grid is
#' @param n_boot number describing how many bootstrap datasets to make
#' @param quantiles two percentiles describing what percentiles of activity points to keep in grid
#' @param alpha number describing what quantile to check for the boostrapped test statistics
#' @param n_cores number describing how many cores to use for the parallelization
#'
#' @return a list contaning the results of the functional ANOVA test
#' @export
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @name %>%
#'
#' @examples
#' par_el = EL_test_parallel(sim_data, n_cores = n_cores)
#' 
EL_test_parallel = function(data, 
                            activity_col,
                            group_col,
                            grid_length = 100, 
                            n_boot = 1000, 
                            quantiles = c(0.05, 0.95), 
                            alpha = 0.95,
                            n_cores = 2) {

  cluster = multidplyr::new_cluster(n_cores)
  multidplyr::cluster_library(cluster, "tidyverse")
  multidplyr::cluster_library(cluster, "purrr")
  multidplyr::cluster_library(cluster, "rootSolve")
  multidplyr::cluster_copy(cluster, "scale_activity_profile")
  multidplyr::cluster_copy(cluster, "neg2logRa")
  multidplyr::cluster_copy(cluster, "initialize_grid")
  multidplyr::cluster_copy(cluster, "ee_model")
  
  # Add the activity profiles for each subject
  processed_data = data %>% 
    add_activity_profile(.data, 
                         activity_col = activity_col,
                         group_col = group_col,
                         quantiles = quantiles,
                         grid_length = grid_length) %>% 
    add_helper_columns(group_col = group_col) %>% 
    tidyr::unnest(c(.data$Ta, .data$ai))
  
  multidplyr::cluster_assign(cluster, processed_data = processed_data)
  
  # Calculate the -2logR(a) statistic for each activity index
  neg2logRa_tests = tibble::tibble(
    ai = 1:grid_length,
    par_group = rep(1:n_cores, length.out = grid_length) 
  ) %>% 
    dplyr::group_by(.data$par_group) %>% 
    multidplyr::partition(cluster) %>% 
    dplyr::mutate(
      # confirmed correct: resulting scaled Ta is sensitive to the range
      # which has been affected by the beta random vars
      scaled_Ta_at_ai = purrr::map(.data$ai, 
                                   ~scale_activity_profile(processed_data, 
                                                           a = .x)), 
      neg2logRa = purrr::map2(.data$scaled_Ta_at_ai, 
                              .data$ai, 
                              ~neg2logRa(data = .x, a = .y))
    ) %>% 
    dplyr::collect() %>% 
    dplyr::pull(neg2logRa) %>% 
    unlist()
  
  # Bootstrap the EL statistic through the uniform approximation U^2
  bs = bootstrap_U_star(processed_data)
  sup_boot = bs %>% dplyr::pull(sup_boot)
  sup_EL_crit = stats::quantile(sup_boot, alpha)

  # Capture the supremum of -2logR(a) across all activity indices
  sup_test = max(neg2logRa_tests)
  
  return(list(
    sup_test = sup_test,
    sup_EL_crit = sup_EL_crit,
    out_sup_pval = mean(sup_boot >= sup_test),
    neg2logRa_tests = neg2logRa_tests,
    bootstrap = bs
  ))
  
}