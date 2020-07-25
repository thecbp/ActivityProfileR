#' Calculates functional ANOVA test on activity profiles
#'
#' @param data dataset containing the activity profiles
#' @param id_col the column that contains the subject ID column
#' @param group_col the column that contains the group ID column
#' @param activity_col the column that contains the activity counts column
#' @param grid_length number describing how long the activity grid is
#' @param n_boot number describing how many bootstrap datasets to make
#' @param quantiles two percentiles describing what percentiles of activity points to keep in grid
#' @param alpha number describing what quantile to check for the boostrapped test statistics
#' @param verbose TRUE/FALSE describing if you want to see the progresss of the functions
#'
#' @return a list contaning the results of the functional ANOVA test
#' @export
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
#' @examples
#' params_by_group = list(
#'   list(n_subj = 70, n_points = 1000, 
#'       scaling = 300, thetas = c(-0.8, 2.0, 2.31)),
#'  list(n_subj = 100, n_points = 1000, 
#'       scaling = 300, thetas = c(-0.402, 1.0, 1.65)),
#'  list(n_subj = 130, n_points = 1000, 
#'       scaling = 300, thetas = c(-0.201, 0.5, 1.18))
#'  )
#' 
#' data = prepare_sim_data(params_by_group)
#' el = EL_test(data, id_col = subject_id, group_col = group, activity_col = Xt)
#' 
EL_test = function(data, 
                   id_col,
                   group_col,
                   activity_col,
                   grid_length = 100, 
                   n_boot = 1000, 
                   quantiles = c(0.05, 0.95), 
                   alpha = 0.95,
                   verbose = FALSE) {
  
  # Quoting the variables so we don't need to quote the variables in code
  id_col = dplyr::enquo(id_col)
  group_col = dplyr::enquo(group_col)
  activity_col = dplyr::enquo(activity_col)
  
  # Add the activity profiles for each subject
  processed_data = data %>% 
    dplyr::select(-!!activity_col) %>% # Remove original activity counts
    tidyr::unnest(c(.data$Ta, .data$ai))

  # Calculate the -2logR(a) statistic for each activity index
  neg2logRa_tests = tibble::tibble(
    ai = 1:grid_length
  ) %>% 
    dplyr::mutate(
      # activity profiles at a given activity index should be scaled
      # such that the max range is 1 
      scaled_Ta_at_ai = purrr::map(.data$ai, 
                                   ~scale_activity_profile(processed_data, 
                                                           a = .x,
                                                           group_col = !!group_col)), 
      # Need to calculate the neg2logR(a) for each activity index
      # which we need to take the supremum over
      neg2logRa = purrr::map2_dbl(.data$scaled_Ta_at_ai, 
                                  .data$ai, 
                                  ~neg2logRa(data = .x, 
                                             a = .y, 
                                             group_col = !!group_col,
                                             verbose = verbose))
    ) %>% 
    dplyr::pull(neg2logRa)
  
  # Capture the supremum of -2logR(a) across all activity indices
  sup_test = max(neg2logRa_tests)
  
  # Bootstrap the EL statistic through the uniform approximation U^2
  bs = bootstrap_U_star(processed_data, 
                        id_col = !!id_col, 
                        group_col = !!group_col)
  sup_boot = bs %>% dplyr::pull(sup_boot)
  sup_EL_crit = quantile(sup_boot, alpha)
  
  list(
    sup_test = sup_test,
    sup_EL_crit = sup_EL_crit,
    out_sup_pval = mean(sup_boot >= sup_test),
    neg2logRa_tests = neg2logRa_tests,
    bootstrap = bs
  )
  
}