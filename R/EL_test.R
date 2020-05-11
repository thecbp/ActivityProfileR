#' Calculates functional ANOVA test on activity profiles
#'
#' @param data dataset containing the activity profiles
#' @param grid_length number describing how long the activity grid is
#' @param n_boot number describing how many bootstrap datasets to make
#' @param quantiles two percentiles describing what percentiles of activity points to keep in grid
#' @param alpha number describing what quantile to check for the boostrapped test statistics
#' @param verbose TRUE/FALSE describing if you want to see the progresss of the functions
#'
#' @return a list contaning the results of the functional ANOVA test
#' @export
#'
#' @examples
#' el = EL_test(sim_data)
#' 
EL_test = function(data, 
                   grid_length = 100, 
                   n_boot = 1000, 
                   quantiles = c(0.05, 0.95), 
                   alpha = 0.95,
                   verbose = FALSE) {
  # Calculates the empirical likelihood test statistic of the data
  # Data assumed to be nested at first
  
  # Create an activity grid based off of the 
  activity_grid = extract_activity_grid(data, 
                                        quantiles = quantiles,
                                        grid_length = grid_length)
  
  processed_data = data %>% 
    dplyr::mutate(
      Ta = purrr::map(Xt, ~create_activity_profile(.x, activity_grid = activity_grid)), 
      ai = purrr::map(group, function(x) { return(1:length(activity_grid)) })
    ) %>% 
    dplyr::group_by(group) %>% 
    dplyr::mutate(
      n = n()
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      # Add a gamma to indicate the proportion each group contributes to data
      gamma = n / nrow(.)
    ) %>% 
    dplyr::select(-Xt)
  
  unnested_processed_data =  processed_data %>% tidyr::unnest(c(Ta, ai))
  
  # Calculate the -2logR(a) statistic for each activity index
  neg2logRa_tests = tibble::tibble(
    ai = 1:grid_length
  ) %>% 
    dplyr::mutate(
      # confirmed correct: resulting scaled Ta is sensitive to the range
      # which has been affected by the beta random vars
      scaled_Ta_at_ai = purrr::map(ai, ~scale_activity_profile(unnested_processed_data, a = .x)), 
      neg2logRa = map2(scaled_Ta_at_ai, ai, ~neg2logRa(data = .x, 
                                                       a = .y, 
                                                       verbose = verbose))
    ) %>% 
    dplyr::pull(neg2logRa) %>% unlist()
  
  # Bootstrap the EL statistic through the uniform approximation U^2
  bs = bootstrap_U_star(processed_data)
  sup_boot = bs %>% dplyr::pull(sup_boot)
  sup_EL_crit = quantile(sup_boot, alpha)

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