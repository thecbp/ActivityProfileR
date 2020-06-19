#' Simulate activity data for a groups of subjects
#'
#' Given parameters for each group, simulates activity data for each group.
#'
#' @param group_params A list of parameters for each group in the simulated sample
#' @param quantiles A vector of two percentages describing how which quantiles of activity to keep
#'
#' @return A tibble of the activity for the sample
#' @export
#'
#' @importFrom magrittr %>%
#' @name %>%
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
#' 
prepare_sim_data = function(group_params, 
                            quantiles = c(0.05, 0.95)) {
  
  n_groups = length(group_params)
  data = NULL
  
  # Simulate data based on parameters for different groups
  for (i in 1:n_groups) {
    data = dplyr::bind_rows(
      data,
      generate_data(n_subj = group_params[[i]]$n_subj, 
                    group_id = i, 
                    n_points = group_params[[i]]$n_points, 
                    scaling = group_params[[i]]$scaling, 
                    thetas = group_params[[i]]$thetas)
    )
  }
  
  data
}
