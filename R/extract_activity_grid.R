#' Given a matrix of activity data, output a standardized equidistant vector of  
#' activity values that is compatible with all values in the data
#' 
#' Takes a dataset containing activity counts for multiple groups
#' and creates a standard grid/vector of activity counts that can 
#' be applied across all of the groups in the dataset. Activity counts
#' assumed to be nested. See output of prepare_sim_data
#'
#' @param data dataset containing groups and activity 
#' @param activity_col the name of the column that contains the activity data
#' @param group_col the name of the column that contains the grouping label
#' @param grid_length number describing how long the activity grid should be
#' @param quantiles two percentiles describing what percentiles of activity 
#' points to keep in grid
#'
#' @return a vector of equidistant activity values from the maximum to the 
#' minimum in the data
#' @export
#'
#' @importFrom rlang .data
#' @importFrom stats quantile
#' @importFrom magrittr %>%
#' @name %>%
#'
#' @examples
#' #' params_by_group = list(
#'   list(n_subj = 70, n_points = 1000, 
#'       scaling = 300, thetas = c(-0.8, 2.0, 2.31)),
#'  list(n_subj = 100, n_points = 1000, 
#'       scaling = 300, thetas = c(-0.402, 1.0, 1.65)),
#'  list(n_subj = 130, n_points = 1000, 
#'       scaling = 300, thetas = c(-0.201, 0.5, 1.18))
#'  )
#'  
#' data = prepare_sim_data(params_by_group)
#' activity_grid = extract_activity_grid(data, activity_col = "Xt", group_col = "group")
extract_activity_grid = function(data, 
                                 activity_col,
                                 group_col, 
                                 grid_length = 100, 
                                 quantiles = c(0.05, 0.95))  {
  
  # Get the 5th and 95th quantiles of activity counts for each group
  activity_quantiles_by_group = data %>% 
    tidyr::unnest(.data[[activity_col]]) %>%
    dplyr::group_by(.data[[group_col]]) %>% 
    dplyr::summarize(
      qt5 = quantile(.data[[activity_col]], quantiles[1]),
      qt95 = quantile(.data[[activity_col]], quantiles[2])
    )
  
  # Filter out activity count data that is outside of the quantiles
  # Then get the new range of activity counts for each group
  filtered_ranges_by_group = dplyr::left_join(data, 
                                              activity_quantiles_by_group, 
                                              by = group_col) %>% 
    tidyr::unnest(.data[[activity_col]]) %>% 
    dplyr::filter(.data[[activity_col]] >= .data$qt5, 
                  .data[[activity_col]] < .data$qt95) %>% 
    dplyr::group_by(.data[[group_col]]) %>% 
    dplyr::summarize(
      min_activity_val = min(.data[[activity_col]]),
      max_activity_val = max(.data[[activity_col]])
    )
  
  # Now extract the minmax and maxmin from these ranges to create the activity grid
  maxmin = max(filtered_ranges_by_group$min_activity_val)
  minmax = min(filtered_ranges_by_group$max_activity_val)
  
  seq(maxmin, minmax, length = grid_length)
}