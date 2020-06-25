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
#'
extract_activity_grid = function(data, 
                                 activity_col,
                                 group_col, 
                                 grid_length = 100, 
                                 quantiles = c(0.05, 0.95))  {
  
  # Quoting the variables so we don't need to quote the variables in code
  activity_col = dplyr::enquo(activity_col)
  group_col = dplyr::enquo(group_col)
  
  # Get the 5th and 95th quantiles of activity counts for each group
  activity_quantiles_by_group = data %>% 
    tidyr::unnest(!!activity_col) %>%
    dplyr::group_by(!!group_col) %>% 
    dplyr::summarize(
      qt5 = quantile(!!activity_col, quantiles[1]),
      qt95 = quantile(!!activity_col, quantiles[2])
    )
  
  # Filter out activity count data that is outside of the quantiles
  # Then get the new range of activity counts for each group
  filtered_ranges_by_group = dplyr::left_join(data, 
                                              activity_quantiles_by_group, 
                                              by = rlang::as_name(group_col)) %>% 
    tidyr::unnest(!!activity_col) %>% 
    dplyr::filter(
      !!activity_col >= .data$qt5, 
      !!activity_col < .data$qt95) %>% 
    dplyr::group_by(!!group_col) %>% 
    dplyr::summarize(
      min_activity_val = min(!!activity_col),
      max_activity_val = max(!!activity_col)
    )
  
  # Extract the minmax and maxmin to create the activity grid
  maxmin = filtered_ranges_by_group %>% 
    dplyr::pull(.data$min_activity_val) %>% max
  minmax = filtered_ranges_by_group %>% 
    dplyr::pull(.data$max_activity_val) %>% min
  
  seq(maxmin, minmax, length = grid_length)
}