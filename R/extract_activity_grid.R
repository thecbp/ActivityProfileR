#' Create a grid of activity values to standardize activity profiles
#' 
#' Takes a dataset containing activity counts for multiple groups
#' and creates a standard grid/vector of activity counts that can 
#' be applied across all of the groups in the dataset. Activity counts
#' assumed to be nested. See output of prepare_sim_data
#'
#' @param data dataset containing groups and activity 
#' @param grid_length number describing how long do we want the activity grid to be
#' @param quantiles two percentiles describing what percentiles of activity points to keep in grid
#'
#' @return a vector of activity values
#' @export
#'
#' @examples
#' activity_grid = extract_activity_grid(data)
extract_activity_grid = function(data, 
                                 grid_length = 100, 
                                 quantiles = c(0.05, 0.95))  {
  # Given a matrix of activity data, output a standardized  equidistant vector of  
  # activity values that is compatible with all values in the data
  
  # Arguments:
  # data: a matrix of activity data,
  #         each row represents the activity data of one subject
  #         should be match the structure of the generate.data function output
  # grid_length: the desired length of the returned vector of activity values
  # quantiles: vector of 2 floats between 0 and 1 describing the 
  #              range of quanitles of activity data that you want kept in the grid
  
  # Returns:
  # a vector of equidistant activity values from the maximum to the minimum in the data
  
  # Get the 5th and 95th quantiles of activity counts for each group
  Xt_quantiles_by_group = data %>% 
    tidyr::unnest(Xt) %>%
    dplyr::group_by(group) %>% 
    dplyr::summarize(
      qt5 = quantile(Xt, quantiles[1]),
      qt95 = quantile(Xt, quantiles[2])
    )
  
  # Filter out activity count data that is outside of the quantiles
  # Then get the new range of activity counts for each group
  filtered_ranges_by_group = dplyr::left_join(data, Xt_quantiles_by_group, by = "group") %>% 
    tidyr::unnest(Xt) %>% 
    dplyr::filter(Xt >= qt5, Xt < qt95) %>% 
    dplyr::group_by(group) %>% 
    dplyr::summarize(
      min_Xt = min(Xt),
      max_Xt = max(Xt)
    )
  
  # Now extract the minmax and maxmin from these ranges to create the activity grid
  maxmin = max(filtered_ranges_by_group$min_Xt)
  minmax = min(filtered_ranges_by_group$max_Xt)
  
  seq(maxmin, minmax, length = grid_length)
}