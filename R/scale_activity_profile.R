#' Scales activity profiles so that max range is 1
#'
#' @param data dataset containing the activity profiles
#' @param a activity index that describes where in the activity profiles to filter
#'
#' @return data with another column containing the scaled activity profiles
#' @export
#'
#' @examples
scale_activity_profile = function(data, a) {
  # Scale the activity profiles so that max range among profiles is 1
  # at a given activity index a
  
  # Arguments: 
  # data: tibble that contains information on group membership and 
  #       data on the activity profiles for each subject
  # a: activity index that describes where in the activity profiles we want to
  #    filter to (should be between 1 and the the length oof the activity grid)
  
  # Returns:
  # the data with another column containing the scaled activity profiles
  
  # NEED TO DO:
  # Create a if-else statement here to catch instances where there is no overlap between groups
  # if (max(min_ups_max_lows[1, ]) >= min(min_ups_max_lows[2, ])) {
  #   if (error_out == 0) {
  #     return(Inf)
  #   } else {
  #     return(list(neg2logR_crit = Inf, 
  #                 still_error = 0   
  #     ))
  #   }
  # } 
  
  
  # Get the minimum and maximum from the activity profiles in each group for a given ai
  range_by_group_at_ai = data %>% 
    dplyr::filter(ai == a) %>% 
    dplyr::group_by(group) %>% 
    dplyr::summarize(
      Ta_min = min(Ta),
      Ta_max = max(Ta)
    )
  
  # Log the extremes of the activity profiles among groups
  minmax_at_ai = range_by_group_at_ai %>% dplyr::pull(Ta_max) %>% min()
  maxmin_at_ai = range_by_group_at_ai %>% dplyr::pull(Ta_min) %>% max()
  
  # Scale the activity profiles at a to have range 1
  scaled_data_at_ai = data %>% 
    dplyr::filter(ai == a) %>% 
    dplyr::mutate(
      scaled_Ta = Ta / (minmax_at_ai - maxmin_at_ai)
    ) %>% 
    dplyr::select(group, n, gamma, scaled_Ta)
  
  return(scaled_data_at_ai)
}