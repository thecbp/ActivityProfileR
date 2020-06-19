#' Scales activity profiles so that max range is 1
#'
#' @param data dataset containing the activity profiles
#' @param a activity index that describes where in the activity profiles to filter
#'
#' @return data with another column containing the scaled activity profiles
#' @export
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @name %>%
#'
scale_activity_profile = function(data, 
                                  a,
                                  group_col) {
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
  
  # Quoting the variables so we don't need to quote the variables in code
  group_col = dplyr::enquo(group_col)
  
  # Get the minimum and maximum from the activity profiles in each group for a given ai
  range_by_group_at_ai = data %>% 
    dplyr::filter(.data$ai == a) %>% 
    dplyr::group_by(!!group_col) %>% 
    dplyr::summarize(
      Ta_min = min(.data$Ta),
      Ta_max = max(.data$Ta)
    )
  
  # Log the extremes of the activity profiles among groups
  minmax_at_ai = range_by_group_at_ai %>% dplyr::pull(.data$Ta_max) %>% min
  maxmin_at_ai = range_by_group_at_ai %>% dplyr::pull(.data$Ta_min) %>% max
  
  # Scale the activity profiles at a to have range 1
  data %>% 
    dplyr::filter(.data$ai == a) %>% 
    dplyr::mutate(
      scaled_Ta = .data$Ta / (minmax_at_ai - maxmin_at_ai)
    ) %>% 
    dplyr::select(.data$group, .data$n, .data$gamma, .data$scaled_Ta)
}