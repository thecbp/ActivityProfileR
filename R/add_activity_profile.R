#' Add an activity profile column to the dataset
#'
#' Takes all of the subject's activity counts, Xt, and transforms it
#' into the activity profile, Ta
#'
#' @param data dataset containing the activity counts for each group in the sample
#' @param activity_col the name of the column that contains the activity data
#' @param group_col the name of the column that contains the grouping label
#'
#' @return the same dataset with an additional columns for the activity profiles and analysis
#' @export
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
add_activity_profile = function(data, 
                                activity_col, 
                                group_col,
                                grid_length = 100, 
                                quantiles = c(0.05, 0.95)) {
  
  # Quoting the variables so we don't need to quote the variables in code
  activity_col = dplyr::enquo(activity_col)
  group_col = dplyr::enquo(group_col)
  
  activity_grid = extract_activity_grid(data, 
                                        activity_col = !!activity_col,
                                        group_col = !!group_col,
                                        grid_length = grid_length,
                                        quantiles = quantiles)
  
  data %>% 
    dplyr::mutate(
      Ta = purrr::map(!!activity_col, 
                      ~create_activity_profile(.x, 
                                               activity_grid = activity_grid)),
      ai = purrr::map(!!group_col, function(x) { 
        tibble::tibble( ai = 1:grid_length )
        })
    )

}
