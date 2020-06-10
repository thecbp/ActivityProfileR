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
#' @name %>%
#'
#' @examples
#' processed_data = data %>% 
#'   add_activity_profile(activity_col = "Xt", group_col = "group")
#' 
add_activity_profile = function(data, 
                                activity_col, 
                                group_col,
                                grid_length = 100, 
                                quantiles = c(0.05, 0.95)) {
  
  activity_grid = extract_activity_grid(data, 
                                        activity_col = activity_col,
                                        group_col = group_col,
                                        grid_length = grid_length,
                                        quantiles = quantiles)
  
  data %>% 
    dplyr::mutate(
      Ta = purrr::map(.data[[activity_col]], 
                      ~create_activity_profile(.x, 
                                               activity_grid = activity_grid)),
      ai = purrr::map(.data[[group_col]], function(x) { 
        tibble::tibble( ai = 1:length(activity_grid) )
        })
    )

}
