#' Creates activity profile from activity counts
#'
#' @param Xt vector of activity counts for a subject
#' @param activity_grid vector ofactivity counts derived from extract_activity_grid
#'
#' @return a vector representing the activity profile for the data
#' @export
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @name %>%
#'
create_activity_profile = function(activity, activity_grid) {
  
  Ta = tibble::tibble(
    ai = activity_grid
  ) %>% 
    dplyr::mutate(
      Ta = unlist(purrr::map(.data$ai, function(a) { 
        mean(activity > a)
      }))
    ) %>% 
    dplyr::select(Ta)
  
  return(Ta)
}