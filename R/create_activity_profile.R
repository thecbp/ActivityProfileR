#' Creates activity profile from activity counts
#'
#' @param Xt vector of activity counts for a subject
#' @param activity_grid vector ofactivity counts derived from extract_activity_grid
#'
#' @return a vector representing the activity profile for the data
#' @export
#'
#' @examples
create_activity_profile = function(Xt, activity_grid) {
  Ta = tibble::tibble(
    ai = activity_grid
  ) %>% 
    mutate(
      Ta = unlist(map(ai, function(a) { 
        mean(Xt > a)
      }))
    ) %>% 
    pull(Ta)
  
  return(Ta)
}