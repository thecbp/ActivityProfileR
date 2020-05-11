#' Add an activity profile column to the dataset
#'
#' Takes all of the subject's activity counts, Xt, and transforms it
#' into the activity profile, Ta
#'
#' @param data dataset containing the activity counts for each group in the sample
#' @param activity_grid A vector derived from extract_activity_grid
#'
#' @return the same dataset with an additional columns for the activity profiles and analysis
#' @export
#'
#' @examples
#' 
#' activity_grid = extract_activity_grid(data)
#' processed_data = add_activity_profile(data, activity_grid)
#' 
add_activity_profile = function(data, activity_grid) {

  create_activity_profile = function(Xt, activity_grid) {
    Ta = tibble::tibble(
      ai = activity_grid
    ) %>% 
      dplyr::mutate(
        Ta = unlist(map(ai, function(a) { 
          mean(Xt > a)
        }))
      ) %>% 
      pull(Ta)
    
    return(Ta)
  }
  
  processed_data = data %>% 
    dplyr::mutate(
      Ta = purrr::map(Xt, ~create_activity_profile(.x, activity_grid = activity_grid)),
      ai = purrr::map(group, function(x) { return(1:length(activity_grid)) })
    ) %>% 
    # These additions are also present in prepare_sim_data, but adding here in case user brings their own data
    dplyr::group_by(group) %>% 
    dplyr::mutate(
      n = n()
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      gamma = n / nrow(.) # indicates the proportion each group contributes to sample
    )
  
  return(processed_data)
}
