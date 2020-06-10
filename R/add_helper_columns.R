#' Add columns that are helpful in calculation of EL test statistic
#'
#' @param data tibble containing the activity data and profiles
#' @param group_col string containing the name of the column with the group 
#' label
#'
#' @return tibble containing extra helper columns "gamma" and "n"
#' @export
#'
#' @examples
#' processed_data = data %>% 
#'   add_activity_profile(activity_col = "Xt", group_col = "group") %>% 
#'   add_helper_columns(group_col = "group")
add_helper_columns = function(data, group_col) {
  
  data %>% 
    dplyr::group_by(.data[[group_col]]) %>% 
    dplyr::mutate(
      n = dplyr::n()
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      gamma = .data$n / nrow(.data) # indicates the proportion each group contributes to sample
    )
}