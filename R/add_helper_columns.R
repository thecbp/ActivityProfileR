#' Add columns that are helpful in calculation of EL test statistic
#'
#' @param data tibble containing the activity data and profiles
#' @param group_col string containing the name of the column with the group 
#' label
#'
#' @return tibble containing extra helper columns "gamma" and "n"
#' @export
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
add_helper_columns = function(data, group_col) {
  
  total_rows = nrow(data)
  group_col = dplyr::enquo(group_col)
  
  data %>% 
    dplyr::group_by(!!group_col) %>% 
    dplyr::mutate(
      n = dplyr::n()
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      gamma = .data$n / total_rows # indicates the proportion each group contributes to sample
    )
}