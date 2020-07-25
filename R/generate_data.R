#' Generate activity data
#'
#' Simulates data for a group of subjects using an Ornstein-Uhlenbeck process
#' with user defined parameters for the process
#'
#' @param n_subj number describing how many subjects are in this particular group
#' @param group_id number describing the index of the group to differentiate it from others
#' @param n_points how much activity count data we want to produce for each person
#' @param scaling how much the activity data should be scaled
#' @param thetas numbers for parameters in creating Ornstein-Uhlenbeck process
#'
#' @return data containing the activity data
#' @export
#' 
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' 
#' @examples
#' data = generate_data(n_subj = 70, group_id = 1, n_points = 1000, 
#'                      scaling = 300, thetas = c(-0.8, 2.0, 2.31))
generate_data = function(n_subj, group_id, n_points, scaling, thetas) {
  
  # Set up tibble for group
  data = tibble::tibble(
    subject_id = 1:n_subj,
    group = rep(group_id, n_subj)
  )
  
  # Simulate a tibble with each row being a subject's activity data
  activity_data = ESGtoolkit::simdiff(n = n_subj, 
                                      horizon = n_points, 
                                      frequency = "annual", 
                                      model = "OU", 
                                      x0 = 0, 
                                      theta1 = thetas[1], 
                                      theta2 = thetas[2], 
                                      theta3 = thetas[3])
  
  activity_data = scaling * pmax(t(activity_data), 0) %>% 
    tibble::as_tibble() %>% 
    dplyr::select(-.data$V1) # We don't need the first point since 0 for everyone
  
  last_col = paste0("V", n_points + 1)
  
  # Construct final dataset from group data and activity
  final_data = dplyr::bind_cols(data, activity_data) %>% 
    tidyr::pivot_longer(
      cols = .data$V2:.data[[last_col]],
      names_to = "time",
      values_to = "Xt"
    ) %>% 
    dplyr::select(-.data$time) %>% 
    tidyr::nest(data = c(.data$Xt)) %>% 
    dplyr::select(.data$subject_id, .data$group, Xt = data)
  
  
  return(final_data)
  
}
