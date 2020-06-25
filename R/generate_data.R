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

  tibble::tibble(
    subject_id = 1:n_subj,
    group = rep(group_id, n_subj),
  ) %>% 
    dplyr::mutate(
      Xt = purrr::map(.data$group, function(i) { 
        
        # Unoptimal: simdiff doesn't allow just one simulation, 
        # so we make 2 and then just take the first one
        sim = ESGtoolkit::simdiff(n = 2, 
                                  horizon = n_points, 
                                  frequency = "annual", 
                                  model = "OU", 
                                  x0 = 0, 
                                  theta1 = thetas[1], 
                                  theta2 = thetas[2], 
                                  theta3 = thetas[3])
        
        # Scaling the data to match a desired range (aka to NHANES)
        # Removing the first item because the process always starts at zero
        tibble::tibble(
          Xt = scaling * c(pmax(t(sim[,1]), 0))[-1]
        )
      })
    )

}
