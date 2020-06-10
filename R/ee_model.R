#' Estimating equation used in optimization of numerator in likelihood ratio
#'
#' @param lambda_gamma vector of parameters in equation, lambda followed by gamma
#' @param data dataset containing the activity profiles
#'
#' @return a vector same length as lambda_gamma that optimizes the equation
#' @export
#'
#' @importFrom magrittr %>%
#' @name %>%
#'
#' @examples
ee_model = function(lambda_gamma, data) {
  # Computes estimating equations in -2logR(a) calculation
  
  # Arguments:
  # data: tibble that contains the data on activity profiles filtered at 
  #       activity index a
  # lg: vector containing the startinglambdas and gammas to be used in multiroot()
  
  # Returns:
  # a vector containing... something
  
  # Compute some values of interest to the problem and prepare data for calculation
  n_groups = data %>% dplyr::pull(.data$group) %>% unique() %>% length()
  n = data %>% dplyr::pull(.data$n) %>% unique() %>% sum()
  lambdas = lambda_gamma[1:(n_groups - 1)]
  gammas = lambda_gamma[-(1:(n_groups - 1))]
  lambdas_kp1 = c(0, lambdas, 0)
  
  # Recalculate the p's based on the given lambda and gamma
  group_ps = tibble::tibble(
    group = 1:n_groups
  ) %>% 
    dplyr::mutate(
      p = purrr::map(.data$group, function(g) {
        group_scaled_Ta = data %>% 
          dplyr::filter(.data$group == g) %>% 
          dplyr::pull(.data$scaled_Ta)
        n = data %>% dplyr::pull(.data$n) %>% unique() %>% sum()
        
        # This is the given calculation for p in the original code
        1 / (n * (gammas[g] + (lambdas_kp1[g] - lambdas_kp1[(g + 1)]) * group_scaled_Ta))
      }),
      scaled_Ta = purrr::map(.data$group, function(g) {
        data %>% 
          dplyr::filter(.data$group == g) %>% 
          dplyr::pull(.data$scaled_Ta)
      }),
      Ta_p_prod = purrr::map2(.data$group, .data$p, function(g, ps) {
        group_scaled_Ta = data %>% 
          dplyr::filter(.data$group == g) %>% 
          dplyr::pull(.data$scaled_Ta)
        ps * group_scaled_Ta
      })
    ) %>% 
    tidyr::unnest(c(.data$p, .data$scaled_Ta, .data$Ta_p_prod))
  
  # Calculate the output of the estimating equations based on the p's and scaled Ta
  out = rep(0, 2 * n_groups - 1)
  for (g in 1:n_groups) {
    
    current_group_p = group_ps %>% 
      dplyr::filter(.data$group == g) %>% 
      dplyr::pull(.data$p)
    
    # Calculates the new gammas
    out[n_groups + g - 1] = sum(current_group_p) - 1
    
    if (g != n_groups) {
      
      current_group_prod = group_ps %>% 
        dplyr::filter(.data$group == g) %>% 
        dplyr::pull(.data$Ta_p_prod)
      
      next_group_prod = group_ps %>% 
        dplyr::filter(.data$group == (g + 1)) %>% 
        dplyr::pull(.data$Ta_p_prod)
      
      out[g] = sum(next_group_prod) - sum(current_group_prod)
    }
  }
  
  return(out)
}