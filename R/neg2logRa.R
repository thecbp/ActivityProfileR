#' Calculates -2logR(a) for a given activity index for data
#'
#' @param data dataset containing activity profiles
#' @param a the activity index that we want to calculate -2logR(a) at 
#' @param verbose TRUE/FALSE describing if you want to see progress of function
#'
#' @return
#' @export
#'
#' @importFrom magrittr %>%
#' @name %>%
#' 
#' @examples
neg2logRa = function(data, a, verbose = FALSE) {
  # Calculates the -2logR(a) test statistic for a given dataset and activity index
  
  # Progress bar since function is slow
  if (verbose) { print(paste0("Activity index: ", a))  }
  
  # Set up parameters for optimization loop
  # if no solution => sup{empty set} = 0 in numerator of -2logR
  n_grid_inc = 0 # times number of grid points for mu and lambda updated
  mu_grid_length = 20
  lambda_grid_length = 20
  gammas_js = unique(data$n) / sum(unique(data$n))
  n_groups = data %>% dplyr::pull(.data$group) %>% unique() %>% length()
  conditions_met = FALSE # F if error in computing -2logR(a), R otherwise
  
  # 6 was the arbirtrary limit set in the original code
  while (!conditions_met & n_grid_inc <= 6) {
    
    n_grid_inc = n_grid_inc + 1
    mu_grid_length = mu_grid_length * n_grid_inc + 1
    lambda_grid_length = lambda_grid_length * n_grid_inc + 1
    
    # Develop the grid of mu and lambda values for use in optimization
    # for the given activity index
    # Confirmed to work 04/06/20
    mu_lambda_grid = initialize_grid(data,
                                     mu_grid_length = mu_grid_length,
                                     lambda_grid_length = lambda_grid_length)
    
    # Find gammas and lambdas that satisfy estimating equations for -2logR(a)
    # Must use for loops since we potentially don't need to look at
    # all of the different mu/lambda combos
    for (i in 1:nrow(mu_lambda_grid)) {
      
      # Manipulate the current candidate mu and lambda to necessary format for multiroot
      current_params = mu_lambda_grid %>% 
        dplyr::slice(i) %>% 
        tidyr::pivot_longer(
          cols = .data$lambda1:.data$lambda2,
          names_to = "lambda_j",
          values_to = "lambdas"
        ) %>% 
        dplyr::mutate(
          lambda_js = cumsum(-.data$lambdas)
        )
      
      current_mu = current_params %>% dplyr::pull(.data$mu) %>% .data[[1]]
      current_lambda_js = current_params %>% dplyr::pull(.data$lambda_js)
      current_gammas = gammas_js - diff(c(0, current_lambda_js, 0)) * -current_mu
      
      # Use the multiroot function to solve for the estimating equation
      # Note to self: Why do we need to watch out for the error?
      # What does the error signify to us?
      
      # The original code has this for some reason?
      # My guess is that this is needed to check if the root actually computes but
      # you still want to store the result of the root elsewhere
      root = list() 
      root$estim.precis = 1  
      min_ps = rep(-1, times = n_groups)  
      
      root_attempt <- tryCatch(
        root <- rootSolve::multiroot(f = ee_model,
                          start = c(current_lambda_js, current_gammas),
                          data = data),
        error = function(e) e)
      
      # # Only move forward with the condition check if no error
      if (!inherits(root_attempt, "error")) {
        
        # Extract the lambdas and gammas from the roots
        root_lambdas = root$root[1:(n_groups - 1)]
        root_gammas = root$root[-(1:(n_groups - 1))]
        root_lambdas_kp1 = c(0, root_lambdas, 0)
        
        # Recalculate the p's based on the root lambdas and gamma
        group_ps = tibble::tibble(
          group = 1:n_groups
        ) %>% 
          dplyr::mutate(
            p = purrr::map(.data$group, function(g) {
              
              group_scaled_Ta = data %>% 
                dplyr::filter(.data$group == g) %>% dplyr::pull(.data$scaled_Ta)
              
              n = data %>% dplyr::pull(.data$n) %>% unique() %>% sum()
              
              # This is the given calculation for p in the original code
              1 / (n * (root_gammas[g] + (root_lambdas_kp1[g] - root_lambdas_kp1[(g + 1)]) * group_scaled_Ta))
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
        
        min_ps = group_ps %>%
          dplyr::group_by(.data$group) %>%
          dplyr::summarize(min_p = min(.data$p)) %>%
          dplyr::pull(.data$min_p)
        
        # Condition checks
        precision_met = root$estim.precis < 10^(-8)
        no_negative_p = sum(min_ps < 0) == 0
        is_not_error = !inherits(root_attempt, "error")
        
        # Condition check if we can break the loop
        # 1) The f.roots are close enough to zero
        # 2) None of the calculated minimum p's are negative
        # 3) There was no error produced by tryCatch
        # If these are met, evaluates to T, breaks for loop
        if (precision_met & no_negative_p & is_not_error) { break }
        
        # Current lambda and gamma did not satisfy condition, try the next combination
      }
    } # end FOR
    
    # Break the while loop with the fulfilled conditions if so
    conditions_met = precision_met & no_negative_p & is_not_error
    
  } # end WHILE
  
  # Return the optimized lambdas and gammas
  # With the optimized lambdas and gammas, calculate the -2logR(a) test statistic
  if (conditions_met) {
    sum_log_ps = group_ps %>% 
      dplyr::mutate( log_p = log(.data$p) ) %>% 
      dplyr::pull(.data$log_p) %>% sum()
    
    subjects_by_group = data %>% dplyr::pull(.data$n) %>% unique()
    
    ts = -2 * sum_log_ps - 2 * sum(subjects_by_group * log(subjects_by_group))
    return(ts)
  } else {
    return("Error occured.")
  }
  
}