#' Creates grid of mu and lambda values to optimize estimation equation
#'
#' @param data dataset containing the activity profiles
#' @param mu_grid_length number describing how long the mu grid should be
#' @param lambda_grid_length number describing how long the lambda grid should be
#' @param mu_tolerance number describing how fine the different mus should be in the grid
#' @param lambda_tolerance number describing how fine the different mus should be in the grid
#'
#' @return a tibble with all of the valid mu-lambda combinations from the data
#' @export
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
initialize_grid = function(data,
                           mu_grid_length = 20,
                           lambda_grid_length = 20,
                           mu_tolerance = 20,
                           lambda_tolerance = 20) {
  # Creates valid mu-lambda combinations for -2logR(a) optimization 
  
  # Arguments:
  # data: tibble containing the subjects and their activity profiles
  #         assumes that the activity profiles are housed in a column named "Ta"
  #         assumes that the patients are grouped by a column named "group"
  #         assumes that the data has already been filtered for data at index ai
  # mu_grid_length: number describing how long you want the mu grid to be
  # lambda_grid_length: number describing how long you want the lambda grid to be
  # mu_tolerance: number describing how fine we want the mu grid to be
  # lambda_tolerance: number describing how fine we want the lambda grid to be
  
  # Get the scaled minimum and maximum in each group
  scaled_minmax_by_group = data %>% 
    dplyr::group_by(.data$group) %>% 
    dplyr::summarize(
      scaled_Ta_min = min(.data$scaled_Ta),
      scaled_Ta_max = max(.data$scaled_Ta)
    )
  
  # Record the extremes of the scaled activity profiles
  scaled_maxmin = scaled_minmax_by_group %>% 
    dplyr::pull(.data$scaled_Ta_min) %>% max
  scaled_minmax = scaled_minmax_by_group %>% 
    dplyr::pull(.data$scaled_Ta_max) %>% min
  
  # Tibble construction for mu and lambda initialization
  mu_res = (scaled_minmax - scaled_maxmin) / mu_tolerance
  mu_grid = seq(scaled_maxmin + mu_res,
                scaled_minmax - mu_res,
                length = mu_grid_length)
  
  n_groups = data %>% dplyr::pull(.data$group) %>% unique() %>% length()
  
  mu_lambda_grid = tibble::tibble(
    mu = mu_grid
  ) %>% 
    dplyr::mutate(
      gj_data = purrr::map(.data$mu, function(u) {
        data %>% 
          dplyr::mutate( gj = (.data$scaled_Ta - u) / .data$gamma )
      }),
      limits = purrr::map(.data$gj_data, function(gd) {
        
        # Get the lower limits based on the positive gj
        ll = gd %>% 
          dplyr::filter(.data$gj > 0) %>% 
          dplyr::group_by(.data$group) %>% 
          dplyr::summarize(
            lower_limit = max((1/.data$n - 1) / .data$gj)
          ) 
        
        # Get the upper limit based on the negative gj
        ul = gd %>% 
          dplyr::filter(.data$gj < 0) %>% 
          dplyr::group_by(.data$group) %>% 
          dplyr::summarize(
            upper_limit = min((1/.data$n - 1) / .data$gj)
          ) %>% 
          dplyr::select(.data$upper_limit)
        
        # Calculate a tolerance based on each group's lmits
        dplyr::bind_cols(ll, ul) %>% 
          dplyr::mutate(
            tol = (.data$upper_limit - .data$lower_limit) / lambda_tolerance
          )
      }),
      lambdas = purrr::map(.data$limits, function(limit_data) {
        grid_data = limit_data %>% 
          dplyr::mutate(
            # Create the lambda grid based on each groups limits and tolerance
            grid = purrr::pmap(list(.data$upper_limit, 
                                    .data$lower_limit, 
                                    .data$tol), function(ul, ll, tol) {
              seq(ll + tol, ul - tol, length = lambda_grid_length)
            })
          ) %>% 
          dplyr::select(.data$group, .data$grid) %>%
          tidyr::unnest(.data$grid) # unravel the lambda grid
      }),
      valid_lambdas = purrr::map2(.data$limits, 
                                  .data$lambdas, function(lims, lams) {
                                    
        # Set up vectors for combination and matching against 3rd group limits
        g1_lambdas = lams %>% 
          dplyr::filter(.data$group == 1) %>% dplyr::pull(.data$grid)
        g2_lambdas = lams %>% 
          dplyr::filter(.data$group == 2) %>% dplyr::pull(.data$grid)
        
        # Expand the limits to match up properly against each 
        g3_lambda_expand = lams %>% 
          dplyr::filter(.data$group == 3) %>% dplyr::pull(.data$grid)
        
        tidyr::expand_grid(lambda1 = g1_lambdas, lambda2 = g2_lambdas) %>%
          dplyr::arrange(.data$lambda2) %>% # This is how it's sorted in the original code
          dplyr::mutate(
            lambda_km = -(.data$lambda1 + .data$lambda2),
            g3_upper = lims %>% 
              dplyr::filter(.data$group == 3) %>% 
              dplyr::pull(.data$upper_limit),
            g3_lower = lims %>% 
              dplyr::filter(.data$group == 3) %>% 
              dplyr::pull(.data$lower_limit)
          ) %>% 
          # We filter based on the limits of the last group, I think?
          dplyr::filter(
            .data$lambda_km > .data$g3_lower,
            .data$lambda_km < .data$g3_upper) %>%
          dplyr::select(.data$lambda1, .data$lambda2) # Up until this point, the grid matches exactly to Delta_grid_expand
      })
    ) %>% 
    dplyr::select(.data$mu, .data$valid_lambdas) %>% 
    tidyr::unnest_wider(.data$valid_lambdas) %>% 
    tidyr::unnest(c(.data$lambda1, .data$lambda2))
  
  return(mu_lambda_grid)

}