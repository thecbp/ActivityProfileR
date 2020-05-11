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
#' @examples
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
    dplyr::group_by(group) %>% 
    dplyr::summarize(
      scaled_Ta_min = min(scaled_Ta),
      scaled_Ta_max = max(scaled_Ta)
    )
  
  # Record the extremes of the scaled activity profiles
  scaled_maxmin = scaled_minmax_by_group %>% dplyr::pull(scaled_Ta_min) %>% max()
  scaled_minmax = scaled_minmax_by_group %>% dplyr::pull(scaled_Ta_max) %>% min()
  
  # Tibble construction for mu and lambda initialization
  mu_res = (scaled_minmax - scaled_maxmin) / mu_tolerance
  mu_grid = seq(scaled_maxmin + mu_res,
                scaled_minmax - mu_res,
                length = mu_grid_length)
  
  n_groups = data %>% dplyr::pull(group) %>% unique() %>% length()
  
  mu_lambda_grid = tibble::tibble(
    mu = mu_grid
  ) %>% 
    mutate(
      gj_data = purrr::map(mu, function(u) {
        data %>% dplyr::mutate( gj = (scaled_Ta - u) / gamma )
      }),
      limits = map(gj_data, function(gd) {
        
        # Get the lower limits based on the positive gj
        ll = gd %>% 
          dplyr::filter(gj > 0) %>% 
          dplyr::group_by(group) %>% 
          dplyr::summarize(
            lower_limit = max((1/n - 1) / gj)
          ) 
        
        # Get the upper limit based on the negative gj
        ul = gd %>% 
          dplyr::filter(gj < 0) %>% 
          dplyr::group_by(group) %>% 
          dplyr::summarize(
            upper_limit = min((1/n - 1) / gj)
          ) %>% 
          dplyr::select(upper_limit)
        
        # Calculate a tolerance based on each group's lmits
        dplyr::bind_cols(ll, ul) %>% 
          dplyr::mutate(
            tol = (upper_limit - lower_limit) / lambda_tolerance
          )
      }),
      lambdas = purrr::map(limits, function(limit_data) {
        grid_data = limit_data %>% 
          dplyr::mutate(
            # Create the lambda grid based on each groups limits and tolerance
            grid = purrr::pmap(list(upper_limit, lower_limit, tol), function(ul, ll, tol) {
              seq(ll + tol, ul - tol, length = lambda_grid_length)
            })
          ) %>% 
          dplyr::select(group, grid) %>%
          tidyr::unnest(grid) # unravel the lambda grid
      }),
      valid_lambdas = purrr::map2(limits, lambdas, function(lims, lams) {
        # Set up vectors for combination and matching against 3rd group limits
        g1_lambdas = lams %>% dplyr::filter(group == 1) %>% dplyr::pull(grid)
        g2_lambdas = lams %>% dplyr::filter(group == 2) %>% dplyr::pull(grid)
        
        # Expand the limits to match up properly against each 
        g3_lambda_expand = lams %>% dplyr::filter(group == 3) %>% dplyr::pull(grid)
        tidyr::expand_grid(lambda1 = g1_lambdas, lambda2 = g2_lambdas) %>%
          dplyr::arrange(lambda2) %>% # This is how it's sorted in the original code
          dplyr::mutate(
            lambda_km = -(lambda1 + lambda2),
            g3_upper = lims %>% dplyr::filter(group == 3) %>% dplyr::pull(upper_limit),
            g3_lower = lims %>% dplyr::filter(group == 3) %>% dplyr::pull(lower_limit)
          ) %>% 
          # We filter based on the limits of the last group, I think?
          dplyr::filter(
            lambda_km > g3_lower,
            lambda_km < g3_upper) %>%
          dplyr::select(lambda1, lambda2) # Up until this point, the grid matches exactly to Delta_grid_expand
      })
    ) %>% 
    dplyr::select(mu, valid_lambdas) %>% 
    tidyr::unnest_wider(valid_lambdas) %>% 
    tidyr::unnest(c(lambda1, lambda2))
  
  return(mu_lambda_grid)

}