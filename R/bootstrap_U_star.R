#' Use bootstrap to calculate approximation of EL statistic
#' 
#' Function uses bootstrap to calculate the uniform approximation
#' of the EL statistic U^2 described in Chang et. al. This function 
#' is used inside EL_test
#'
#' @param data the dataset containing the activity profiles
#' @param group_col the column that contains the group ID column
#' @param id_col the column that contains the subject ID column
#' @param n_boot how many bootstrap datasets we want to use
#'
#' @return a tibble with the results of all of the bootstrapping
#' @export
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @name %>%
#'
bootstrap_U_star = function(data, 
                            group_col,
                            id_col,
                            n_boot = 1000) {
  
  # Quoting the variables so we don't need to quote the variables in code
  group_col = dplyr::enquo(group_col)
  id_col = dplyr::enquo(id_col)
  
  # Calculate some useful values for the bootstrapping
  subjects_by_group = data %>% dplyr::pull(.data$n) %>% unique()
  n_groups = data %>% dplyr::pull(!!group_col) %>% unique() %>% length()

  # Calculate the mean and variance profiles for each group and activity index
  observed_profile_stats = data %>% 
    dplyr::group_by(!!group_col, .data$ai) %>%
    dplyr::summarize(
      n = unique(.data$n),
      gamma = unique(.data$gamma),
      mu_at_ai = mean(.data$Ta),
      S2_at_ai = stats::var(.data$Ta) * ((.data$n - 1) / .data$n),
      wjs = 1 / (.data$S2_at_ai / .data$gamma)
    )
  
  # Normalize weights across groups to sum to 1 for each activity index
  observed_sum_weights = observed_profile_stats %>% 
    dplyr::group_by(.data$ai) %>% 
    dplyr::summarize( sum_wjs = sum(.data$wjs) )
  
  observed_profile_stats = observed_profile_stats %>% 
    dplyr::left_join(observed_sum_weights, by = "ai") %>% 
    dplyr::mutate(
      normalized_wjs = .data$wjs / .data$sum_wjs,
      sq_normalized_wjs = sqrt(.data$normalized_wjs)
    )
  
  # Create an overall tibble to house the bootstrapped datasets
  # Each bootstrap sample needs the same amount of people, taking the same 
  # amount from each group, without uncoupling a subject from their data
  bootstrap_data_by_subject = tibble::tibble(
    boot_idx = rep(1:n_boot, each = n_groups),
    group_idx = rep(1:n_groups, times = n_boot)
  ) %>% 
    dplyr::mutate(
      # Re-nest the data to make sure each subject's data stays with that
      # person. Then filter for the group, to randomly pick from just that group
      group_data = purrr::map(.data$group_idx, function(gi) {
        gdata = data %>% 
          dplyr::group_by(!!id_col, !!group_col) %>% 
          tidyr::nest(data = c(.data$Ta, .data$ai)) %>% 
          dplyr::ungroup() %>% 
          dplyr::filter(!!group_col == gi)
        
        dplyr::sample_n(
          gdata, 
          size = subjects_by_group[gi], 
          replace = TRUE)
      })
    )
  
  bootstrap_data = tibble::tibble( 
    boot_idx = 1:n_boot 
    ) %>% 
    dplyr::mutate(
      bs_data = purrr::map(.data$boot_idx, function(bi) {
        
        # Recombine the individual bootstrapped group data 
        # into an overall bootstrap dataset
        filtered_bs = bootstrap_data_by_subject %>% 
          dplyr::filter(.data$boot_idx == bi) %>% 
          dplyr::pull(.data$group_data)
        
        full_data = Reduce(function(x, y) { 
          dplyr::bind_rows(x, y) 
          }, filtered_bs)
        
        full_data %>% 
          tidyr::unnest(data)
        
      }),
      U2_boot = purrr::map(.data$bs_data, function(bd) {
        
        # Calculate the mean activity profile by group in bootstrap data
        bs_stat_star = bd %>% 
          dplyr::group_by(!!group_col, .data$ai) %>% 
          dplyr::summarize(
            mu_star_at_ai = mean(.data$Ta)
          )
        
        # Calculate U* using the bootstrap mean and the observed mean profiles
        U_hat_star = bs_stat_star %>% 
          dplyr::left_join(observed_profile_stats, 
                           by = c(rlang::as_name(group_col), "ai")) %>% 
          dplyr::mutate(
            U_hat_star = sqrt(.data$n) * (.data$mu_star_at_ai - .data$mu_at_ai) / sqrt(.data$S2_at_ai)
          )
        
        # Calculate the weighted average U* across the groups
        wavg_U_hat_star = U_hat_star %>% 
          dplyr::group_by(.data$ai) %>% 
          dplyr::summarize(
            wavg_U_hat_star = sum(.data$sq_normalized_wjs * .data$U_hat_star)
          )
        
        dplyr::left_join(U_hat_star, wavg_U_hat_star, by = "ai") %>%
          dplyr::group_by(.data$ai) %>%
          dplyr::summarize(
            U2_boot = sum(.data$normalized_wjs * (.data$U_hat_star / .data$sq_normalized_wjs - .data$wavg_U_hat_star)^2)
          )
      }),
      sup_boot = unlist(purrr::map(.data$U2_boot, function(uboot) {
        uboot %>% 
          dplyr::summarize(
            sup_boot = max(.data$U2_boot)
          )
      }))
    )
  
  return(bootstrap_data)
}