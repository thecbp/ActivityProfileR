#' Use bootstrap to calculate approximation of EL statistic
#' 
#' Function uses bootstrap to calculate the uniform approximation
#' of the EL statistic U^2 described in Chang et. al. This function 
#' is used inside EL_test
#'
#' @param data the dataset containing the activity profiles
#' @param n_boot how many bootstrap datasets we want to use
#'
#' @return a tibble with the results of all of the bootstrapping
#' @export
#'
#' @examples
bootstrap_U_star = function(data, n_boot = 1000) {
  # Bootstraps the calculation of the empricial likelihood test statistic
  #   i.e the value of U^* ("U-hat-star") from the data
  
  # Arguments:
  # data: tibble containing the subjects and their activity profiles
  #         assumes that the activity profiles are housed in a column named "Ta"
  # n_boot: how many bootstrap mean profiles should we create?
  
  # Returns:
  # A single vector representing the mean activity profile of the given data
  
  # Reconverting the list-col into a tibble to easily get the mean
  # for each activity level
  
  subjects_by_group = data %>% dplyr::pull(n) %>% unique()
  n_groups = data %>% dplyr::pull(group) %>% unique() %>% length()
  
  # Need to implement parallel processing for this eventually
  
  # Calculate the mean and variance profiles for each group
  observed_profile_stats = data %>% 
    tidyr::unnest(c(Ta, ai)) %>% 
    dplyr::group_by(group, ai) %>% # maybe ai h
    dplyr::summarize(
      n = unique(n),
      gamma = unique(gamma),
      mu_at_ai = mean(Ta),
      S2_at_ai = var(Ta) * ((n - 1) / n),
      wjs = 1 / (S2_at_ai / gamma)
    )
  
  # Normalize weights across groups to sum to 1 for each activity index
  observed_sum_weights = observed_profile_stats %>% 
    dplyr::group_by(ai) %>% 
    dplyr::summarize( sum_wjs = sum(wjs) )
  
  observed_profile_stats = observed_profile_stats %>% 
    dplyr::left_join(., observed_sum_weights, by = "ai") %>% 
    dplyr::mutate(
      normalized_wjs = wjs / sum_wjs,
      sq_normalized_wjs = sqrt(normalized_wjs)
    )
  
  # Create an overall tibble to house the bootstrapped datasets
  bootstrap_data_by_subject = tibble::tibble(
    boot_idx = rep(1:n_boot, each = n_groups),
    group_idx = rep(1:n_groups, times = n_boot)
  ) %>% 
    dplyr::mutate(
      # Boostrap the data by group and then recombine
      group_data = purrr::map(group_idx, function(gi) {
        dplyr::sample_n(data %>% dplyr::filter(group == gi), size = subjects_by_group[gi], replace = TRUE)
      })
    )
  
  bootstrap_data = tibble::tibble( boot_idx = 1:n_boot ) %>% 
    dplyr::mutate(
      bs_data = purrr::map(boot_idx, function(bi) {
        # Recombine the individual bootstrapped group data into an overall bootstrap dataset
        filtered_bs = bootstrap_data_by_subject %>% 
          dplyr::filter(boot_idx == bi) %>% 
          dplyr::pull(group_data)
        Reduce(function(x, y) { dplyr::bind_rows(x, y) }, filtered_bs)
      }),
      U2_boot_full = purrr::map(bs_data, function(bd) {
        # Calculate the mean activity profile by group in bootstrap data
        bs_stat_star = bd %>% 
          tidyr::unnest(c(Ta, ai)) %>% 
          dplyr::group_by(group, ai) %>% 
          dplyr::summarize(
            mu_star_at_ai = mean(Ta)
          )
        
        # Calculate U* using the bootstrap mean and the observed mean profiles
        U_hat_star = bs_stat_star %>% 
          dplyr::left_join(., observed_profile_stats, by = c("group", "ai")) %>% 
          dplyr::mutate(
            U_hat_star = sqrt(n) * (mu_star_at_ai - mu_at_ai) / sqrt(S2_at_ai)
          )
        
        # Calculate the weighted average U* across the groups
        wavg_U_hat_star = U_hat_star %>% 
          dplyr::group_by(ai) %>% 
          dplyr::summarize(
            wavg_U_hat_star = sum(sq_normalized_wjs * U_hat_star)
          )
        
        dplyr::left_join(U_hat_star, wavg_U_hat_star, by = "ai")
        # Calculate U^2 (uniform approximation of the EL statistic)
        # left_join(U_hat_star, wavg_U_hat_star, by = "ai") %>% 
        #   group_by(ai) %>% 
        #   summarize(
        #     U2_boot = sum(normalized_wjs * (U_hat_star / sq_normalized_wjs - wavg_U_hat_star)^2)
        #   )
      }),
      U2_boot = purrr::map(bs_data, function(bd) {
        # Calculate the mean activity profile by group in bootstrap data
        bs_stat_star = bd %>% 
          tidyr::unnest(c(Ta, ai)) %>% 
          dplyr::group_by(group, ai) %>% 
          dplyr::summarize(
            mu_star_at_ai = mean(Ta)
          )
        
        # Calculate U* using the bootstrap mean and the observed mean profiles
        U_hat_star = bs_stat_star %>% 
          dplyr::left_join(., observed_profile_stats, by = c("group", "ai")) %>% 
          dplyr::mutate(
            U_hat_star = sqrt(n) * (mu_star_at_ai - mu_at_ai) / sqrt(S2_at_ai)
          )
        
        # Calculate the weighted average U* across the groups
        wavg_U_hat_star = U_hat_star %>% 
          dplyr::group_by(ai) %>% 
          dplyr::summarize(
            wavg_U_hat_star = sum(sq_normalized_wjs * U_hat_star)
          )
        
        dplyr::left_join(U_hat_star, wavg_U_hat_star, by = "ai") %>%
          dplyr::group_by(ai) %>%
          dplyr::summarize(
            U2_boot = sum(normalized_wjs * (U_hat_star / sq_normalized_wjs - wavg_U_hat_star)^2)
          )
      }),
      sup_boot = unlist(purrr::map(U2_boot, function(uboot) {
        uboot %>% 
          dplyr::summarize(
            sup_boot = max(U2_boot)
          )
      }))
    )
  
  return(bootstrap_data)
}