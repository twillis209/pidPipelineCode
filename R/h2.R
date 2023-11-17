##' Transform heritability (h2) estimates from the observed scale to the liability scale
##'
##' @param h2.obs estimated h2 on the observed scale
##' @param P case proportion in study
##' @param K disease prevalence in population
##' ##' @importFrom stats dnorm qnorm
##' @export
obs_scale_to_liab_scale <- function(h2.obs, P, K) {
  z_2 <- dnorm(qnorm(1-K))^2

  h2.obs * ((K*(1-K))^2/(P*(1-P)))/z_2
}
