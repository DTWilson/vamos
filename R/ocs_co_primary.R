ocs_co_primary <- function(n, target, cor_z, alpha_nom) {

  crit <- opt_crit_co_primary(alpha_nom)

  # Type II error is maximised at the point alternative, when both effect
  # sizes are equal to "target"
  beta <- 1 - pmvnorm(mean = c(target, target)/sqrt(2/n),
                      lower = c(crit, crit), upper = c(Inf, Inf),
                      sigma = matrix(c(1, cor_z, cor_z, 1), ncol = 2))[1]

  return(c(alpha_nom, beta))
}

opt_crit_co_primary <- function(alpha_nom) {

  # Type I error is maximised at the margins, so the critical value corresponds
  # to what we would use in a univariate test
  crit <- qnorm(1 - alpha_nom)

  return(crit)
}
