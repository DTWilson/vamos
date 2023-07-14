ocs_mult_primary <- function(n, target, cor_z, alpha_nom) {

  crit <- opt_crit_mult_primary(cor_z, alpha_nom)

  # Type II maximised at the margins
  beta <- pnorm(crit, mean = 0.3/sqrt(2/n))

  return(c(alpha_nom, beta))
}

opt_crit_mult_primary <- function(cor_z, alpha_nom) {

  # Type I error is maximised at the point null
  crit <- optim(0, mult_primary_obj, lower = -5, upper = 5, method = "Brent",
                cor_z = cor_z, alpha_nom = alpha_nom)$par

  return(crit)
}

mult_primary_obj <- function(crit, cor_z, alpha_nom) {

  # Objective function to minimise when determining the critical value to use
  # in the multiple primary endpoint case
  alpha <- 1 - pmvnorm(mean = c(0, 0),
                       lower = c(-Inf, -Inf), upper =  c(crit, crit),
                       sigma = matrix(c(1, cor_z, cor_z, 1), ncol = 2))[1]

  return((alpha - alpha_nom)^2)
}
