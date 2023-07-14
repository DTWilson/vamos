
ocs_value <- function(n, target, cor_z, alpha_nom, variance) {

  crit <- opt_crit_value(cor_z, alpha_nom, a, c)

  tI <- sqrt(value_obj(crit, cor_z, alpha_nom, a, c)) + alpha_nom

  # Find the distance from the null of 0 to the critical region
  # at the margins
  b_crit <- crit/(1 + a*(c - crit)) + crit

  # Find distance of alternative from null at the margins
  b <- target/sqrt(2*variance^2/n)
  # Get the shift parameter corresponding to the alternative hypothesis
  aa <- a; bb <- (-2 - a*c - b*a); cc <- b + b*a*c
  sh <- (-bb - sqrt(bb^2 - 4*aa*cc))/(2*aa)

  # Start search at the point where z_1 = z_2
  z_1_min <- (sqrt((2 - 2*a*sh)^2 - 4*a*(a*sh^2 - 2*sh - c)) - (2 - 2*a*sh))/(2*a)
  z_1 <- seq(z_1_min, 10, 0.1)
  # Find corresponding z_2 on the value-based boundary
  z_2 <- (c - (z_1 - sh))/(1 + a*(z_1 - sh)) + sh
  # Correct for when we are in the no-trade-off region
  z_2 <- pmax(z_2, b)

  # Run a a simple exhaustive search
  tII <- NULL
  for(i in 1:length(z_1)){
    tII <- c(tII, 1 - integrate(int_cond_z_1, -Inf, Inf,
                                mean = c(z_1[i], z_2[i]),
                                cor_z=cor_z, a=a, c=c, b_crit=b_crit, crit=crit)$value)
  }

  return(c(tI, max(tII)))
}


opt_crit_value <- function(cor_z, alpha_nom, a, c) {

  crit <- optim(0, value_obj, lower = -5, upper = 5, method = "Brent",
                cor_z = cor_z, alpha_nom = alpha_nom,
                a = a, c = c)$par

  return(crit)
}


value_obj <- function(crit, cor_z, alpha_nom, a, c) {

  b_crit <- crit/(1 + a*(c - crit)) + crit

  # Search over the null hypothesis boundary to find the maximum type I error
  # rate.
  # Start search at the point where z_1 = z_2
  z_1_min <- (sqrt(a*c + 1) - 1)/a
  z_1 <- seq(z_1_min, 10, 0.1)
  # Find corresponding z_2 on the value-based boundary
  z_2 <- (c - z_1)/(1 + a*z_1)
  # Correct for when we are in the no-trade-off region
  z_2 <- pmax(z_2, 0)

  # Run a a simple exhaustive search
  tI <- NULL
  for(i in 1:length(z_1)){
    tI <- c(tI, integrate(int_cond_z_1, -Inf, Inf,
                          mean = c(z_1[i], z_2[i]),
                          cor_z=cor_z, a=a, c=c, b_crit=b_crit, crit=crit)$value)
  }

  # Return the penalised objcetive to be minimised
  return((max(tI)  - alpha_nom)^2)
}


int_cond_z_1 <- function(z_1, mean, cor_z, a, c, b_crit, crit) {
  # For a given z_1, calculate the conditional probability of landing in the
  # critical region based on the conditional distribution z_2 | z_1, when
  # the joint dist of (z_1, z_2) has a given mean and correlation

  # Get mean m and sd s of the conditional distribution
  m <- mean[2] + cor_z*(z_1 - mean[1])
  s <- sqrt(1 - cor_z^2)

  # Find the lower limit for the integration over z_2
  low <- ifelse(z_1 > c, b_crit,
                ifelse(z_1 > b_crit, (c - (z_1- crit))/(1 + a*(z_1 - crit)) + crit,
                       Inf))

  # Calculate the conditional prob and weight by the marginal
  # prob of z_1
  (1 - pnorm(low, m, s))*dnorm(z_1, mean[1], 1)
}
