
#sigma <- matrix(c(1, 0, 0, 1), ncol = 2)
#sigma <- matrix(c(1, 0.8, 0.8, 1), ncol = 2)
#problem <- "multiple"
#problem <- "co-primary"
# target <- 0.3

vamos <- function(target, sigma, alpha_nom, beta_nom) {
  powered <- FALSE
  n <- 0
  while(!powered) {
    n <- n + 1
    ocs <- vamos_n(n, target, sigma, alpha_nom, sigma)
    beta <- ocs[2]
    powered <- beta < beta_nom
  }
  n

  return(c(n, ocs))
}

vamos_n <- function(n, target, sigma, alpha_nom) {

  # For a given n, find the critical value which gives a type I error
  # rate of alpha_nom for the given problem and testing procedure

  # Get correlation of z statistics
  cor_z <- sigma[1,2]/(sqrt(sigma[1,1]*sigma[2,2]))

  if(problem == "co-primary"){

    ocs <- ocs_co_primary(n, target, cor_z, alpha_nom)


  } else if (problem == "multiple") {

    ocs <- ocs_mult_primary(n, target, cor_z, alpha_nom)

  } else {

    # value-based implementation
    ocs <- ocs_value(n, target, cor_z, alpha_nom, a, c, variance = sigma[1,1])
  }

  return(beta)
}
