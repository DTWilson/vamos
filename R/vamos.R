
vamos <- function(problem, alpha_nom, beta_nom, sigma, test = NULL) {

  sigma <- matrix(c(1, 0, 0, 1), ncol = 2)
  sigma <- matrix(c(1, 0.8, 0.8, 1), ncol = 2)
  problem <- "multiple"
  problem <- "co-primary"
  powered <- FALSE
  n <- 0
  while(!powered) {
    n <- n + 1
    beta <- vamos_n(n, problem, alpha_nom, sigma, test)
    powered <- beta < beta_nom
  }
  n

  return(n)
}

vamos_n <- function(n, problem, alpha_nom, sigma, test = NULL) {

  # For a given n, find the design (i.e. critical value) which gives a type I
  # error rate of alpha for the given problem and testing procedure

  # Get covariance matrix of z statistics
  cov_z <- sigma[1,2]/(2*sqrt(sigma[1,1]*sigma[2,2]))
  sigma_z <- matrix(c(1, cov_z, cov_z, 1),
                    ncol = 2)

  if(problem == "co-primary"){
    # Type I error maximised at the margins
    crit <- qnorm(1 - alpha_nom)
    # Type II error maximised at the point alternative
    beta <- 1 - pmvnorm(mean = c(0.3, 0.3)/sqrt(2/n),
                    lower = c(crit, crit), upper = c(Inf, Inf),
                    sigma = sigma_z)[1]

  } else if (problem == "multiple") {
    # Type I error maximised at the point null
    crit <- optim(0, multiple_f, lower = -5, upper = 5, method = "Brent",
                  sigma_z = sigma_z, alpha_nom = alpha_nom)$par
    # Type II maximised at the margins
    beta <- pnorm(crit, mean = 0.3/sqrt(2/n))

  } else {
    # value-based implementation
    beta <- vamos_n_value()
  }

  return(beta)
}

multiple_f <- function(crit, sigma_z, alpha_nom) {
  # Objective function to minimise when determining the critical value to use
  # in the multiple primary endpoint case
  alpha <- 1 - pmvnorm(mean = c(0, 0),
                   lower = c(-Inf, -Inf), upper =  c(crit, crit),
                   sigma = sigma_z)[1]
  return((alpha - alpha_nom)^2)
}

vamos_n_value <- function() {

}
