
#sigma <- matrix(c(1, 0, 0, 1), ncol = 2)
#sigma <- matrix(c(1, 0.8, 0.8, 1), ncol = 2)
#problem <- "multiple"
#problem <- "co-primary"
# target <- 0.3

vamos <- function(null, alternative, sigma, alpha_nom, beta_nom, problem, max_n, ...) {
  powered <- FALSE
  min_n <- 0
  final_n <- NULL
  while(max_n - min_n > 1) {
    n <- floor((min_n + max_n)/2)
    ocs <- vamos_n(n, null, alternative, sigma, alpha_nom, problem, ...)
    beta <- ocs[2]
    powered <- beta < beta_nom
    if(powered) {
      final_n <- n
      max_n <- n
    } else {
      min_n <- n
    }
  }
  n

  return(c(n, ocs))
}

vamos_n <- function(n, null, alternative, sigma, alpha_nom, problem, ...) {

  # For a given n, find the critical value which gives a type I error
  # rate of alpha_nom for the given problem and testing procedure

  # Get correlation of z statistics
  cor_z <- sigma[1,2]/(sqrt(sigma[1,1]*sigma[2,2]))
  variance <- sigma[1,1]

  if(problem == "co-primary"){

    ocs <- ocs_co_primary(n, target, cor_z, alpha_nom)


  } else if (problem == "multiple") {

    ocs <- ocs_mult_primary(n, target, cor_z, alpha_nom)

  } else {
    # value-based implementation
    vp <- list(...)
    ocs <- ocs_value_cp(n, null, alternative, variance, cor_z, alpha_nom,
                        vp$.a, vp$.c_x, vp$.c_y, vp$.n_y)
  }

  return(ocs)
}

test <- function(x, ...) {
  vp <- list(...)
  print(vp$b)
  for(i in vp) {
    print(x * i)
  }
}
