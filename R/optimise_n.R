
optimise_n <- function(null, alternative, sigma, alpha_nom, beta_nom, max_n, a, c_x, c_y, b_y) {

  # Get correlation of z statistics
  cor_z <- sigma[1,2]/(sqrt(sigma[1,1]*sigma[2,2]))
  variance <- sigma[1,1]

  # A bisection search for the smallest n which satisfies error rate constraints
  powered <- FALSE
  min_n <- 0
  final_n <- NULL
  cat("Searching for optimal n")
  while(max_n - min_n > 1) {
    cat(".")
    n <- floor((min_n + max_n)/2)
    ocs <- ocs_value_cp(n, null, alternative, variance, cor_z, alpha_nom,
                        a, c_x, c_y, b_y)
    beta <- ocs[2]
    powered <- beta < beta_nom
    if(powered) {
      final_n <- n
      max_n <- n
    } else {
      min_n <- n
    }
  }
  cat("\n")
  n

  return(c(n, ocs))
}
