
ocs_value_cp_wrong <- function(n, null, alternative, sigma, alpha_nom, a, c_x, n_y) {

  # Optimal OCs for value-based hypotheses of the co-primary nature, but when
  # using co-primary style critical region. This uses the same basic approach
  # as the value based case, but now it's generally simpler to estimate error
  # rates - we can use pmvnorm to do the integrations, rather than integrate()

  # Get correlation of z statistics
  cor_z <- sigma[1,2]/(sqrt(sigma[1,1]*sigma[2,2]))
  variance <- sigma[1,1]
  c_y <- c_x

  b_null <- null/sqrt(2*variance/n)
  b_alt <- alternative/sqrt(2*variance/n)

  opt <- optimise(value_cp_wrong_obj, lower = b_null, upper = 5,
                  cor_z = cor_z, alpha_nom = alpha_nom,
                  a = a, c_x = c_x, c_y = c_y, n_y = n_y, b_null = b_null)

  crit <- opt$minimum
  tI <- sqrt(opt$objective) + alpha_nom

  df <- search_points(a, c_x, c_y, n_y, b_alt)

  # Run a a simple exhaustive search
  tII <- NULL
  if(a >= 0) {
    for(i in 1:nrow(df)){
      tII <- c(tII, 1 - pmvnorm(lower = c(crit, crit), upper = c(Inf, Inf),
                                mean = as.numeric(df[i,]), sigma = matrix(c(1, cor_z, cor_z, 1), ncol = 2)))
    }
  } else {
    for(i in 1:nrow(df)){
      tII <- c(tII, pmvnorm(lower = c(-Inf, -Inf), upper = c(crit, crit),
                            mean = as.numeric(df[i,]), sigma = matrix(c(1, cor_z, cor_z, 1), ncol = 2)))
    }
  }
  max(tII)

  return(c(tI, max(tII)))
}

value_cp_wrong_obj <- function(crit, cor_z, alpha_nom, a, c_x, c_y, n_y, b_null) {
  # Objective junction to use when searching for a critical value that will give
  # us the nominal type I error rate.

  # Search over the null hypothesis boundary to find the maximum type I error
  # rate.
  df <- search_points(a, c_x, c_y, n_y, b_null)

  # Run a a simple exhaustive search
  tI <- NULL
  if(a >= 0) {
    for(i in 1:nrow(df)){
      tI <- c(tI, pmvnorm(lower = c(crit, crit), upper = c(Inf, Inf),
                          mean = as.numeric(df[i,]), sigma = matrix(c(1, cor_z, cor_z, 1), ncol = 2)))
    }
  } else {
    for(i in 1:nrow(df)){
      tI <- c(tI, 1 - pmvnorm(lower = c(-Inf, -Inf), upper = c(crit, crit),
                              mean = as.numeric(df[i,]), sigma = matrix(c(1, cor_z, cor_z, 1), ncol = 2)))
    }
  }
  max(tI)

  # Return the penalised objective to be minimised
  return((max(tI)  - alpha_nom)^2)
}


