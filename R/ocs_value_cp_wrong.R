
ocs_value_cp_wrong <- function(n, null, alternative, cor_z, alpha_nom, variance, a, c_x, c_y, n_y) {

  # Optimal OCs for value-based hypotheses of the co-primary nature, but when
  # using co-primary style critical region.

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

value_cp_obj <- function(crit, cor_z, alpha_nom, a, c_x, c_y, n_y, b_null) {

  # Search over the null hypothesis boundary to find the maximum type I error
  # rate.
  df <- search_points(a, c_x, c_y, n_y, b_null)

  # Run a a simple exhaustive search
  tI <- NULL
  for(i in 1:nrow(df)){

    tI <- c(tI, integrate(f = int_cond_z_1_cp, lower = crit, upper = Inf,
                          mean = as.numeric(df[i,]),
                          cor_z=cor_z, a=a, c_x=c_x, c_y=c_y, n_y=n_y, crit=crit)$value)
  }
  which.max(tI)

  # Return the penalised objective to be minimised
  return((max(tI)  - alpha_nom)^2)
}
