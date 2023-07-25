
ocs_value_cp <- function(n, null, alternative, cor_z, alpha_nom, variance, a, c_x, c_y, n_y) {

  # Optimal OCs for value-based approach of the co-primary nature

  b_null <- null/sqrt(2*variance/n)
  b_alt <- alternative/sqrt(2*variance/n)

  crit <- opt_crit_value_cp(cor_z, alpha_nom, a, c_x, c_y, n_y, b_null)

  tI <- sqrt(value_cp_obj(crit, cor_z, alpha_nom, a, c_x, c_y, n_y, b_null)) + alpha_nom

  df <- search_points(a, c_x, c_y, n_y, b_alt)

  # Run a a simple exhaustive search
  tII <- NULL
  for(i in 1:nrow(df)){
    tII <- c(tII, 1 - integrate(int_cond_z_1_cp, -Inf, Inf,
                          mean = as.numeric(df[i,]),
                          cor_z=cor_z, a=a, c_x=c_x, c_y=c_y, crit=crit)$value)
  }
  max(tII)

  return(c(tI, max(tII)))
}


opt_crit_value_cp <- function(cor_z, alpha_nom, a, c_x, c_y, n_y, b_null) {

  crit <- optim(0, value_cp_obj, lower = -5, upper = 5, method = "Brent",
                cor_z = cor_z, alpha_nom = alpha_nom,
                a = a, c_x = c_x, c_y = c_y, n_y = n_y, b_null = b_null)$par

  return(crit)
}


value_cp_obj <- function(crit, cor_z, alpha_nom, a, c_x, c_y, n_y, b_null) {

  # Search over the null hypothesis boundary to find the maximum type I error
  # rate.
  df <- search_points(a, c_x, c_y, n_y, b_null)

  # Run a a simple exhaustive search
  tI <- NULL
  for(i in 1:nrow(df)){
    tI <- c(tI, integrate(int_cond_z_1_cp, -Inf, Inf,
                          mean = as.numeric(df[i,]),
                          cor_z=cor_z, a=a, c_x=c_x, c_y=c_y, crit=crit)$value)
  }
  max(tI)

  # Return the penalised objective to be minimised
  return((max(tI)  - alpha_nom)^2)
}

int_cond_z_1_cp <- function(z_1, mean, cor_z, a, c_x, c_y, crit) {
  # For a given z_1, calculate the conditional probability of landing in the
  # critical region based on the conditional distribution z_2 | z_1, when
  # the joint dist of (z_1, z_2) has a given mean and correlation

  # Get mean m and sd s of the conditional distribution
  m <- mean[2] + cor_z*(z_1 - mean[1])
  s <- sqrt(1 - cor_z^2)

  w <- (c_x)/(c_y - n_y - a*c_x*n_y)
  s <- (c_x - crit)/c_x

  crit_y <- c_y - s*(c_y - (w*c_y - c_x)/(w + w*a*c_x))

  # Find the lower limit for the integration over z_2
  low <- ifelse(z_1 > c_x, crit_y,
                ifelse(z_1 > crit, c_y - s*(c_y - (w*c_y - (c_x - (c_x - z_1)/s))/(w + w*a*(c_x - (c_x - z_1)/s))),
                       Inf))

  # Calculate the conditional prob and weight by the marginal
  # prob of z_1
  (1 - pnorm(low, m, s))*dnorm(z_1, mean[1], 1)
}

ocs_value_cp_wrong <- function(n, null, alternative, cor_z, alpha_nom, variance, a, c_x, c_y, n_y) {

  # Optimal OCs for value-based hypotheses of the co-primary nature, but when
  # using co-primary style critical region.

  b_null <- null/sqrt(2*variance/n)
  b_alt <- alternative/sqrt(2*variance/n)

  opt <- optim(0, value_cp_wrong_obj, lower = -5, upper = 5, method = "Brent",
                cor_z = cor_z, alpha_nom = alpha_nom,
                a = a, c_x = c_x, c_y = c_y, n_y = n_y, b_null = b_null)

  crit <- opt$par
  tI <- sqrt(opt$value) + alpha_nom

  df <- search_points(a, c_x, c_y, n_y, b_alt)

  # Run a a simple exhaustive search
  tII <- NULL
  for(i in 1:nrow(df)){
    tII <- c(tII, 1 - pmvnorm(lower = c(crit, crit), upper = c(Inf, Inf),
                              mean = as.numeric(df[i,]), sigma = matrix(c(1, cor_z, cor_z, 1), ncol = 2)))
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
  for(i in 1:nrow(df)){
    tI <- c(tI, pmvnorm(lower = c(crit, crit), upper = c(Inf, Inf),
                        mean = as.numeric(df[i,]), sigma = matrix(c(1, cor_z, cor_z, 1), ncol = 2)))
  }
  max(tI)

  # Return the penalised objective to be minimised
  return((max(tI)  - alpha_nom)^2)
}


search_points <- function(a, c_x, c_y, n_y, b, size = 50) {
  # Build a data frame of points along a hypothesis boundary
  w <- (c_x)/(c_y - n_y - a*c_x*n_y)
  s <- (c_x - b)/c_x

  z_1 <- (seq(sqrt(0.001), sqrt(10), length.out = size))^2 + b
  z_2 <- c_y - s*(c_y - (w*c_y - (c_x - (c_x - z_1)/s))/(w + w*a*(c_x - (c_x - z_1)/s)))

  df <- data.frame(z_1 = z_1, z_2 = z_2)
  df <- df[df$z_1 <= c_x & df$z_2 <= c_y & df$z_2 >= n_y,]

  df <- rbind(data.frame(z_1 = rep(b, size), z_2 = seq(15, c_y, length.out = size)),
              df,
              data.frame(z_1 = seq(c_x, 15, length.out = size), z_2 = rep(df[nrow(df), "z_2"], size)))

  return(df)
}
