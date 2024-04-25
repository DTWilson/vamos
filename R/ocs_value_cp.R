
ocs_value_cp <- function(n, null, alternative, variance, cor_z, alpha_nom, a, c_x, c_y, n_y) {

  # For a given sample size, find the critical value which controls type I
  # error and find the corresponding worst case type II error

  # Translate hypothesis parameters to the z scale
  b_null <- null/sqrt(2*variance/n)
  b_alt <- alternative/sqrt(2*variance/n)

  # Find the optimal critical value
  r <- opt_crit_value_cp(cor_z, alpha_nom, a, c_x, c_y, n_y, b_null, b_alt)
  crit <- r[1]
  z_null <- r[2:3]

  # Get the corresponding tI error rate, reversing the transformation of the objective function
  tI <- sqrt(value_cp_alpha(crit, z_null, cor_z, alpha_nom, a, c_x, c_y, n_y, b_null)) + alpha_nom

  df <- search_points(a, c_x, c_y, n_y, b_alt)

  # Run a simple exhaustive search
  tII <- NULL
  if(a >= 0) {
    for(i in 1:nrow(df)){

      tII <- c(tII, 1 - integrate(f = int_cond_z_1_cp, lower = crit, upper = Inf,
                            mean = as.numeric(df[i,]),
                            cor_z=cor_z, a=a, c_x=c_x, c_y=c_y, n_y=n_y, crit=crit)$value)
    }
  } else {
    for(i in 1:nrow(df)){

      tII <- c(tII, 1 - integrate(f = int_cond_z_1_cp, lower = -Inf, upper = Inf,
                                  mean = as.numeric(df[i,]),
                                  cor_z=cor_z, a=a, c_x=c_x, c_y=c_y, n_y=n_y, crit=crit)$value)
    }
  }
  z_alt <- as.numeric(df[which.max(tII),])

  return(c(tI, max(tII), crit, z_null, z_alt))
}


opt_crit_value_cp <- function(cor_z, alpha_nom, a, c_x, c_y, n_y, b_null, b_alt) {

  # A starting critical value mid-way between the hypotheses
  crit <- (b_null + b_alt)/2
  error <- 1
  while(error > 10^-5) {
    # For the current critical value, find the point where tI error is maximised
    z <- z_1_to_maximise(crit, cor_z, alpha_nom, a, c_x, c_y, n_y, b_null)[[1]]
    #
    new_crit <- crit_to_minimise(z, cor_z, alpha_nom, a, c_x, c_y, n_y, b_null)
    error <- (crit - new_crit)^2
    crit <- new_crit
  }

  return(c(crit, z))
}

z_1_to_maximise <- function(crit, cor_z, alpha_nom, a, c_x, c_y, n_y, b_null, size = 100) {

  # A set of points in the parameter space on the null boundary
  df <- search_points(a, c_x, c_y, n_y, b_null, size)

  # At each point, calculate the tI error by integrating over the marginal
  # distribution in one dimension, calculating the conditional probs wrt to other
  tI <- NULL
  if(a >= 0){
    for(i in 1:nrow(df)){
      tI <- c(tI, integrate(f = int_cond_z_1_cp, lower = crit, upper = Inf,
                            mean = as.numeric(df[i,]),
                            cor_z=cor_z, a=a, c_x=c_x, c_y=c_y, n_y=n_y, crit=crit)$value)
    }
  } else {
    for(i in 1:nrow(df)){
      tI <- c(tI, integrate(f = int_cond_z_1_cp, lower = -Inf, upper = Inf,
                            mean = as.numeric(df[i,]),
                            cor_z=cor_z, a=a, c_x=c_x, c_y=c_y, n_y=n_y, crit=crit)$value)
    }
  }

  # Find the point where tI error is maximised
  z <- as.numeric(df[which.max(tI),])
  return(list(z, tIs = cbind(df, tI)))
}

crit_to_minimise <- function(z, cor_z, alpha_nom, a, c_x, c_y, n_y, b_null) {

  # Optimise the critical value so we get the nominal tI error at a specific
  # point in the parameter space z
  crit <- optimise(value_cp_alpha, lower = b_null, upper = 15, mean = z,
                   cor_z = cor_z, alpha_nom = alpha_nom,
                   a = a, c_x = c_x, c_y = c_y, n_y = n_y, b_null = b_null)$minimum

  return(crit)

}

value_cp_alpha <- function(crit, mean, cor_z, alpha_nom, a, c_x, c_y, n_y, b_null) {

  if(a >= 0) {
    alpha <- integrate(f = int_cond_z_1_cp, lower = crit, upper = Inf,
                       mean = mean,
                       cor_z=cor_z, a=a, c_x=c_x, c_y=c_y, n_y=n_y, crit=crit)$value
  } else {
    alpha <- integrate(f = int_cond_z_1_cp, lower = -Inf, upper = Inf,
                       mean = mean,
                       cor_z=cor_z, a=a, c_x=c_x, c_y=c_y, n_y=n_y, crit=crit)$value
  }

  return((alpha - alpha_nom)^2)
}

int_cond_z_1_cp <- function(z_1, mean, cor_z, a, c_x, c_y, n_y, crit) {
  # For a given z_1, calculate the conditional probability of landing in the
  # critical region based on the conditional distribution z_2 | z_1, when
  # the joint dist of (z_1, z_2) has a given mean and correlation

  # Get mean m and sd s of the conditional distribution
  m <- mean[2] + cor_z*(z_1 - mean[1])
  sd <- sqrt(1 - cor_z^2)

  w <- (c_x)/(c_y - n_y - a*c_x*n_y)
  s <- (c_x - crit)/c_x

  crit_y <- c_y - s*(c_y - (w*c_y - c_x)/(w + w*a*c_x))

  # Find the lower limit for the integration over z_2
  if(a >= 0) {
    low <- ifelse(z_1 > c_x, crit_y,
                  ifelse(z_1 > crit, c_y - s*(c_y - (w*c_y - (c_x - (c_x - z_1)/s))/(w + w*a*(c_x - (c_x - z_1)/s))),
                         Inf))
  } else {
    low <- ifelse(z_1 < c_x, crit_y,
                  ifelse(z_1 < crit, c_y - s*(c_y - (w*c_y - (c_x - (c_x - z_1)/s))/(w + w*a*(c_x - (c_x - z_1)/s))),
                         -Inf))
  }

  # Calculate the conditional prob and weight by the marginal
  # prob of z_1
  (1 - pnorm(low, m, sd))*dnorm(z_1, mean[1], 1)
}





check_errors <- function(v) {
  # To check our error rates are maximised properly, construct plots of:
  # - the null boundary with the maximising point
  # - the type I error over the null boundary; and
  # - the alternative boundary with the maximising point
  # - the type II error over the null boundary.
  # Note, we use a large number of search points here.

  df <- z_1_to_maximise(v$crit, v$sigma[1,2], v$alpha_nom, v$a, v$c_x, v$c_y, v$n_y, v$b_null, size = 1000)[[2]]
  df$i <- 1:nrow(df)

  p_tI_1 <- ggplot(df, aes(z_1, z_2)) + geom_line(aes(colour = tI)) +
    geom_point(data = df[which.max(df[,3]),]) +
    theme_minimal()

  p_tI_2 <- ggplot(df, aes(i, tI)) + geom_line(aes(colour = tI)) +
    geom_point(data = df[which.max(df[,3]),]) +
    theme_minimal()

  df_alt <- z_1_to_maximise(v$crit, v$sigma[1,2], v$alpha_nom, v$a, v$c_x, v$c_y, v$n_y, v$b_alt, size = 1000)[[2]]
  df_alt$tII <- 1 - df_alt$tI
  df_alt$i <- 1:nrow(df_alt)

  p_tII_1 <- ggplot(df_alt, aes(z_1, z_2)) + geom_line(aes(colour = tII)) +
    geom_point(data = df_alt[which.max(df_alt[,4]),]) +
    theme_minimal()

  p_tII_2 <- ggplot(df_alt, aes(i, tII)) + geom_line(aes(colour = tII)) +
    geom_point(data = df_alt[which.max(df_alt[,4]),]) +
    theme_minimal()

  return(list(p_tI_1, p_tI_2, p_tII_1, p_tII_2))
}


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
