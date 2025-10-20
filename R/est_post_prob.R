# Given a sample of points from a bivariate posterior distribution
# and a critical region, estimate the prosterior probability of being
# in that critical region.

s <- data.frame(mu_x = runif(10^5, -2, 8),
                mu_y = runif(10^5, -5, 11))

est_post_prob <- function(crit, post_samples, a, c_x, c_y, b_y)
{
  #a = 1000; c_x = 6; c_y = 9; b_y = -2

  mu_x <- post_samples[,1]; mu_y <- post_samples[,2]

  # Get weight for value function
  w <- c_x/(c_y - b_y - a*c_x*b_y)

  # Critical region
  s <- (c_x - crit)/c_x

  crit_y <- c_y - s*(c_y - (w*c_y - c_x)/(w + w*a*c_x))

  crit_xy <- ifelse(mu_x > c_x, crit_y,
         ifelse(mu_x > crit, c_y - s*(c_y - (w*c_y - (c_x - (c_x - mu_x)/s))/(w + w*a*(c_x - (c_x - mu_x)/s))),
                Inf))

  post_samples$go <- mu_y > crit_xy

  return(list(mean(post_samples$go), post_samples))
}

get_post_samples <- function(prior_point)
{
  # sample data
  s_x <- 1; s_y <-1; tau <- 0.5
  cov_zs <- tau/(s_x*s_y)
  Sig <- matrix(c(1, cov_zs, cov_zs, 1), ncol=2)
  d <- rmvnorm(1, mean = as.numeric(prior_point),
               sigma = Sig)

  # posterior hyper-parameters
  S <- matrix(c(0.01, 0.005, 0.005, 0.1), ncol=2)
  m <- c(0,0)
  post <- rmvnorm(10^4,
                  mean =
                  sigma = Sig)
}

a = 10; c_x = 6; c_y = 9; b_y = -2; b_null <- 0; b_alt <- 2
crit <- seq(0,6,0.1)

s_null <- est_post_prob(crit = b_null, post_samples = s, a, c_x, c_y, b_y)[[2]]
s_null <- s_null[!s_null$go, 1:2]

s_alt <- est_post_prob(crit = b_alt, post_samples = s, a, c_x, c_y, b_y)[[2]]
s_alt <- s_alt[s_alt$go, 1:2]
