# Given a sample of points from a bivariate posterior distribution
# and a critical region, estimate the prosterior probability of being
# in that critical region.

# est_post_prob <- function(crit, post_samples, a, c_x, c_y, b_y)
# {
#   #a = 1000; c_x = 6; c_y = 9; b_y = -2
#
#   mu_x <- post_samples[,1]; mu_y <- post_samples[,2]
#
#   # Get weight for value function
#   w <- c_x/(c_y - b_y - a*c_x*b_y)
#
#   # Critical region
#   s <- (c_x - crit)/c_x
#
#   crit_y <- c_y - s*(c_y - (w*c_y - c_x)/(w + w*a*c_x))
#
#   crit_xy <- ifelse(mu_x > c_x, crit_y,
#          ifelse(mu_x > crit, c_y - s*(c_y - (w*c_y - (c_x - (c_x - mu_x)/s))/(w + w*a*(c_x - (c_x - mu_x)/s))),
#                 Inf))
#
#   post_samples$go <- mu_y > crit_xy
#
#   return(list(mean(post_samples$go), post_samples))
# }
#
# get_post_samples <- function(prior_point, m0, S0)
# {
#   # sample data
#   s_x <- 1; s_y <-1; tau <- 0.005
#   n <- 150
#   #cov_zs <- tau/(s_x*s_y)
#   #Sig <- matrix(c(1, cov_zs, cov_zs, 1), ncol=2)
#   Sig <- matrix(c(2*s_x^2/n, 2*tau/n, 2*tau/n, 2*s_y^2/n), ncol=2)
#   d <- t(rmvnorm(1, mean = as.numeric(prior_point),
#                sigma = Sig))
#
#   # posterior samples
#   post <- rmvnorm(10^5,
#                   mean = solve(solve(S0) + solve(Sig)) %*% (solve(S0) %*% m0 + solve(Sig) %*% d),
#                   sigma = solve(solve(S0) + solve(Sig)))
#   post
# }
#
# m0 <- c(0.3, 0.3)
# S0 <- matrix(c(0.3^2, 0.02, 0.02, 0.3^2), ncol = 2)
#
# s <- rmvnorm(10^5, mean = m0, sigma = S0)
# s <- as.data.frame(s)
#
# a = 1; c_x = 0.6; c_y = 0.6; b_y = 0
# b_null <- 0; b_alt <- 0.3
#
# nulls <- est_post_prob(crit = b_null, post_samples = s, a, c_x, c_y, b_y)[[2]]
# nulls <- !nulls$go
#
# alts <- est_post_prob(crit = b_alt, post_samples = s, a, c_x, c_y, b_y)[[2]]
# alts <- alts$go
#
# s$h <- ifelse(nulls, "null",
#               ifelse(alts, "alt", "mid"))
#
# crits <- seq(0.05, 0.25, 0.05)
#
# ptm <- proc.time()
# gos <- matrix(rep(0, length(crits)*2), ncol = 2)
# for(i in 1:1000){#nrow(s_null)){
#   prior_point <- s_null[i,]
#   post <- get_post_samples(prior_point, m0, S0)
#   go_b <- go_f <- NULL
#   for(j in 1:length(crits)){
#     p_b <- est_post_prob(crits[j], as.data.frame(post), a=a, c_x, c_y, b_y)[[1]]
#     p_f <- est_post_prob(crits[j], as.data.frame(post), a=10^5, c_x, c_y, b_y)[[1]]
#     gos[j,] <- gos[j,] + c(p_b > 0.5, p_f > 0.5)
#   }
# }
# gos <- gos/1000
#
# gos_a <- matrix(rep(0, length(crits)*2), ncol = 2)
# for(i in 1:1000){#nrow(s_alt)){
#   prior_point <- s_alt[i,]
#   post <- get_post_samples(prior_point, m0, S0)
#   go_b <- go_f <- NULL
#   for(j in 1:length(crits)){
#     p_b <- est_post_prob(crits[j], as.data.frame(post), a=a, c_x, c_y, b_y)[[1]]
#     p_f <- est_post_prob(crits[j], as.data.frame(post), a=10^5, c_x, c_y, b_y)[[1]]
#     gos_a[j,] <- gos_a[j,] + c(p_b > 0.5, p_f > 0.5)
#   }
# }
# gos_a <- gos_a/1000
# proc.time() - ptm
#
# df <- data.frame(tI = c(gos[,1], gos[,2]),
#                  tII = c(gos_a[,1], gos_a[,2]),
#                  t = rep(c("b", "f"), each = length(crits)))
#
# ggplot(df, aes(tI, tII, colour = t)) + geom_point()
#
#
#
#
# #########################
#
#
#
# m0 <- c(0.3, 0.3)
# S0 <- matrix(c(0.3^2, 0.02, 0.02, 0.3^2), ncol = 2)
#
# s <- rmvnorm(10^5, mean = m0, sigma = S0)
# s <- as.data.frame(s)
#
# a = 3; c_x = 0.6; c_y = 0.6; b_y = 0
# b_null <- 0; b_alt <- 0.3
#
# nulls <- est_post_prob(crit = b_null, post_samples = s, a, c_x, c_y, b_y)[[2]]
# nulls <- !nulls$go
#
# alts <- est_post_prob(crit = b_alt, post_samples = s, a, c_x, c_y, b_y)[[2]]
# alts <- alts$go
#
# s$h <- ifelse(nulls, "null",
#               ifelse(alts, "alt", "mid"))
#
# crits <- seq(0.01, 0.29, 0.025)
# ptm <- proc.time()
# N <- 10^4
# r_b <- matrix(rep(0, N*length(crits)), nrow = N)
# r_f <- matrix(rep(0, N*length(crits)), nrow = N)
# for(i in 1:N){
#   prior_point <- s[i,1:2]
#   post <- get_post_samples(prior_point, m0, S0)
#   go_b <- go_f <- NULL
#   for(j in 1:length(crits)){
#     p_b <- est_post_prob(crits[j], as.data.frame(post), a=a, c_x, c_y, b_y)[[1]]
#     p_f <- est_post_prob(crits[j], as.data.frame(post), a=10^5, c_x, c_y, b_y)[[1]]
#     r_b[i, j] <- p_b > 0.5
#     r_f[i, j] <- p_f > 0.5
#   }
# }
#
# rr <- NULL
# for(j in 1:length(crits)){
#   df_b <- data.frame(x = s[1:N, 1],
#                      y = s[1:N, 2],
#                      g = r_b[,j])
#
#   fit <- gam(g ~ ti(x) + ti(y) + ti(x,y), family = binomial(), data = df_b)
#
#   #ss <- s[1:N,]
#   ss <- data.frame(x = s[,1], y = s[,2], h = s[,3])
#   ss$p <- predict(fit, newdata = ss, type = "response")
#
#   ss_n <- ss[ss$h == "null",]
#   #ss_n[which.max(ss_n$p),]
#   tI_b <- max(ss_n$p)
#
#   ss_a <- ss[ss$h == "alt",]
#   #ss_a[which.min(ss_a$p),]
#   tII_b <- 1 - min(ss_a$p)
#
#   df_f <- data.frame(x = s[1:N, 1],
#                      y = s[1:N, 2],
#                      g = r_f[,j])
#
#   fit <- gam(g ~ ti(x) + ti(y) + ti(x,y), family = binomial(), data = df_f)
#
#   #ss <- s[1:N,]
#   ss <- data.frame(x = s[,1], y = s[,2], h = s[,3])
#   ss$p <- predict(fit, newdata = ss, type = "response")
#
#   ss_n <- ss[ss$h == "null",]
#   #ss_n[which.max(ss_n$p),]
#   tI_f <- max(ss_n$p)
#
#   ss_a <- ss[ss$h == "alt",]
#   #ss_a[which.min(ss_a$p),]
#   tII_f <- 1 - min(ss_a$p)
#
#   rr <- rbind(rr, c(tI_b, tII_b, 1))
#   rr <- rbind(rr, c(tI_f, tII_f, 2))
# }
#
# df <- data.frame(x = rr[,1],
#                  y = rr[,2],
#                  t = factor(rr[,3]))
#
# ggplot(df, aes(x, y, colour = t)) + geom_point()
# proc.time() - ptm
#
#
# df <- data.frame(x = ss_n[,1],
#                  y = ss_n[,2],
#                  p = ss_n[,4])
#
# ggplot(df, aes(x, y, colour = p)) + geom_point()
#
# df <- data.frame(x = ss_a[,1],
#                  y = ss_a[,2],
#                  p = ss_a[,4])
#
# ggplot(df, aes(x, y, colour = p)) + geom_point()
