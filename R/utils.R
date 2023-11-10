
search_points <- function(a, c_x, c_y, b_y, b, size = 100) {

  # Build a data frame of points along a hypothesis boundary defined by b, the
  # marginal alternative for the x outcome
  w <- (c_x)/(c_y - b_y - a*c_x*b_y)

  # Get the canonical hypothesis passing through (0, c_y)
  df <- data.frame(x = seq(0, c_x, length.out = size))
  df$y <- (w*c_y - df$x)/(w + a*w*df$x)

  s <- (c_x - b)/c_x
  df <- as.data.frame(t( s*(t(df) - c(c_x, c_y)) + c(c_x, c_y) ))

  # Two cases, co_primary and multiple endpoints
  if(a >= 0){



    z_1 <- (seq(sqrt(0.001), sqrt(10), length.out = size))^2 + b
    z_2 <- c_y - s*(c_y - (w*c_y - (c_x - (c_x - z_1)/s))/(w + w*a*(c_x - (c_x - z_1)/s)))

    crit_y <- c_y - s*(c_y - (w*c_y - c_x)/(w + w*a*c_x))

    df <- data.frame(z_1 = z_1, z_2 = z_2)
    df <- df[df$z_1 <= c_x & df$z_2 <= c_y,]

    df <- rbind(data.frame(z_1 = rep(b, size), z_2 = seq(15, max(c_y, crit_y), length.out = size)),
                df,
                data.frame(z_1 = seq(max(c_x, b), 15, length.out = size), z_2 = rep(crit_y, size)))
  } else {
    z_1 <- -(seq(sqrt(0.001), sqrt(10), length.out = size))^2 + b
    z_2 <- c_y - s*(c_y - (w*c_y - (c_x - (c_x - z_1)/s))/(w + w*a*(c_x - (c_x - z_1)/s)))

    crit_y <- c_y - s*(c_y - (w*c_y - c_x)/(w + w*a*c_x))

    df <- data.frame(z_1 = z_1, z_2 = z_2)
    df <- df[df$z_1 >= c_x & df$z_2 >= c_y,]

    df <- rbind(data.frame(z_1 = rep(b, size), z_2 = seq(-15, min(c_y, crit_y), length.out = size)),
                df,
                data.frame(z_1 = seq(min(c_x, b), -15, length.out = size), z_2 = rep(crit_y, size)))

  }

  return(df)
}
