

new_vamos <- function(n, null, alternative, sigma, alpha_nom, beta_nom,
                      a, c_x, c_y, n_y, sol, b_null, b_alt) {
  structure(list(n = n, null = null, alternative = alternative, sigma = sigma,
                 alpha_nom = alpha_nom, beta_nom = beta_nom,
                 a = a, c_x = c_x, c_y = c_y, n_y = n_y,
                 crit = sol[1], z_null = sol[2:3], z_alt = sol[4:5],
                 b_null = b_null, b_alt = b_alt),
            class = "vamos")
}

#' @export
vamos <- function(null = 0, alternative, sigma = matrix(c(1, 0, 0, 1), nrow = 2),
                  alpha_nom = 0.05, beta_nom = 0.2,
                  a = 1, c_x = 5, c_y = c_x, n_y = 0,
                  max_n = 1000) {

  r <- optimise_n(null, alternative, sigma, alpha_nom, beta_nom, max_n,
                  a, c_x, c_y, n_y)

  sol <- r[4:8]

  b_null <- null/sqrt(2*sigma[1,1]/r[1])
  b_alt <- alternative/sqrt(2*sigma[1,1]/r[1])

  return(new_vamos(n=r[1], null, alternative, sigma, alpha_nom, beta_nom,
                   a, c_x, c_y, n_y, sol, b_null, b_alt))

}


#' @export
#'
print.vamos <- function(x, ...) {
  x$n
}



#' @export
#'
plot.vamos <- function(x, y, ...) {
  # Grid of points for plotting
  df <- expand.grid(x = seq(-5,15,0.3),
                    y = seq(-5,15,0.3))

  # Get weight for value function
  a <- x$a; c_x <- x$c_x; c_y <- x$c_y; n_y <- x$n_y
  w <- (c_x)/(c_y - n_y - a*c_x*n_y)

  b_null <- x$null/sqrt(2*x$sigma[1,1]/x$n)
  b_alt <- x$alternative/sqrt(2*x$sigma[1,1]/x$n)
  crit <- x$crit

  # Get the required shift parameters to transform the basic hypothesis (set up
  # assuming n_x = 0) to null and alternatives
  s_n <- (c_x - b_null)/c_x
  s_a <- (c_x - b_alt)/c_x
  s_c <- (c_x - crit)/c_x

  # Null hypothesis
  x_n <- seq(b_null, c_x, length.out = 100)
  y_n <- c_y - s_n*(c_y - (w*c_y - (c_x - (c_x - x_n)/s_n))/(w + w*a*(c_x - (c_x - x_n)/s_n)))

  # Alternative hypothesis
  x_a <- seq(b_alt, c_x, length.out = 100)
  y_a <- c_y - s_a*(c_y - (w*c_y - (c_x - (c_x - x_a)/s_a))/(w + w*a*(c_x - (c_x - x_a)/s_a)))

  # Critical region
  x_c <- seq(crit, c_x, length.out = 100)
  y_c <- c_y - s_c*(c_y - (w*c_y - (c_x - (c_x - x_c)/s_c))/(w + w*a*(c_x - (c_x - x_c)/s_c)))

  df <- data.frame(x = c(x_n, x_a),
                   y = c(y_n, y_a),
                   h = rep(c("null", "alt"), each = 100))

  df_crit <- data.frame(x = x_c, y = y_c)

  p <- ggplot(df, aes(x, y)) + geom_line(aes(colour = h)) +
    xlim(c(-15, 15)) + ylim(c(-15, 15)) +
    scale_color_manual(name = "Hypothesis", breaks = c("null", "alt"),
                       values = c(2,4), labels = c("N", "A")) +
    geom_line(data = df_crit, linetype = 2) +
    xlab(expression(z[1])) + ylab(expression(z[2])) +
    coord_fixed() +
    theme_minimal()

  if(x$a >= 0) {
    p <- p + geom_segment(aes(x = b_null, xend = b_null, y =  15, yend = c_y, colour = "null")) +
      geom_segment(aes(x = c_x, xend = 15, y =  n_y, yend = n_y, colour = "null")) +
      geom_segment(aes(x = b_alt, xend = b_alt, y =  15, yend = c_y, colour = "alt")) +
      geom_segment(aes(x = c_x, xend = 15, y =  y_a[length(y_a)], yend = y_a[length(y_a)], colour = "alt")) +
      geom_segment(aes(x = crit, xend = crit, y =  15, yend = c_y), linetype = 2) +
      geom_segment(aes(x = c_x, xend = 15, y =  y_c[length(y_c)], yend = y_c[length(y_c)]), linetype = 2)
  } else {
    p <- p + geom_segment(aes(x = b_null, xend = b_null, y =  -15, yend = c_y, colour = "null")) +
      geom_segment(aes(x = c_x, xend = -15, y =  n_y, yend = n_y, colour = "null")) +
      geom_segment(aes(x = b_alt, xend = b_alt, y =  -15, yend = c_y, colour = "alt")) +
      geom_segment(aes(x = c_x, xend = -15, y =  y_a[length(y_a)], yend = y_a[length(y_a)], colour = "alt")) +
      geom_segment(aes(x = crit, xend = crit, y =  -15, yend = c_y), linetype = 2) +
      geom_segment(aes(x = c_x, xend = -15, y =  y_c[length(y_c)], yend = y_c[length(y_c)]), linetype = 2)
  }

  p <- p + geom_point(data = data.frame(x = x$z_null[1], y = x$z_null[2])) +
    geom_point(data = data.frame(x = x$z_alt[1], y = x$z_alt[2]))
  p
}
