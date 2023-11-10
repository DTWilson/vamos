
new_vamos <- function(n, null, alternative, sigma, alpha_nom, beta_nom,
                      a, c_x, c_y, n_y, sol, b_null, b_alt) {
  # Constuctor

  structure(list(n = n, null = null, alternative = alternative, sigma = sigma,
                 alpha_nom = alpha_nom, beta_nom = beta_nom,
                 a = a, c_x = c_x, c_y = c_y, b_y = b_y,
                 crit = sol[1], z_null = sol[2:3], z_alt = sol[4:5],
                 b_null = b_null, b_alt = b_alt),
            class = "vamos")
}

#' @export
vamos <- function(null = 0, alternative, sigma = matrix(c(1, 0, 0, 1), nrow = 2),
                  alpha_nom = 0.05, beta_nom = 0.2,
                  a = 1, c_x = 5, c_y = c_x, n_y = 0,
                  max_n = 1000) {
  # Helper

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

  # Get weight for value function
  a <- x$a; c_x <- x$c_x; c_y <- x$c_y; b_y <- x$b_y
  w <- (c_x)/(c_y - b_y - a*c_x*b_y)

  b_null <- x$null/sqrt(2*x$sigma[1,1]/x$n)
  b_alt <- x$alternative/sqrt(2*x$sigma[1,1]/x$n)
  crit <- x$crit

  if(a >= 0){
    # co-primary
    x_up <- c_x + 2; y_up <- c_y + 2
    x_lo <- b_null - 2; y_lo <- b_y -2
  } else {
    # mutliple primary
    x_up <- b_alt + 2; y_up <- b_y + 2
    x_lo <- c_x - 2; y_lo <- c_y - 2
  }

  # Get weight for value function
  w <- (c_x)/(c_y - b_y - a*c_x*b_y)

  # Get the canonical hypothesis passing through (0, c_y)
  df <- data.frame(x = seq(0, c_x, length.out = 100))
  df$y <- (w*c_y - df$x)/(w + a*w*df$x)

  # Null
  s <- (c_x - b_null)/c_x
  df_null <- as.data.frame(t( s*(t(df) - c(c_x, c_y)) + c(c_x, c_y) ))

  # Alternative
  s <- (c_x - b_alt)/c_x
  df_alt <- as.data.frame(t( s*(t(df) - c(c_x, c_y)) + c(c_x, c_y) ))

  # Critical region
  s <- (c_x - crit)/c_x
  df_crit <- as.data.frame(t( s*(t(df) - c(c_x, c_y)) + c(c_x, c_y) ))

  df <- rbind(df_null, df_alt)
  df$h <- rep(c("null", "alt"), each = 100)


  ggplot(df, aes(x, y)) + geom_line(aes(colour = h)) +
    xlim(c(x_lo, x_up)) + ylim(c(y_lo, y_up)) +

    # null
    geom_segment(aes(x = b_null, xend = b_null, y = (a >= 0)*y_up + (a < 0)*y_lo, yend = c_y, colour = "null")) +
    geom_segment(aes(x = c_x, xend = (a >= 0)*x_up + (a < 0)*x_lo, y =  df_null$y[nrow(df_null)-1], yend = df_null$y[nrow(df_null)-1], colour = "null")) +

    # alternative
    geom_segment(aes(x = b_alt, xend = b_alt, y = (a >= 0)*y_up + (a < 0)*y_lo, yend = c_y, colour = "alt")) +
    geom_segment(aes(x = c_x, xend = (a >= 0)*x_up + (a < 0)*x_lo, y =  df_alt$y[nrow(df_alt)-1], yend = df_alt$y[nrow(df_alt)-1], colour = "alt")) +

    # Critical region
    geom_line(data = df_crit, linetype = 2) +
    #geom_segment(aes(x = b_alt, xend = b_alt, y = (a >= 0)*y_up + (a < 0)*y_lo, yend = c_y, colour = "alt")) +
    #geom_segment(aes(x = c_x, xend = (a >= 0)*x_up + (a < 0)*x_lo, y =  df_alt$y[nrow(df_alt)-1], yend = df_alt$y[nrow(df_alt)-1], colour = "alt")) +

    scale_color_manual(name = "Hypothesis", breaks = c("null", "alt"),
                       values = c(2,4), labels = c("N", "A")) +

    xlab(expression(z[1])) + ylab(expression(z[2])) +
    coord_fixed() +
    theme_minimal()
}
