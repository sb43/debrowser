#' Prepares all2all scatter plots for given datasets. 
#' 
#' @param data, data that have the sample names in the header.
#' @return all2all scatter plots
#' @examples
#'   plot<-all2all(mtcars)
#' 
#' @export
#' 
all2all <- function(data) {
    pairs(log10(data[1:1000, ]), pch = 19, cex = 0.25,
            diag.panel = panel.hist, lower.panel = panel.cor)
}

#' Prepares the historgram for the all2all plot. 
#' 
#' @param x, a vector of values for which the histogram is desired
#' @param ..., any additional params
#' @return all2all histogram plots
#' @examples
#'   panel.hist(1)
#' 
#' @export
#' 
panel.hist <- function(x, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5))
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks
    nb <- length(breaks)
    y <- h$counts
    y <- y / max(y)
    rect(breaks[-nb], 0, breaks[-1], y, col = "red", ...)
}

#' Prepares the correlations for the all2all plot. 
#' 
#' @param x, numeric vector x
#' @param y, numeric vector y
#' @param ..., additional parameters
#' @return all2all correlation plots
#' @examples
#'   panel.cor(c(1,2,3), c(4,5,6))
#' 
#' @export
#' 
panel.cor <- function(x, y, ...) {
    par(new = TRUE)
    cor_val <- cor.test(x, y, method = "spearman",
                        na.rm = TRUE, exact = FALSE)$estimate
    cor_val <- round(cor_val, digits = 2)

    legend("center", cex = 1.5, bty = "n", paste("rho=", cor_val))
}

#' Smoothes the scatterplot for the all2all plot. 
#' 
#' @param x, x coordinates
#' @param y, y coordinates
#' @param ..., any additional params
#' @return all2all smoothed scatter plots
#' @examples
#'   n   <- 10000
#'   x1  <- matrix(rnorm(n), ncol = 2)
#'   x2  <- matrix(rnorm(n, mean = 3, sd = 1.5), ncol = 2)
#'   x   <- rbind(x1, x2)
#'   y   <- rbind(x2, x1)
#'   panel.smoothScatter(x, y)
#' 
#' @export
#' 
panel.smoothScatter <- function(x, y, ...) {
    par(new = TRUE)
    smoothScatter(x, y, nrpoints = 0)
}
