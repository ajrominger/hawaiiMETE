#' @title Bootstrap the mean squared error of a `meteR::meteDist` object
#'  
#' @description Produce bootstrap sampling distribution of the MSE of a `meteDist` object
#' 
#' @details Currently only SAD and IPD classes are supported
#' 
#' @param x a `meteDist` object
#' @param B number of bootstrap replicates
#' @param type should the error be calculated on the rank distribution of cumulative distribution
#' @param relative logical, should the error be calculated on a relative scale
#' @param log logical, should the error be calculated on a log scale
#' 
#' @examples
#' x <- rpois(100, 1) + 1 # this won't look very METEish
#' esf <- meteESF(1:length(x), x)
#' bootSAD <- bootMSE.meteDist(sad(esf))
#'
#' @author A.J. Rominger <ajrominger@@gmail.com>
#' @export

bootMSE.meteDist <- function(x, B = 999, type = 'rank', relative = TRUE, log = FALSE) {
    o <- replicate(B, {
        newx <- sample.meteDist(x, x$state.var['N0'], replace = TRUE)
        return(mse(newx, type = type, relative = relative, log = log))
    })
    
    return(c(o, mse(x, type = type, relative = relative, log = log)))
}
