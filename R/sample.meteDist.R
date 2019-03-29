#' @title Sample from a `meteR::meteDist` object
#'  
#' @description Draw a random sample of individuals from a `meteDist` object, with or 
#' without replacement.
#' 
#' @details Currently only SAD and IPD classes are supported
#' 
#' @param x a `meteDist` object
#' @param size number of individuals to sample
#' @param replace logical, whether sampling should be with or without replacement
#' 
#' @examples
#' x <- rpois(100, 1) + 1 # this won't look very METEish
#' esf <- meteESF(1:length(x), x)
#' subSAD <- sample.meteDist(sad(esf), size = round(0.5 * sum(x)), replace = TRUE)
#'
#' @author A.J. Rominger <ajrominger@@gmail.com>
#' @export

sample.meteDist <- function(x, size, replace = FALSE) {
    if('sad' %in% class(x)) {
        newx <- rep(1:length(x$data), x$data)
        subx <- sample(newx, size, replace = replace)
        outx <- as.integer(table(subx))
        
        return(sad(meteESF(1:length(outx), outx)))
    } else if('ipd' %in% class(x)) {
        outx <- sample(x$data, size, replace = replace)
        
        return(ipd(meteESF(S0 = x$state.var['S0'], N0 = size, power = outx, minE = x$emin)))
    } else {
        stop('only classes `sad` and `ipd` are supported')
    }
    
}
