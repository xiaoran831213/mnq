`%||%` <- function(x, v) if(is.null(x)) v else x

#' Create data frame
#'
#' A wrapper  for R \code{data.frame},  which ensures that string  variables not
#' turn into factors.
#' @noRd
.d <- function(...) data.frame(..., stringsAsFactors=FALSE)

#' Flood objects in a container to an environment
#'
#' \code{flood} assign named  items in a list to the  calling environment, thus,
#' syntax like:
#'
#' a <-  result$a
#' b <-  result$b
#' ...
#' z  <- result$z
#'
#' is simplified to:
#'
#' \code{flood(result)}
#'
#' Be careful with the silent overwriting of existing object.
#'
#' @param x a container with named objects, typically a R-list
#' @param e environment to flood the objects, default to the calling environment
#' @noRd
flood <- function(x, e=parent.frame())
{
    for(n in names(x)) assign(n, x[[n]], e)
    invisible(NULL)
}

#' Collect function calling arguments
#'
#' Put \code{get.arg} inside a function's  body to capture the calling arguments
#' in a single  row \code{data.frame}.  The \code{get.arg} utility  is useful in
#' documenting  the  simulation  configurations  which  usually  are  passed  as
#' function calling arguments, such as the  sample size, the number of variants,
#' the size of noise or heritability, etc.
#'
#' \code{get.arg}  is   not  recommended  for  functions   accepting  non-scalar
#' arguments such as genotype matrix or vector of effects.
#'
#' @return a data.frame of function arguments
#' @noRd
get.arg <- function()
{
    a <- as.list(match.call(sys.function(1), sys.call(1), expand.dots=TRUE))
    a <- lapply(a[-1], function(.)
    {
        switch(class(.), call=deparse(.), name=as.character(.), .)
    })

    ## drop NULL(s) and return
    a <- a[!sapply(a, is.null)]
    do.call(data.frame, c(a, list(stringsAsFactors=FALSE)))
}
