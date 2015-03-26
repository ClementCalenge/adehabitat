"asc2im" <- function(x)
{
    ## Verifications
    if (!inherits(x, "asc"))
        stop("should be an object of class \"asc\"")
    if (attr(x, "type")=="factor")
        stop("function not yet implemented for factors")

    ## spatstat needed
    if (!require(spatstat))
        stop("the package spatstat should be available for this function")

    ## Results
    xy<-getXYcoords(x)
    sorties<-spatstat::im(t(unclass(x)), xy$x,xy$y)
    return(sorties)
  }
