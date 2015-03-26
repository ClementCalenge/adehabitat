"ltraj2spdf" <- function(ltr)
{
    ## Verifications
    if (!inherits(ltr, "ltraj"))
      stop("ltr should be of class \"ltraj\"")

    ## Conversion
    tr <- do.call("rbind", ltr)
    class(tr) <- "data.frame"
    xy <- tr[!is.na(tr$x),c("x","y")]
    tr <- tr[!is.na(tr$x),]
    tr$y <- tr$x <- NULL
    res <- SpatialPointsDataFrame(xy, tr)

    ## Output
    return(res)
  }

