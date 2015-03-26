\name{hist.kasc}
\alias{hist.kasc}
\title{Histograms of Mapped Variables}
\description{
\code{hist.kasc} performs histograms of the variables mapped in objects
of class \code{kasc}.
}
\usage{
\method{hist}{kasc}(x, type = c("h", "l", "b"), adjust = 1, col,
          border, lwd = 1, \dots)
}
\arguments{
  \item{x}{a raster map of class \code{kasc} }
  \item{type}{what type of plot should be drawn.  Possible types are:\cr
    * \code{"h"} for histograms,\cr
    * \code{"l"} for kernel density estimates (see \code{?density}).\cr
    * \code{"b"} for both histograms and kernel density estimates (see \code{?density}).\cr
    By default, \code{type = "h"} is used.  If \code{type = "l"} is used,
    the position of the mean of each distribution is indicated by dotted
    lines}
  \item{adjust}{if \code{type = "l"}, a parameter used to control the
    bandwidth of the density estimate (see \code{?density})}
  \item{col}{color for the histogram}
  \item{border}{color for the border of the histogram}
  \item{lwd}{if \code{type = "l"}, line width for the density estimate}
  \item{\dots}{further arguments passed to or from other methods}
}
\author{Mathieu Basille \email{basille@ase-research.org}}
\seealso{\code{\link{kasc}}}
\examples{

## Example with factors and numeric variables
data(puechabon)
hist(puechabon$kasc, type = "h")

## Aspect is a factor, then it's not possible to use 
## kernel density estimates for it :
hist(puechabon$kasc)

## Removing the factor Aspect, and smoothing gives :
hist(puechabon$kasc[-2], type = "l")

}
\keyword{hplot}

