\name{histniche}
\alias{histniche}
\title{Histograms of the Ecological Niche}
\description{
  \code{histniche} draws histograms of the niche-environment system: an
  histogram of the available resource units (environment) is drawn on
  the same graph as an histogram of the used resource units (i.e. the
  niche), for comparison.
}
\usage{
histniche(x, pr, type = c("h", "l"), adjust = 1,
          Acol, Ucol, Aborder, Uborder, Alwd = 1,
          Ulwd = 1, ylim, ncla = 15, ...)
}
\arguments{
  \item{x}{a data frame giving the value of environmental variables
    (columns) in resource units (rows, e.g. the pixels of a
    raster map)}
  \item{pr}{a vector of integers with the same length as \code{nrow(x)}
    (giving for example the number of detections in the pixels)}
  \item{type}{what type of plot should be drawn.  Possible types are:\cr
    * \code{"h"} for histograms,\cr
    * \code{"l"} for kernel density estimates (see \code{?density}).\cr
    By default, \code{type = "h"} is used.  If \code{type = "l"} is used,
    the position of the mean of each distribution is indicated by dotted
    lines}
  \item{adjust}{if \code{type = "l"}, a parameter used to control the
    bandwidth of the density estimate (see \code{?density})}
  \item{Acol}{color for the histograms of the available pixels}
  \item{Ucol}{color for the histograms of the used pixels}
  \item{Aborder}{if \code{type = "h"}, color for the border of the
    histograms of the available pixels (see \code{help(hist.default)})}
  \item{Uborder}{if \code{type = "h"}, color for the border of the
    histograms of the used pixels (see \code{help(hist.default)})} 
  \item{Alwd}{if \code{type = "l"}, line width for the density estimate
    of the available pixels} 
  \item{Ulwd}{if \code{type = "l"}, line width for the density estimate
    of the used pixels} 
  \item{ylim}{the limits for the y axis}
  \item{ncla}{The number of classes of the histogram} 
  \item{\dots}{further arguments passed to or from other methods}
}
\author{
  Mathieu Basille \email{basille@ase-research.org}
}
\examples{
\dontrun{

data(puechabon)
cp <- count.points(puechabon$locs[,c("X","Y")], puechabon$kasc)
puechabon$kasc
li <- kasc2df(puechabon$kasc)
cpi <- c(cp)[li$index]

histniche(li$tab, cpi)
histniche(li$tab, cpi, ty="l")

}
}
\keyword{multivariate}
