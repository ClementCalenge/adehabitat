\name{plot.ltraj}
\alias{plot.ltraj}
\title{Graphical Display of an Object of Class "ltraj"}
\description{
  \code{plot.ltraj} allows various graphical displays of the
  trajectories.
}
\usage{
\method{plot}{ltraj}(x, id = unique(unlist(lapply(x, attr, which = "id"))),
           burst = unlist(lapply(x, attr, which = "burst")), asc = NULL,
           area = NULL, xlim = NULL, ylim = NULL, colasc =
           gray((240:1)/256), colpol = "green",  addpoints = TRUE,
           addlines = TRUE, perani = TRUE, final = TRUE, ...)
}
\arguments{
  \item{x}{an object of class \code{ltraj}}
  \item{id}{a character vector containing the identity of the
    individuals of interest}
  \item{burst}{a character vector containing the burst levels of
    interest}
  \item{asc}{an object of class \code{asc}}
  \item{area}{an object of class \code{area}}
  \item{xlim}{the ranges to be encompassed by the x axis}
  \item{ylim}{the ranges to be encompassed by the y axis}
  \item{colasc}{a character vector giving the colors of the map of class
    \code{asc}}
  \item{colpol}{a character vector giving the colors of the polygon contour
    map, when \code{area} is not \code{NULL}}
  \item{addpoints}{logical.  If \code{TRUE}, points corresponding to
    each relocation are drawn}
  \item{addlines}{logical.  If \code{TRUE}, points corresponding to each
    relocation are drawn}
  \item{perani}{logical.  If \code{TRUE}, one plot is drawn for each
    value of \code{id}, and the several bursts are superposed on the
    same plot for a given animal.  If \code{FALSE}, one plot is drawn
    for each value of \code{burst} }
  \item{final}{logical. If \code{TRUE}, the initial and final
    relocations of each burst are indicated in blue and red,
    respectively}
  \item{\dots}{arguments to be passed to the generic
    function \code{plot} }
}
\author{  Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{For further information on the class \code{ltraj},
  \code{\link{ltraj}}. }
\examples{

data(puechcirc)

plot(puechcirc)
plot(puechcirc, perani = FALSE)
plot(puechcirc, id = "JE93", perani = FALSE)

data(puechabon)
ele <- getkasc(puechabon$kasc, "Elevation")
plot(puechcirc, perani = FALSE, asc = ele)

cont <- getcontour(ele)
plot(puechcirc, area = cont)

}
\keyword{spatial}
\keyword{hplot}
