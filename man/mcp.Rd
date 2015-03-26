\name{mcp}
\alias{mcp}
\alias{plot.hrsize}
\alias{mcp.area}
\title{Estimation of the Home Range Using the Minimum Convex Polygon
  Estimator}
\description{
  \code{mcp} computes the home range of several
  animals using the Minimum Convex Polygon estimator.\cr
  \code{mcp.area} is used for home-range size estimation. \cr
  \code{plot.hrsize} is used to display the home-range size estimated at
  various levels.
}
\usage{
mcp(xy, id, percent = 95)
mcp.area(xy, id, percent = seq(20,100, by = 5),
         unin = c("m", "km"),
         unout = c("ha", "km2", "m2"), plotit = TRUE)
\method{plot}{hrsize}(x, \dots)
}
\arguments{
  \item{xy}{a data frame with two columns containing the coordinates of
    the relocations of the monitored animals}
  \item{id}{a factor giving the identity of the
    animal for each relocation}
  \item{percent}{100 minus the proportion of outliers to be excluded
    from the computation }
  \item{unin}{the units of the relocations coordinates.  Either
    \code{"m"} (default) for
    meters or \code{"km"} for kilometers}
  \item{unout}{the units of the output areas.  Either \code{"m2"} for
    square meters, \code{"km2"} for square kilometers or \code{"ha"} for
    hectares (default)}
  \item{plotit}{logical.  Whether the plot should be drawn.}
  \item{x}{an objet of class \code{hrsize} returned by the function
    \code{mcp.area}, or \code{kernel.area} (see \code{kernelUD()})}
  \item{\dots}{additionnal arguments to be passed to the function
    \code{plot}}
}
\details{
  This function computes the Minimum Convex Polygon estimation after the
  removal of \code{(100 minus percent)} percent of the relocations the
  farthest away from the barycenter of the home range (computed by the
  arithmetic mean of the coordinates of the relocations for each
  animal).\cr
}
\value{
  \code{mcp} returns an object of class \code{area}, with one polygon
  per level of the factor ID.\cr
  \code{mcp.area} returns a data frame of class \code{hrsize}, 
  with one column per animal and one row per level of
  estimation of the home range.
}
\references{
  Mohr, C.O. (1947) Table of equivalent populations of north
  american small mammals. \emph{The American Midland Naturalist},
  \bold{37}, 223-249.}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\link[grDevices]{chull}}, \code{\link{plot.area}} to have a
  graphical display of the home ranges, \code{\link{area}} for
  additionnal information on the class \code{area}, and
  \code{\link{area2dxf}} for further exportation toward a GIS.
  \code{\link[ade4]{s.chull}} for another way to display MCP.}
\examples{
data(puechabon)
locs <-  puechabon$locs

cp <- mcp(locs[,4:5], locs[,1])

## Plot the home ranges
opar <- par(mar = c(0,0,0,0))
area.plot(cp)

## ... And the relocations
points(locs[,4:5], pch = 16, col = as.numeric(locs[,1]))
par(opar)

## Computation of the home-range size:
cuicui1 <- mcp.area(locs[,4:5], locs[,1])
plot(cuicui1)

}

\keyword{hplot}
