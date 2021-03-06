\name{redisltraj}
\alias{redisltraj}
\title{Rediscretization of a Trajectory With Regular Step Length}
\description{
  This functions rediscretizes one or several trajectories in an object
  of class \code{ltraj}.
}
\usage{
redisltraj(l, u, burst = NULL, samplex0 = FALSE, addbit = FALSE,
           nnew = 5)
}
\arguments{
  \item{l}{an object of class \code{ltraj}}
  \item{u}{the new step length}
  \item{burst}{The burst identity of trajectories to be rediscretized.}
  \item{samplex0}{Whether the first relocation of the trajectory should
    be sampled }
  \item{addbit}{logical. Whether the line segment linking the last
    relocation of the rediscretized trajectory and the last relocation
    of the raw trajectory should be added to the result (can be useful
    for computation of fractal dimension)}
  \item{nnew}{optionnally, you may specify the maximum ratio between number of
    relocations of the new trajectory.  If not specified, this maximum
    is equal to 5 times the number of relocations of the raw
    trajectory.}
}
\details{
  The rediscretization of trajectory has been advocated by several
  authors in the literature (Turchin 1998, Bovet & Benhamou 1988).  It
  is also the first step of the computation of the fractal dimension of
  the path (Sugihara & May 1990).  
}
\value{
  An object of class \code{"ltraj"}
}
\references{
  Bovet, P., & Benhamou, S. (1988) Spatial analysis of animal's
  movements using a correlated random walk model. \emph{Journal of
    Theoretical Biology} \bold{131}: 419--433.
  
  Turchin, P. (1998) \emph{Quantitative analysis of movement},
  Sunderland, MA.

  Sugihara, G., & May, R. (1990) Applications of fractals in
  Ecology. \emph{Trends in Ecology and Evolution} \bold{5}: 79--86.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{ltraj}} for further information on objects of
  class \code{ltraj}}
\examples{
data(puechcirc)

puechcirc

## before rediscretization
plot(puechcirc, perani = FALSE)

## after rediscretization
toto <- redisltraj(puechcirc, 100)
plot(toto, perani = FALSE)

}
\keyword{spatial}
