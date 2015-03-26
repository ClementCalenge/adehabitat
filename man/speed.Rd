\name{speed}
\alias{speed}
\title{Computes the Speed Between Successive Relocations of an Animal
  - Deprecated}
\description{
  \code{speed} measures the speed between successive relocations of
  animals, using objects of class \code{traj}.
}
\usage{
speed(x, id = levels(x$id), burst = levels(x$burst), date = NULL,
      units = c("seconds", "hours", "days"))
}
\arguments{
  \item{x}{an object of class \code{traj}}
  \item{id}{a character vector giving the identity of the animals for
    which the speed is to be computed}
  \item{burst}{a character vector giving the identity of the circuits
    for which the speed is to be computed (see \code{traj})}
  \item{date}{a vector of class \code{POSIXct} of length 2 (beginning, end)
    delimiting the period of interest}
  \item{units}{a character string. It determines how the
    speeds are computed.  For example, if the coordinates are given in
    meters, and if \code{units = "seconds"}, speeds are returned in
    meters per second.}
}
\value{
  Returns a data frame with the following components:
  \item{id}{the identity of the animal}
  \item{x}{if the speed is computed between the relocation 1 and 2,
    the x coordinate of the relocation 1.}
  \item{y}{if the speed is computed between the relocation 1 and 2,
    the y coordinate of the relocation 1.}
  \item{date}{a vector of class \code{POSIXct}, giving the date of
    relocation 1.}
  \item{burst}{the identity of the circuit}
  \item{sp.x}{the computed speed of the animal in the x direction}
  \item{sp.y}{the computed speed of the animal in the y direction}
  \item{speed}{the computed speed of the animal on the plane.}
  \item{dt}{the duration between the two relocations (in the units given
  by the parameter \code{units}).}
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\note{
  The function \code{speed} is deprecated.  The class \code{ltraj}
  computes the speeds automatically (see \code{ltraj}).
}
\seealso{ \code{\link{traj}}}
\examples{
\dontrun{
#### Computes the speed for each wild boar
#### monitored at Puechabon
data(puechcirc)
puechcirc <- ltraj2traj(puechcirc)
puechcirc

plot(puechcirc)
sp <- speed(puechcirc)
sp[1:4,]
}
}
\keyword{spatial}

