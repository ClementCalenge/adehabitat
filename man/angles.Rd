\name{angles}
\alias{angles}
\title{Compute Turning Angles - Deprecated}
\description{
  \code{angles} computes the turning angles (in radians) between
  consecutive moves from an object of class \code{traj}. See examples
  for a clearer definition.
}
\usage{
angles(x, id = levels(x$id), burst = levels(x$burst),
       date = NULL, slsp =  c("remove", "missing"))
}
\arguments{
  \item{x}{an object of class \code{traj} }
  \item{id}{a character vector giving the identity of the animals for which
    the angles are to be computed}
  \item{burst}{a character vector giving the identity of the circuits for which
    the angles are to be computed (see \code{traj})}
  \item{date}{a vector of class \code{POSIXct} of length 2 (beginning, end)
    delimiting the period of interest}
  \item{slsp}{a character string.  If \code{"remove"}, successive
    relocations located at the same place are replaced by a single
    relocation, allowing the computation of the angles.  If
    \code{"missing"}, a missing value is returned for the angles when
    successive relocations located at the same place. }
}
\value{
  Returns a data frame with the following components:
  \item{id}{the identity of the animal}
  \item{x}{the x coordinate of the relocation at which the angle is computed}
  \item{y}{the y coordinate of the relocation at which the angle is computed}
  \item{date}{a vector of class \code{POSIXct}, giving the date at which
    the relocation has been taken}
  \item{burst}{the id of the circuit (see \code{help(traj)})}
  \item{angles}{the turning angles between the successive moves.}
}
\references{
  Turchin, P. (1998) \emph{Quantitative analysis of movement. Measuring
    and modeling population redistribution in animals and plants.}
  Sunderland, Massachusetts: Sinauer Associates.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\note{
  The function \code{angles} is deprecated.  The class \code{ltraj}
  computes the turning angles automatically (see \code{ltraj}).
}
\seealso{ \code{\link{speed}} for computation of movement speeds,
  \code{\link{traj}} for additional information about objects of
  class \code{traj}  }
\examples{

\dontrun{
## loads an object of class "traj"
data(puechcirc)
puechcir <- ltraj2traj(puechcirc)
puechcir

## Gets a part of the trajectory of the wild boar named
## CH93 and draws it
## Also displays the turning angles
toto <- puechcir[2:5,]
plot(toto$x, toto$y, asp = 1, ylim = c(3158300, 3158550),
     pch = 16,
     main = "Turning angles between\nthree consecutive moves",
     xlab="X", ylab="Y")
lines(toto$x, toto$y)
lines(c(toto$x[2], 700217.6),
      c(toto$y[2], 3158310), lty=2)
lines(c(toto$x[3],700289),
      c(toto$y[3],3158546), lty=2)
ang1x <- c(700234.8, 700231.9, 700231, 700233.7, 700238.8, 700243.2)
ang1y <- c(3158332, 3158336, 3158341, 3158347, 3158350, 3158350)
ang2x <- c(700283.3, 700278.8, 700275.4, 700272.4, 700271.2, 700271.6,
           700274.7)
ang2y <- c(3158522, 3158522, 3158520, 3158517, 3158514, 3158508, 3158504)
lines(ang1x, ang1y)
lines(ang2x, ang2y)
text(700216.1, 3158349, expression(theta[1]), cex=2)
text(700247.7, 3158531, expression(theta[2]), cex=2)
text(c(700301, 700231), c(3158399, 3158487),
     c("Beginning", "End"), pos=4)



## Computation of the turning angles with real data
## on wild boars
plot(puechcir)
ang <- angles(puechcir)

## The angles are in the column angles:
ang[1:4,]
}

}
\keyword{spatial}
