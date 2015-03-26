\name{convnum}
\alias{convnum}
\title{Conversion from Factor to Numeric for Raster Map}
\description{
The objects of class \code{kasc} may contain maps of type \code{"numeric"}
(e.g. the elevation) or of type \code{"factor"} (e.g. the type of
vegetation).  With \code{convnum}, factor maps are transformed into a
collection of k maps of type "numeric" (where k is the number of levels
of the factor), with zero means and unit variance (see dudi.mix in the
ade4 package for further details). Additionaly, maps of type "numeric"
are also centered and scaled (see \code{dudi.mix} in the
\code{ade4} package for further details).
}
\usage{
convnum(kasc)
}
\arguments{
  \item{kasc}{an object of class \code{kasc}}
}
\value{
  Returns a list with the following components:
  \item{kasc}{an object of class \code{kasc}}
  \item{weight}{the weights associated with
    each map of the object kasc, so that the collection of maps defining
    a categorical variable has the same weight as a numeric variable
    (see examples).} 
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link[ade4]{dudi.mix}} (package ade4), \code{\link{kasc}}}
\examples{

## loads the data
data(puechabon)
kasc <- puechabon$kasc
image(kasc)

## scales all the variables
toto <- convnum(kasc)
image(toto$kasc)

## the Aspect has four levels:
## four variables have been defined.
toto$weight

## The sum of the weights given to each
## level of aspect is equal to:
toto$weight[2:5]
sum(toto$weight[2:5])

## The same weight is therefore given to the factor variable Aspect and
## to any continuous variable, e.g. the elevation.

}
\keyword{manip}

