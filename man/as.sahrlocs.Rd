\name{as.sahrlocs}
\alias{as.sahrlocs}
\alias{sahrlocs}
\alias{getsahrlocs}
\alias{print.sahrlocs}
\title{Exploratory Analysis of Habitat Selection}
\description{
  \code{as.sahrlocs} creates objects of class \code{sahrlocs}.
  This class has a central place in habitat selection studies relying on
  radio-tracking data.  Niche analysis and K-select analysis can be
  performed using this class of objects.  This class may also be used
  for exploratory purposes.  This class of
  object has three main components: an object of class \code{kasc}
  describing the study area, an object of class \code{kasc} describing
  what is available to the animals, and an object of class \code{kasc}
  describing the relocations of the animals (\code{sahrlocs} = Study
  Area - Home Range - reLOCationS).  \code{getsahrlocs} converts one of the
  components of an object of class \code{sahrlocs} into an object of class
  \code{kasc}.
}
\usage{
as.sahrlocs(mlocs, mhr, msa, descan = NULL)
getsahrlocs(x, what = c("sa", "hr", "locs"))
\method{print}{sahrlocs}(x, \dots)
}

\arguments{
  \item{mlocs}{an object of class \code{kasc} returned by the
    function \code{count.points.id}}
  \item{mhr}{an object of class \code{kasc} returned by the function
    \code{hr.rast}, by the function \code{buffer.ani}, or by any other
    user-defined function}
  \item{msa}{an object of class \code{kasc} describing the study area}
  \item{descan}{an optionnal data frame with the number of rows equal to
    the number of monitored animals.  Each column of this data frame
    gives any type of information on the monitored animals (e.g. sex,
    age, and so on)}
  \item{what}{a character string giving the component of the
    \code{sahrlocs} object to be
    converted.  Should be either \code{"sa"} (study area), \code{"hr"}
    (home ranges) or \code{"locs"} (relocations)}
  \item{x}{an object of class \code{sahrlocs}}
  \item{\dots}{additionnal parameters to be passed to the generic
    function \code{print}}
}
\details{
  The different maps in the input list of objects \code{kasc} should describe
  the same area (same attributes: same \code{cellsize}, same \code{xll},
  same \code{yll} and same dimensions).
}
\value{
  Returns one object of class \code{sahrlocs}, which is a list
  containing the input arguments for the function
  \code{as.sahrlocs}. Objects of class \code{sahrlocs} have the same
  attributes as objects of class \code{kasc} (\code{xll}, \code{yll},
  \code{cellsize}, \code{nrow}, and \code{ncol}).\cr
  
  \code{getsahrlocs} returns
  an object of class \code{kasc} (see \code{\link{kasc}}).
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{\code{\link{kasc}} for additionnal information on objects
  of class \code{kasc}, \code{\link{hr.rast}} and \code{\link{buffer.ani}}
  for the creation of the \code{"hr"}
  component of this object, \code{\link{count.points.id}} for the
  creation of the \code{"locs"} component of this
  object, \code{\link{plot.sahrlocs}} and \code{\link{image.sahrlocs}}
  for a graphical display of such objects. }

\examples{
data(puechabon)
kasc <- puechabon$kasc
locs <-  puechabon$locs

## Computes the home ranges of the animals...
cp <- mcp(locs[,4:5], locs[,1])

## ... and converts it to raster
cprast <- hr.rast(cp, kasc)
locrast <- count.points.id(locs[,4:5], locs[,1], kasc)

## Creation of the sahrlocs object
(sahr <- as.sahrlocs(locrast, cprast, kasc))

## adds information on the monitored animals
age <- factor(tapply(locs[,2], locs[,1], mean))
sex <- factor(tapply(locs[,3], locs[,1], mean))
info <- as.data.frame(cbind(sex, age))
(sahr <- as.sahrlocs(locrast, cprast, kasc, info))


## Gets the "study area" component of the object
toto <- getsahrlocs(sahr)
image(toto)
}
\keyword{spatial}
\keyword{multivariate}
