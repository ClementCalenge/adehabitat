\name{getascattr}
\alias{getascattr}
\alias{getkascattr}
\title{Copy the Attributes of an Object of Class 'asc' or 'kasc'
  to another Object }
\description{
  \code{getascattr} copies the attributes of an object of class \code{asc}
  to another matrix of the same size.\cr
  \code{getkascattr} performs the same operation for objects of class
  \code{kasc}.
}
\usage{
getascattr(xfrom, xto, type = c("numeric", "factor"), lev = NULL)
getkascattr(xkfrom, xkto)
}

\arguments{
  \item{xfrom}{an object of class \code{asc}}
  \item{xto}{a matrix with the same number of rows and columns
    as \code{xfrom}}
  \item{type}{a character string giving the type of the map
    (\code{"factor"} for maps of categorical 
    variables, and \code{"numeric"} otherwise)}
  \item{lev}{if \code{type = "factor"},
    a character vector giving the levels of the mapped variable
    (see \code{help(asc)})}
  \item{xkfrom}{an object of class \code{kasc}}
  \item{xkto}{a data frame with the same number of rows and columns
    as \code{xkfrom}}
}
\value{
  \code{getascattr} returns a raster matrix of class \code{asc},\cr
  \code{getkascattr} returns a data frame of class \code{kasc}
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\link{kasc}} for additionnal information on objects of
  class \code{kasc}, \code{\link{asc}} for additionnal information on
  objects of class \code{asc}}

\examples{
data(puechabon)

## my.map is a map of elevation
my.map <- getkasc(puechabon$kasc, "Elevation")
sl <- getkasc(puechabon$kasc, "Slope")
attributes(sl) <- NULL
sl <- matrix(sl, ncol = ncol(my.map))

## sl is a matrix with the same size as my.map
toto <- getascattr(my.map, sl)
image(toto)

## Same rationale with aspect
asp <- getkasc(puechabon$kasc, "Aspect")
le <- levels(asp)
attributes(asp) <- NULL
asp <- matrix(asp, ncol = ncol(my.map))

## asp is now a matrix with the same size as my.map
tutu <- getascattr(my.map, asp, typ = "factor", lev = le)
cl <- colasc(tutu, NorthEast = "blue", SouthEast = "red", 
           SouthWest = "orange", NorthWest = "green")
image(tutu, clfac = cl)

}
\keyword{programming}
