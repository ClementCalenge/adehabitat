\name{setmask}
\alias{setmask}
\title{Applies a Mask on Objects of Class 'asc' or 'kasc'}
\description{
  Applies a mask on objects of class \code{asc} or \code{kasc}. In other
  words, the function creates an object of class \code{asc} or
  \code{kasc}, with \code{NA} for all pixels \code{NA} on the masking map.
}
\usage{
setmask(x, mask)
}
\arguments{
  \item{x}{an object of class \code{asc} or \code{kasc}}
  \item{mask}{an object of class \code{asc}}
}
\value{
  Returns an object of class \code{asc} or \code{kasc}
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{ \code{\link{asc}} for additionnal information on
  objects of class \code{asc}. }

\examples{

data(puechabon)
kasc <- puechabon$kasc
image(kasc)
elev <- getkasc(kasc, "Elevation")
slope <- getkasc(kasc, "Slope")

## ma is the mask: only areas with elevation > 250 m
## are kept
ma <- elev
ma[ma < 250] <- NA
ma <- getascattr(elev, ma)
image(ma)

## The mask is applied on maps of slope
slp <- setmask(slope, ma)
image(slp)

## The mask is applied on all maps in kasc
im <- setmask(kasc, ma)
image(im)

}
\keyword{spatial}

