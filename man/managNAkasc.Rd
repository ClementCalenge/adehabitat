\name{managNAkasc}
\alias{managNAkasc}
\alias{managNAkasc}

\title{"Cleaning" Objects of Class 'kasc'}
\description{
  An object of class \code{kasc} stores several maps in a data frame
  (one column per variable, and one row per pixel of the raster map).
  However, the features mapped are rarely rectangle-shaped, whereas the
  map are inevitably rectangles.  Therefore, a lot of pixels of the maps
  do not contain data.  The pixels of the map that do not contain data
  are \code{NA} in this data frame. \cr
  It often occurs that several variables are not mapped for exactly the
  same area (that is, some pixels are \code{NA} for some variables, and
  not for others). \code{managNAkasc} will set to \code{NA} all pixels
  that are not mapped for all variables.
}
\usage{
managNAkasc(x)
}

\arguments{
  \item{x}{an object of class \code{kasc}}
}
\value{
  Returns an object of class \code{kasc}.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{\code{\link{kasc}} for additional information on
  objects of class \code{kasc}.
}
\keyword{multivariate}
\keyword{NA}
