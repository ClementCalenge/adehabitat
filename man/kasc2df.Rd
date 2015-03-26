\name{kasc2df}
\alias{kasc2df}
\alias{df2kasc}
\title{Conversion of Objects of Class kasc}
\description{
  An object of class \code{kasc} stores several maps in a data frame
  (one column per variable, and one row per pixel of the raster map).
  However, the features mapped are rarely rectangle-shaped, whereas the
  map are inevitably rectangles.  Therefore, a lot of pixels of the maps
  do not contain data.  The pixels of the map that do not contain data
  are \code{NA} in this data frame.  \code{kasc2df} will "clean" the
  object of class \code{kasc} from these NAs, and will return a data
  frame containing only mapped values that can be used in subsequent
  analysis. \cr
  After these analyses, \code{df2kasc} may be used to convert the
  modified data frame to an object of 
  class \code{kasc} for mapping (e.g. for maps of factorial axes,
  using \code{dudi} analyses, see \code{help(dudi.pca)}).
}
\usage{
kasc2df(x, var = names(x))
df2kasc(df, index, x)
}

\arguments{
  \item{x}{an object of class \code{kasc} in \code{kasc2df}\cr
  an object of class \code{kasc} or \code{mapattr} in \code{df2kasc}}
  \item{var}{a character vector.  The names of the variables in the
    \code{kasc} that are to be kept in the output}
  \item{df}{a data frame resulting from a computation of the component
    \code{tab} of the list previously returned by the function \code{kasc2df}
    (see section \code{Value} below).  This computation may be any form
    of analysis (Principal component analysis, modelling techniques,
    etc.) }
  \item{index}{an integer vector giving the position of the rows of
    \code{df} in the returned \code{kasc} (such an index can be computed
    using \code{kasc2df})}
}
\value{
  \code{kasc2df} returns a list with the following components:
  \item{tab}{a data frame without NAs, with a number of variables equals
    to \code{length(var)}.}
  \item{index}{a vector of indices of the rows of the \code{kasc} kept
    for the analyses (that is, not \code{NA}).}
  \code{df2kasc} returns an object of class \code{kasc}.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{\code{\link{kasc}} for additional information on
  objects of class \code{kasc}.
}

\examples{
data(puechabon)
kasc <- puechabon$kasc

# Display the kasc object
image(kasc)

# Preparation for Principal component analysis
x <- kasc2df(kasc)
x$tab <- x$tab[, (names(x$tab) != "Aspect")]

# Principal component analysis
ana <- dudi.pca(x$tab, scannf = FALSE)

s.corcircle(ana$co)
s.label(ana$li, clab = 0)


## Map of the scores of the rows
scores <- df2kasc(ana$li, x$index, kasc)
image(scores)
}
\keyword{multivariate}
\keyword{NA}
