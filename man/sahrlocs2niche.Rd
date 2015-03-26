\name{sahrlocs2niche}
\alias{sahrlocs2niche}
\title{OMI Analysis of Radio-Tracking Data}
\description{
  \code{sahrlocs2niche} converts an object of class \code{sahrlocs} into a
  list that contains the arguments needed for a niche analysis
  (function \code{\link[ade4]{niche}} in package \code{ade4}).
}
\usage{
sahrlocs2niche(x, ani = names(x$hr), var = names(x$sa), used = c("hr",
"locs")) 
}
\arguments{
  \item{x}{an object of class \code{sahrlocs}.}
  \item{ani}{a character vector giving the name of the animals in the
    analyses}
  \item{var}{a character vector giving the name of the variables in the
    analyses}
  \item{used}{a character string.  If \code{"hr"}, the cells of the raster
    map that are considered used by the animals are taken from the
    \code{hr} component of the object \code{sahrlocs}, if \code{"locs"},
    the cells of the raster map that are considered used by the animals
    are taken from the \code{locs} component (see
    \code{\link{as.sahrlocs}}))
  }
}
\value{
  Returns a list with three components:
  \item{index}{a vector of integer giving the index of the rows of
    \code{x$sa} kept for the analysis (this component may then be used
    with the function \code{\link{df2kasc}}).}
  \item{tab}{the table on which the \code{dudi} analysis will be processed.}
  \item{y}{a table giving the weights of the niche analysis.}
}
\references{Doledec, S., Chessel, D. and Gimaret, C. (2000) Niche
  separation in community analysis: a new method. \emph{Ecology}, \bold{81},
  2914--1927.}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{\code{\link{as.sahrlocs}} for additional information on the
  objects of class \code{sahrlocs}, \code{\link[ade4]{niche}} for an
  OMI analysis, and \code{\link{kasc2df}} for additional information on
  the \code{index} component of the output. }

\examples{
data(puechabon)
sahr <- puechabon$sahr

s <- sahrlocs2niche(sahr, used = "locs")


# 1. dudi.mix of the dataset
pc <- dudi.mix(s$tab, scannf = FALSE)

# 2. niche analysis
n <- niche(pc, s$y, scannf = FALSE)

plot(n)
}
\keyword{multivariate}
