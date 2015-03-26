\name{sahrlocs2kselect}
\alias{sahrlocs2kselect}
\title{Preparation of K-Select Analysis}
\description{
  \code{sahrlocs2kselect} converts an object of class sahrlocs into a
  list that contains the arguments needed for a kselect analysis (see
  \code{help(kselect)}).
}
\usage{
sahrlocs2kselect(sahr)
}
\arguments{
  \item{sahr}{an object of class \code{sahrlocs}}
}
\value{
  Returns a list with the following components:
  \item{tab}{a data frame giving the habitat composition of the home range of
    animals}
  \item{factor}{a factor giving for each row of \code{tab} the name of the
    corresponding animal}
  \item{weight}{a numeric vector giving for each row of \code{tab} the number
    of relocations numbered in this cell of the raster map}
}
\references{
  Calenge, C., Dufour, A.B. and Maillard,
  D. (\emph{submitted}). K-select analysis, a new method to analyse
  habitat selection in radio-tracking studies.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{\code{\link{as.sahrlocs}} for additional information on the
  objects of class \code{sahrlocs}, \code{\link{kselect}} for a K-select
  analysis, and \code{\link{kasc2df}} for additional information on
  the index component of the output.} 

\examples{
data(puechabon)
sahr <- puechabon$sahr
s <- sahrlocs2kselect(sahr)

# 1. PCA of the dataset
pc <- dudi.mix(s$tab, scannf = FALSE)

# 2. k-select analysis
kn <- kselect(pc, s$factor, s$weight, scannf = FALSE)


scatter(kn)
}
\keyword{multivariate}


