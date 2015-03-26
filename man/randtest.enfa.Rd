\name{randtest.enfa}
\alias{randtest.enfa}
\title{Randomisation Test for the Ecological Niche Factor Analysis }
\description{
\code{randtest.enfa} performs a randomisation test for the Ecological
Niche Factor analysis (ENFA).
}
\usage{
\method{randtest}{enfa}(xtest, nrepet = 999, \dots)
}
\arguments{
  \item{xtest}{an object of class \code{enfa}}
  \item{nrepet}{the number of iterations for the randomisation test}
  \item{\dots}{further arguments to be passed to the generic function
    \code{randtest}}

}
\details{
  This test is carried out by simulating a random distribution of the
  species occurrences in the pixels of a map.

  At each step  of the randomisation procedure, the test randomly
  allocates the nk occurrences (where nk is the 
  sum of the occurrence vector \code{pr} of the object of class \code{enfa})
  in the Ik pixels of the focus area (where Ik is the length of this occurrence
  vector).

  At each step of the procedure, the first eigenvalue of the ENFA
  performed on the randomised data set is recomputed. This
  value provides a criterion to test the pertinence of the ENFA
  analysis.
}
\value{
  returns a list of class \code{randtest}
}
\references{
  Manly, B.F.J. (1997) \emph{Randomization, Bootstrap and Monte Carlo
    Methods in Biology}. London: Chapman & Hall.
  
  Hirzel, A.H., Hausser, J., Chessel, D. and Perrin, N. (2002)
  Ecological-niche factor analysis: How to compute habitat suitability
  maps without absence data? \emph{Ecology}, \bold{83}, 2027--2036.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\seealso{\code{\link{enfa}}}
\examples{
\dontrun{
data(chamois)
locs <- chamois$locs
map <- chamois$map
map$Vegetation <- NULL
da <- data2enfa(map, locs)
en <- enfa(dudi.pca(da$tab, scannf=FALSE), da$pr, scannf = FALSE)
(tutu <- randtest(en, nrepet = 100))
plot(tutu)
}
}
\keyword{multivariate}
