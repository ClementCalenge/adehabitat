\name{rand.kselect}
\alias{rand.kselect}
\alias{print.rand.kselect}
\title{Test of the Third-Order Habitat Selection }
\description{
  \code{rand.kselect} tests whether the marginality vector of animals
  is significantly larger than what is expected under the hypothesis of
  random habitat use (third-order habitat selection:
  selection by the animals of the relocations within their home range;
  the habitat availability is measured for each animal).  The effect of
  each variable  on individual marginality is also tested. 
  Finally, the pertinence of a K-select analysis is also tested.  This is
  a randomisation test.  The alpha-level of the tests is
  ajusted using the Bonferroni inequality.
}
\usage{
rand.kselect(dudi, factor, weight, nrep = 200, alpha = 0.05, ewa = FALSE)
\method{print}{rand.kselect}(x, \dots)
}
\arguments{
  \item{dudi}{an object of class \code{dudi}.}
  \item{factor}{a factor defining the animals identity}
  \item{weight}{a weight vector of integer values (number of relocations
    counted in each resource unit in row of the object dudi)}
  \item{nrep}{the number of repetitions of the test}
  \item{alpha}{the alpha level for the tests.}
  \item{ewa}{logical. If \code{TRUE}, uniform weights are given
    to all animals in the analysis. If \code{FALSE}, animal weights are
    given by the proportion of relocations of each animal (i.e. an
    animal with 10 relocations has a weight 10 times lower than an
    animal with 100 relocations)}
  \item{x}{an object of the class \code{rand.kselect}.}
  \item{\dots}{further arguments to be passed to the generic function 
     \code{print}}
}
\details{
  This test is carried out by simulating a random use of space by
  animals. \code{rand.kselect} is closely related to the function
  \code{kselect} (same  arguments). 

  At each step  of the randomisation procedure, and for each
  animal, the test randomly allocates the nk relocations (where nk is the
  sum of the weight vector for the animal k) in the Ik pixels available
  to this animal (where Ik is the length of the weight vector for animal
  k).

  The length of the marginality vector is 
  recomputed at each step of the randomisation procedure and for each
  animal. The effect of each variable on the use of pixels by each
  animal is measured by the criterion "(average habitat variable j
  used by animal i) minus (average habitat variable j available to
  animal i)". Finally the value of the first eigenvalue of the 
  K-select analysis provides a criterion to test the pertinence of the
  K-select analysis. 
  
  All these values are then compared to the observed values to assess
  the significance of theses effects.
}
\value{
  Returns an object of class \code{rand.kselect}. This list has three
  components: 
  \item{ global }{a vector of length 2 
    giving the results of the randomisation procedure for the first
    eigenvalue of the K-select analysis.}
  \item{ marg }{a matrix giving the significance of the marginality of
    each animal.} 
  \item{ per.ind }{a list giving the results of the randomisation test for the 
    coordinates of the marginality vector for each animal
    on each habitat variable.}
  \item{ alpha }{the alpha level of the tests.}
}
\references{
  Calenge, C., Dufour, A.B. and Maillard,
  D. (\emph{submitted}). K-select analysis, a new method to analyse
  habitat selection in radio-tracking studies.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{ \code{\link{kselect}} to perform a K-select analysis. }

\examples{
\dontrun{
## Loads the data
data(puechabon)
sahr <- puechabon$sahr

## prepares the data for the kselect analysis: only two variables are kept
x <- sahrlocs2kselect(sahr)
tab <- x$tab
dud <- dudi.mix(tab, scannf = FALSE, nf = 2)

## the randomisation tests
## be patient, this can be very long on some machines
(te <- rand.kselect(dud, x$factor, x$weight, nrep = 50))

}
}
\keyword{multivariate}
