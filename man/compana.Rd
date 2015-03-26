\name{compana}
\alias{compana}
\alias{print.compana}
\title{Compositional Analysis of Habitat Use}
\description{
\code{compana} performs a classical compositional analysis of habitat
use (Aebischer et al., 1993).
}
\usage{
compana(used, avail, test = c("randomisation", "parametric"),
        rnv = 0.01, nrep = 500, alpha = 0.1)
}
\arguments{
  \item{used}{a matrix or a data frame describing the percentage of use
    of habitats (in columns) by animals (in rows). }
  \item{avail}{a matrix or a data frame describing the percentage of
    availability of habitats (in columns) by animals (in rows). }
  \item{test}{a character string.  If \code{"randomisation"},
    randomisation tests are performed for
    both the  habitat ranking and
    the test of habitat selection.
    If \code{"parametric"}, usual
    parametric tests are performed (chi-square for the test of habitat
    selection and t-tests for habitat ranking).}
  \item{rnv}{the number replacing the 0 values occurring in the
    matrix \code{used}. }
  \item{nrep}{the number of repetitions in the
    randomisation tests. }
  \item{alpha}{the alpha level for the tests. }
}
\details{
  The compositional analysis of habitat use has been recommended by
  Aebischer et al. (1993) for the analysis of habitat selection by
  several animals, when the resources are defined by several categories
  (e.g. vegetation types).\cr

  This analysis is carried out in two steps: first the significance of
  habitat selection is tested (using a Wilks lambda). Then, a ranking
  matrix is built, indicating whether the habitat type in row is
  significantly used more or less than the habitat type in column.
  When this analysis is performed on radio-tracking data, Aebischer et
  al. recommend to study habitat selection at two levels: (i) selection
  of the home range within the study area, and (ii) selection of the
  relocations within the home range.  The first level is termed
  second-order habitat selection on Johnson's scale (1980), and the
  second one, third-order habitat selection.\cr

  When zero values are found
  in the matrix of used habitats, they are replaced by a small value (by
  default, 0.01), according to the recommendations of Aebischer et al.
  (1993).\cr

  When zero values are found in the matrix of available habitats, the
  function \code{compana} uses the procedure termed
  "weighted mean lambda" described in Aebischer et al.
  (1993: Appendix 2), instead of the usual lambda (see examples). Zero
  values can be found in the matrix of available habitats when the
  third-order habitat selection is under focus. In this case, it may
  occur that some habitat types are available to some animals and not to
  the others.\cr

  Note that this method rely on the following hypotheses: (i)
  independence between animals, and (ii) all animals are selecting
  habitat in the same way (in addition to "traditional" hypotheses in
  these kinds of studies: no territoriality, all animals having equal
  access to all available resource units, etc.).  The function
  \code{eisera} can be used as a preliminary to identify whether this is
  indeed the case (see examples). 
}
\value{
  Returns a list of the class \code{compana}:
  \item{used}{the matrix of used habitats}
  \item{avail}{the matrix of available habitats}
  \item{type.test}{a character string.  Either \code{"randomisation"} or
    \code{"parametric"}}
  \item{test}{the results of the test of habitat selection}
  \item{rm}{the ranking matrix: a square matrix with \code{nh} rows and
    \code{nh} columns, where \code{nh} is the number of habitat types
    under study.  At the intersection of the row \code{i} and of the
    column \code{j}, there is a \code{"+"} when the habitat \code{i} is
    more used than the habitat in column, and \code{"-"} otherwise.
    When the difference is significant, the sign is tripled.}
  \item{rmnb}{the matrix containing the number of animals used to
    perform the tests in \code{rm}.}
  \item{rank}{the rank of the habitat types.  It is equal to the number
    of \code{"+"} for each habitat type in row of \code{rm}.}
  \item{rmv}{the matrix of statistics used to build \code{rm}.
    If \code{(test = "parametric")}, the matrix contains the
    values of \code{t}, in the t-test comparing the row and the column
    habitat.  If \code{(test = "randomisation")}, the matrix contains
    the mean difference between the used and available log-ratios (see
    Aebischer et al., 1993).}
  \item{profile}{the profile of preferences: resource types are sorted
    so that the left type is the most preferred and the right type is
    the most avoided.  Habitats for which the intensity of habitat
    selection is similar (no significant difference) are connected by a
    line.}
}
\references{
Aebischer, N. J. and Robertson, P. A. (1992) Practical aspects of
compositional analysis as applied to pheasant habitat
utilisation. pp. 285--293 In: Priede, G. and Swift, S. M. \emph{Wildlife
  telemetry, remote monitoring and tracking of animals}. 

Aebischer, N. J., Robertson, P. A. and Kenward, R. E. (1993)
Compositional analysis of habitat use from animal radiotracking
data. \emph{Ecology}, \bold{74}, 1313--1325.

Johnson, D. H. (1980) The comparison of usage and availability
measurements for evaluating resource preference. \emph{Ecology}, \bold{61},
65--71.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\note{
In the examples below, the results differ from those published in
Aebischer et al. (squirrel example, selection of the relocations within
the home range). In fact, there has been a confusion in the column names
in the paper. Actually, Aebischer (pers. com.) indicated that the
ranking matrix given in this example is correct.
}
\seealso{ \code{\link{eisera}} to perform an eigenanalysis of selection
  ratios, preliminary to the use of compositional analysis. }

\examples{

## The examples presented here
## are the same as those presented in
## the paper of Aebischer et al. (1993)


#############################
## Pheasant dataset: first
## example in Aebischer et al.

data(pheasant)

## Second order habitat selection
## Selection of home range within the
## study area (example of parametric test)
pheana2 <- compana(pheasant$mcp, pheasant$studyarea,
                   test = "parametric")
pheana2

## The ranking matrix:
print(pheana2$rm, quote = FALSE)

## Third order habitat selection
## (relocation within home range)
## We remove the first pheasant of the analysis
## (as in the paper of Aebischer et al.)
## before the analysis
pheana3 <- compana(pheasant$locs[-1,], pheasant$mcp[-1,c(1,2,4)])
pheana3

## The ranking matrix:
print(pheana3$rm, quote = FALSE)



#############################
## Squirrel data set: second
## example in Aebischer et al.

data(squirrel)

## Second order habitat selection
## Selection of home range within the
## study area
squiana2 <- compana(squirrel$mcp, squirrel$studyarea)
squiana2

## The ranking matrix:
print(squiana2$rm, quote = FALSE)


## However, note that here, the hypothesis of identical use 
## on which this analysis relies is likely to be false.
## Indeed, an eisera indicates:

us <- round(30 * squirrel$locs / 100)
av <- squirrel$studyarea
ii <- eisera(us, av, scannf = FALSE)
scatter(ii, grid = FALSE, clab = 0.7)

## There are clearly two groups of animals.  In such cases,
## compositional analysis is to be avoided in this case.




## Third order habitat selection
## (relocation within home range)
## We remove the second column
## (as in the paper of Aebischer et al.)
squiana3 <- compana(squirrel$locs[,-2], squirrel$mcp[,-2])
squiana3

## The ranking matrix:
print(squiana3$rm, quote = FALSE)


}
\keyword{math}
