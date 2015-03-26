\name{wi}
\alias{widesI}
\alias{widesII}
\alias{widesIII}
\alias{print.wiI}
\alias{print.wiII}
\alias{print.wiIII}
\alias{plot.wi}
\alias{wi}
\alias{wiI}
\alias{wiII}
\alias{wiIII}
\title{Computation of Selection Ratios for Habitat Selection Studies.}
\description{
  These functions compute the resource selection ratios (wi) for design I, II
  and III data types, with resources defined by several categories.
  Basic tests are also provided.
}
\usage{
widesI(u, a, avknown = TRUE, alpha = 0.05)
widesII(u, a, avknown = TRUE, alpha = 0.05)
widesIII(u, a, avknown = TRUE, alpha = 0.05)
\method{print}{wiI}(x, \dots)
\method{print}{wiII}(x, \dots)
\method{print}{wiIII}(x, \dots)
\method{plot}{wi}(x, caxis = 0.7, clab = 1, ylog = FALSE, errbar = c("CI", "SE"),
        main = "Manly selectivity measure", noorder = TRUE, \dots)
}
\arguments{
  \item{u}{for \code{widesI}, a vector with named elements describing
    the sample of used resource units.  For \code{widesII} and
    \code{widesIII} a matrix or a data frame giving the number of used
    resource units for each animal (in rows) in each resource category
    (in columns)}
  \item{a}{for \code{widesI} and \code{widesII}, a vector with named
    elements describing the sample or the proportion of available
    resource units.  For \code{widesIII} a matrix or a
    data frame giving the number or the proportion of available resource
    units for each animal (in rows) in each resource category (in
    columns)}
  \item{avknown}{logical. \code{TRUE} if the available
    proportions are known, and \code{FALSE} if they are estimated}
  \item{alpha}{the threshold value for the tests and confidence
    intervals}
  \item{x}{an object of class \code{wi}}
  \item{caxis}{character size on axes to be passed to
    \code{par("cex.axis")}}
  \item{clab}{character size of axes labels to be passed to
    \code{par("cex.lab")}}
  \item{ylog}{logical.  If \code{TRUE}, the selection ratios are plotted
    on a log scale}
  \item{errbar}{a character string.  Type of error bars: either \code{"CI"}
    for confidence intervals or \code{"SE"} for standard errors}
  \item{main}{a character string. The title of the graph}
  \item{noorder}{logical. If \code{TRUE}, the habitat categories are
    ordered on the graph in decreasing order of their preference. If
    FALSE, they are not ordered (i.e. they are in the same order as the
    columns in \code{used} and \code{available}}
  \item{\dots}{additionnal arguments to be passed to the function
    \code{plot}} 
}
\details{
  \code{widesI} may be used to explore resource selection by
  animals, when designs I are involved (habitat use and availability are
  measured at the population level - animals are not identified).  The
  function tests habitat selection with the Khi2 of Pearson
  and log-likelihood Khi2 (recommended, see Manly et al. 2003).  The
  Manly selectivity measure (selection ratio = used/available) is
  computed, the preference / avoidance is tested for
  each habitat, and the differences between selection ratios are computed
  and tested. \cr

  \code{widesII} computes the selection ratios with design II
  data (same availability for all animals, but use measured for each one).
  Tests of identical habitat use for all animals, and of habitat
  selection are also provided.\cr

  \code{widesIII} computes the
  selection ratios for design III data (when the use and the
  availability are measured for each animal - 
  see examples on the wild boar below).  Habitat selection is tested
  using a Chi-square for each animal, and the overall habitat selection
  is also tested.\cr

  Note that all these methods rely on the following hypotheses: (i)
  independence between animals, and (ii) all animals are selecting
  habitat in the same way (in addition to "traditional" hypotheses in
  these kinds of studies: no territoriality, all animals having equal
  access to all available resource units, etc., see Manly et al. 2002
  for further details).
}
\value{
  \code{widesI} returns a list of the class \code{wiI}.  \code{widesII}
  returns a list of class \code{wiII}.  \code{widesIII}
  returns a list of class \code{wiIII}.  These objects are all
  inheriting from the class \code{wi}.  They have the following components:
  \item{used.prop}{the proportion of use for each resource type.}
  \item{avail.prop}{the proportion of available resources.}
  \item{wi}{the Manly selectivity measure (selection ratio:
    used/available).}
  \item{se.wi}{the standard error of the selection ratios.}
  \item{comparisons}{a list with the following components:\describe{
    \item{\code{diffwi}}{a matrix with the differences of the selection
      ratios for each pair of resource type.}
    \item{\code{ICdiffupper}}{a matrix containing the upper limit of confidence
      interval on the differences of the selection ratios between each
      pair of resource type.}
    \item{\code{ICdifflower}}{a matrix containing the lower limit of confidence
      interval on the differences of the selection ratios between each
      pair of resource type.}
    \item{\code{signif}}{the ranking matrix, with the sign of the differences
      between the resource type in row and the resource type in
      column.  When the difference is significant, the sign is
      tripled.}
    }
  }
  \item{profile}{the profile of preferences: resource types are sorted
    so that the left type is the most preferred and the right type is
    the most avoided.  Habitats for which the selection ratios are not
    significant are connected by a line.}
  \item{alpha}{the parameter \code{alpha} of this function.}
  \item{avknown}{the parameter \code{avknown} of this function.}
  \item{se.used}{only for designs I, the standard error of
    the proportion of use.}
  \item{se.avail}{only for designs I, the standard error of
    the available proportion.}
  \item{chisquwi}{only for designs I, the results of Chi-Square tests of the
    hypothesis that the selection ratios are in average equals to zero.}
  \item{Bi}{only for designs I, equals to \code{wi/sum(wi)}.}
  \item{Khi2P}{only for designs I, test of random resource use (Pearson
    statistic).}
  \item{Khi2L}{For designs I, test of random resource use
    (Log-likelihood statistic).  For design III, global test of random 
    resource use (Log-likelihood statistic)}
  \item{Khi2L1}{only for designs II, test of identical use of habitat by all
    animals.}
  \item{Khi2L2}{only for designs II, test of overall habitat selection.}
  \item{Khi2L2MinusL1}{only for designs II, test of hypothesis that animals
    are on average using resources in proportion to availability,
    irrespective of whether they are the same or not (= Khi2L2 -
    Khi2L1).}
  \item{wij}{only for designs II and III, a matrix with the selection
    ratios for all animals and all resource categories.}
  \item{ICwiupper}{only for designs II and III, the upper limit of the
    confidence intervals on the selection ratios.}
  \item{ICwilower}{only for designs II and III, the lower limit of the
    confidence intervals on the selection ratios.}
  \item{Khi2Lj}{only for designs III, the test of habitat selection for each
    animal.}
}
\references{
  Manly B.F.J., McDonald L.L., Thomas, D.L., McDonald, T.L. & Erickson,
  W.P. (2003) \emph{Resource selection by animals - Statistical design and
    Analysis for field studies. Second edition}
  London: Kluwer academic publishers.

  Thomas D. L. and Taylor E. J. (1990) Study designs and tests for
  comparing resource use and availability. \emph{Journal of Wildlife
  Management}, \bold{54}, 322--330.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}

\seealso{\code{\link{compana}} for compositional analysis, and
  \code{\link{eisera}} to perform an eigenanalysis of selection
  ratios.}

\examples{

############################
## Example of moose (Manly et al., 2003, p.52)
## Known available proportions on design I data
moose.avail <- c(0.34, 0.101, 0.104, 0.455)
moose.used <- c(25, 22, 30, 40)
names(moose.used) <- c("InBurnInterior",
                     "InBurnEdge",
                     "OutOfBurnEdge",
                     "OutOfBurnFurther")
names(moose.avail) <- names(moose.used)
## Computation of wi
(wiRatio <- widesI(moose.used, moose.avail))

## plot the values of the selection ratios
opar <- par(mfrow=c(2,2))
plot(wiRatio)

par(opar)




############################
## Example of Elk (Manly et al., 2003, p.62)
## Estimated available proportions on design I data
elk.avail <- c(15, 61, 84, 40)
elk.used <- c(3, 90, 181, 51)
names(elk.used) <- c("0\%", "1-25\%", "26-75\%", ">75\%")
names(elk.avail) <- names(elk.used)
## Computation of wi
(wiRatio <- widesI(elk.used, elk.avail, avknown=FALSE))

## plot the values of the selection ratios
opar <- par(mfrow=c(2,2))
plot(wiRatio)

par(opar)




############################
## Example of Bighorn (Manly et al., 2003, p.67)
## Known available proportions on design II data
data(bighorn)
## Computation of wi
(wi <- widesII(bighorn$used, bighorn$availT, alpha = 0.1))

## plot the values of the selection ratios
opar <- par(mfrow=c(2,2))
plot(wi)



############################
## Example of Bighorn (Manly et al., 2003, p.74)
## Estimated available proportions on design II data
## Computation of wi
(wi <- widesII(bighorn$used, bighorn$availE, avknown = FALSE, alpha = 0.1))

## plot the values of the selection ratios
plot(wi)

par(opar)


############################
## Example of Wild boar
## Estimated available proportions on design III data
data(puechdesIII)
used <- puechdesIII$used
available <- puechdesIII$available

## calculation of the selectio ratios
## with sampled availability
(i <- widesIII(used,available, avknown = FALSE, alpha = 0.1))

opar <- par(mfrow = c(2,2))
plot(i)

par(opar)


}
\keyword{hplot}
\keyword{math}
