\name{engen2008II}
\alias{engen2008II}
\alias{engen2008I}
\alias{print.engenetalI}
\alias{print.engenetalII}
\title{Measuring Habitat Selection Using the Method of Engen et al. (2008) }
\description{
  These functions implements the method described by Engen et al. to
  measure the preference of animals for habitat variables in habitat
  selection studies.
}
\usage{
engen2008II(us, av, id, nsim = 500, nsimra = 500)

engen2008I(us, av, nsimra=500)

\method{print}{engenetalI}(x, \dots)

\method{print}{engenetalII}(x, \dots)

}
\arguments{
  \item{us}{a data frame containing the value of numeric habitat
    variables (columns) in each site (rows) used by the animals.}
  \item{av}{a data frame containing the value of numeric habitat
    variables (columns) in each site (rows) available to the animals.}
  \item{id}{a factor with as many elements as there are rows in
    \code{us}, indicating the ID of the animal that used the
    corresponding rows in \code{us}.}
  \item{nsim}{the number of randomizations used in
    the calculation of the total variance.}
  \item{nsimra}{the number of random allocation of ranks used in the
    calculation of the normal score (see details). }
  \item{x}{an object of class \code{engenetalI} or \code{engenetalII}}
  \item{\dots}{additional arguments to be passed to other functions
    (currently unused)}
}
\details{
  Engen et al. (2008) proposed an original approach to measure the
  preference of animals for values of each particular variable of a
  multivariate set of environmental variables.  Their approach was
  originally developed for the case where there is a sample of used site
  is for each animal in a sample of identified animals (e.g. using
  radiotelemetry or GPS), with several sites per animal (i.e., design II
  according to the classification of Thomas and Taylor, 1990).  However,
  we extended this approach to also include the case where habitat use
  is described by a sample of used site, with one site per unidentified
  animal (i.e., design I).\cr

  The original approach is the following: first, a normal score
  transformation of each habitat variable is performed: for each
  variable, the empirical cumulative distribution is computed, by
  dividing the rank of the value of each available site by the number of
  observations.  Note that the ties are ranked randomly.  Then, the
  inverse of the standard normal integral (see \code{?qnorm}) of the
  cumulative distribution function is computed for all available sites:
  this results into a perfectly normal distribution of the habitat
  variables for the available sites.  Then, the value of the cumulative
  distribution -- estimated from the available sites -- is computed for
  each used site.  Then, the inverse of the standard normal integral is
  computed for each one.\cr

  Engen et al. (2008) suppose the following model describing how habitat
  use results from habitat availability.  Let \deqn{Z_{ij}}{Z_{ij}} be
  the value of a given habitat variable (transformed according to the
  normal score) for the j-th site used by the i-th animal.  Then this
  value can be described by the model:

  \deqn{Z_{ij} = \mu + U_i + V_{ij}}{Z_{ij} = mu + U_i + V_{ij}}

  where \deqn{\mu}{mu} is the preference for the habitat variable (0
  indicates a non-preference), \deqn{U_i}{U_i} and \deqn{V_{ij}}{V_{ij}}
  are normal distributions with means equal to zero and variances equal
  to \deqn{\sigma^2}{sigma^2} and \deqn{\tau^2}{tau^2} respectively.
  Engen et al. give fomula for the estimation of these parameters.
  Their estimation is done by first estimating the total variance
  \deqn{\sigma^2+ \tau^2}{sigma^2 + tau^2} (this variance is estimated
  by sampling randomly one observation per animal -- the parameter
  \code{nsim} controls the number of samples used in this computation;
  see Engen et al. 2008).  Note that the correlation between the value
  observed for two used units sampled from all the units used by a given
  animal is \deqn{\rho = \sigma^2/(\sigma^2 + \tau^2)}{rho = sigma^2 /
  (sigma^2 + tau^2)}.  A large value of rho indicates a large variation
  in the habitat used between animals (or a small within-animal
  variation).  The main parameter of concern is here the preference.
  The function \code{engen2008II} allows to estimate these
  parameters.\cr

  The function \code{engen2008I} extends this model for design I
  studies (a sample of used sites and a sample of available sites,
  animals not identified), by considering the following model for these
  studies:

  \deqn{Z_{ij} = \mu + V_{ij}}{Z_{ij} = mu + V_{ij}}

  where \deqn{\mu}{mu} is the preference for the habitat variable (0
  indicates a non-preference), and \deqn{V_{ij}}{V_{ij}}
  are normal distributions with means equal to zero and variances equal
  to \deqn{\sigma^2}{sigma^2}.\cr
  
  Note that the habitat variables may be correlated on the study
  area. In this case, observed preference for a given variable may be
  an artefact of other variables prefered by animals.  Supposing that
  the data.frame containing the \deqn{Z_{ij}}{Z_{ij}}'s is a realization
  of a multivariate normal distribution, we can compute, for each
  habitat variable and each used site, the *conditional* mean
  \deqn{m_{ij}}{m_{ij}} and *conditional* standard deviation
  \deqn{s_{ij}}{s_{ij}} of this variable at this site, *given* the values
  of the other habitat variables at this site (this is done using the
  algorithm described by Ripley, 1987, p.98).  We then compute the
  standardized values \deqn{P_{ij} = (Z_{ij} - m_{ij})/s_{ij}}{P_{ij} =
  Z_{ij} - m_{ij})/s_{ij}}.  The preference is then computed using these
  standardized values.\cr
  
  Because there may be ties in the distribution of values of habitat
  variables, the results may vary depending on the random order chosen
  for ties when computing the normal scores.  Engen et al. recommended
  to repeat the function a large number of times, and to use the mean
  values as estimates of the parameters.  This is what the function
  does, and the number of randomization is controlled by the parameter
  \code{nsimra}.\cr

  Note that all these methods rely on the following hypotheses: (i)
  independence between animals, (ii) independence between sites, and
  (iii) all animals are selecting habitat in the same way (in addition
  to "traditional" hypotheses in these kinds of studies: no
  territoriality, all animals having equal access to all available
  resource units, etc., see Manly et al. 2002 for further details).
  
  That the examples below provide an illustration and discussion of
  interesting and (at first sight) surprising properties of this
  method.
}
\value{
  \code{engen2008I} returns a list of class \code{engenetalI}, and
  \code{engen2008II} returns a list of class \code{engenetalII}.  Both
  types of list contain two elements:
  \item{raw}{this is a list containing one data frame per habitat
    variable, containing the value of the correlation rho (for
    \code{engenetalII} objects), mean preferences and standard error of
    these preferences (columns) for each randomization performed (rows);}
  \item{results}{a data frame containing the mean values over all
    the randomizations, of these parameters (columns) for each habitat
    variable (rows).}
}
\references{
  Engen, S., Grotan, V., Halley, D. and Nygard, T. (2008) An efficient
  multivariate approach for estimating preference when individual
  observations are dependent. \emph{Journal of Animal Ecology},
  \bold{77}: 958--965.

  Thomas, D. and Taylor, E. (1990) Study designs and tests for comparing
  resource use and availability. \emph{Journal of Wildlife Management},
  \bold{54}, 322--330.

  Ripley, B. (1987) Stochastic Simulation. John Wiley and Sons.   
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\note{
  Be patient! these functions can be very long (depending on the number of
  sites and on the value of \code{simra})  
}
\seealso{ \code{\link[ade4]{niche}}, \code{\link{madifa}},
  \code{\link{gnesfa}} for another approach to tackle the study of
  habitat selection.  For categorical variables, see
  \code{\link{kselect}}
}
\examples{

\dontrun{

#################################
#################################
#################################


## Practical use of engen2008II

data(puechabon)
kasc <- puechabon$kasc

## Removes the aspect (no factor allowed in the function)
kasc$Aspect <- NULL

## engen2008II:
avail <- kasc2df(kasc)$tab
use <- join.kasc(puechabon$locs[,c("X","Y")], kasc)
id <- puechabon$locs$Name

## This function can be very long:
engen2008II(use, avail, id, nsimra=10)




## Practical use of engen2008I

data(lynxjura)
ka <- lynxjura$map
lo <- lynxjura$locs[,1:2]
av <- kasc2df(ka)$tab
us <- join.kasc(lo, ka)

## Idem, be patient here:
engen2008I(us, av, nsimra=10)




#################################
#################################
#################################
##
##  For a deeper discussion on
##  this method... a simulation:



#################################
##
## First, simulation of a dataset
## copy and paste this part into R,
## but skip the reading of the
## comments if you are not interested
## into this simulation

## simulate the available points

set.seed(235)
av <- cbind(rnorm(1000, mean=0, sd=3), rnorm(1000, mean=0, sd=0.5))
tt <- cbind(c(cos(-pi/4),sin(-pi/4)),c(cos(pi/4), sin(pi/4)))
av <- as.data.frame(as.matrix(av)\%*\%tt)

## simulate the used points: we simulate a selection on the first
## principal component of the PCA of the data.frame describing the
## availability. In other words, we simulate the case where the
## habitat selection occurs on the "common part" of the two habitat
## variable (no preference for one particular variable).

us <- do.call("rbind", lapply(1:5, function(i) {
    us1 <- cbind(rnorm(30, mean=rnorm(1, -4, 1), sd=0.5),
                 rnorm(30, mean=rnorm(1, 0, 0.5), sd=0.2))

    return(us1\%*\%tt)
}))
colnames(us) <- colnames(av) <- c("var1", "var2")
id <- gl(5,30)


#################################
##
## Study of the habitat selection on these data
## The data are:
## - us: a matrix containing the used sites for two
##       habitat variables
## - av: a matrix containing the available sites for two
##       habitat variables
## - id: a vector containing the id of 5 animals

## First illustrate the use and availability of the two variables:

plot(av, xlab="Habitat variable 1", ylab="Habitat variable 2",
     col="grey", pch=16)
lius <- split(as.data.frame(us), id)
junk <- lapply(1:5, function(i) points(lius[[i]], pch=16, col=i))

## -----> ***It is very clear that there is a selection***:
## the animals select the low values of both habitat variables.
## (this is what we actually simulated


## Now perform the method of Engen et al. (2008):
engen2008II(us, av, id)


## Surprisingly, the method seems to fail to identify the clear
## habitat selection identified graphically...
##
## In fact, it does not fail:
## this method identifies the part of habitat selection that is clearly
## attributable to a given variable.  Here the animals select the
## the common factor expressed in the two variables, and it is impossible
## to identify whether the selection is due only to the variable 1 or to
## the variable 2: it is caused by both variable simultaneously.
## Once the selection on the variable 2 (including the common part)
## has been removed, there is no longer appearant selection on
## variable 1.  Once the selection caused by the variable 1
## (including the common part) has been removed, there is no
## longer selection on variable 2...
##
## For this reason, Engen et al. recommended to use this method
## concurrently with other factor analyses of the habitat selection
## such as madifa, kselect, niche (in ade4 package), etc.
##
## Note also the strong correlation between the value of two random
## points used by a given animal. This indicates a strong variability
## among animals...



}


}
\keyword{multivariate}
