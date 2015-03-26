\name{gnesfa}
\alias{gnesfa}
\alias{print.gnesfa}
\title{General Niche-Environment System Factor Analysis}
\description{
  The function \code{gnesfa} allows to perform a general
  niche-environment system factor analysis.
}
\usage{
gnesfa(dudi, Focus, Reference,
       centering = c("single", "twice"),
       scannf = TRUE, nfFirst = 2, nfLast = 0)
\method{print}{gnesfa}(x, \dots)
}
\arguments{
  \item{dudi}{an object of class \code{dudi}}
  \item{Focus}{a vector containing the focus weights}
  \item{Reference}{a vector containing the reference weights}
  \item{centering}{a character string indicating the type of centering
    (see details)}
  \item{scannf}{a logical value indicating whether the eigenvalues bar
    plot should be displayed}
  \item{nfFirst}{the number of first axes to be kept}
  \item{nfLast}{the number of last axes to be kept}
  \item{x}{an object of class GNESFA}
  \item{\dots}{further arguments to be passed to other functions}
}
\details{
  The GNESFA is an algorithm which generalises several factor analyses
  of the ecological niche.  A table X gives the values of
  P environmental variables in N resource units (e.g. the pixels of a
  raster map).  A distribution of weights D describes the availability
  of the resource units to the species (if not specified, these weights
  are considered to be uniform).  Another distribution of weights Dp
  describes the use of the resource units by the species (for example
  the proportion of relocations in each pixel of a raster map).\cr

  Each environmental variable defines a dimension in a multidimensional
  space, the ecological space. The N resource units define a cloud of
  points in this space.  Each point is associated to two weights.
  The GNESFA finds, in the ecological space, the directions on which
  these two distributions of weights are the most different.\cr
  
  The GNESFA relies on a choice of the analyst, followed by three
  steps.  Before all, the analyst has to choose one distribution of
  weights as the Reference distribution, and the other one as the Focus
  distribution;  (i) The first table X is centred on the centroid of the
  Reference distribution;  (ii) a principal component analysis of this
  Reference distribution is performed;  (iii) the cloud of points is
  distorted, so that the Reference distribution takes a standard
  spherical shape;  (iv) a non centred principal component analysis of
  the Focus distribution allows to identify the directions of the
  ecological space where the two distributions are the most different.\cr

  Depending on the distribution chosen as Reference, this algorithm
  returns results with different meanings (see examples).  This
  algorithm is closely related to several common analyses of habitat
  selection/niche (ENFA, MADIFA, Mahalanobis distances, selection
  ratios, etc.).  The examples below give some examples of the
  mathematical properties of this algorithm.\cr

  
  
  Note that the function takes a parameter named \code{centering}.
  Indeed, two types of centering can be performed prior to the GNESFA.
  The choice \code{"single"} consists in the centering of the cloud of
  point in the ecological space on the centroid of the Reference
  distribution.   The choice \code{"twice"} consist to center the cloud
  of points on both the centroid of the Reference distribution and the
  centroid of the Focus distribution.  This is done by projecting the
  cloud of points on the hyperplane orthogonal to the marginality vector
  (the vector connecting the two centroids.  If this choice is done, the
  GNESFA is identical to the commonly used Ecological Niche Factor
  Analysis (see examples).
}
\value{
  \code{gnesfa} returns a list of class \code{gnesfa} containing the
  following components: 
  \item{call}{original call.}
  \item{centering}{The type of centering required.}
  \item{tab}{a data frame with n rows and p columns.}
  \item{Reference}{a vector of length n containing the Reference
    weights.}
  \item{Focus}{a vector of length n containing the Focus
    weights.}
  \item{nfFirst}{the number of kept first axes.}
  \item{nfLast}{the number of kept last axes.}
  \item{eig}{a vector with all the eigenvalues of the analysis.}
  \item{li}{row coordinates, data frame with n rows and nf columns.}
  \item{l1}{row normed coordinates, data frame with n rows and nf columns.}
  \item{co}{column scores, data frame with p rows and nf columns.}
  \item{cor}{the correlation between the GNESFA axes and the original
    variables}
}
\references{
  Calenge, C. and Basille, M. (in prep.) A General Framework for the
  Exploration of the Ecological Niche.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{madifa}}, \code{\link{mahasuhab}},
  \code{\link{enfa}}, \code{\link{wi}} for closely related methods (see
  Examples)}
\examples{

\dontrun{

################################################################
##
## Study of the habitat selection by the chamois in the French
## mountains of Les Bauges


## Loads the data
data(bauges)
names(bauges)
kasc <- bauges$kasc
locs <- bauges$locs


## displays the data
image(kasc)
image(getkasc(kasc,1))
points(locs, col = "red", pch = 16)


## Prepares the data for the GNESFA:
litab <- kasc2df(kasc)
pc <- dudi.pca(litab$tab, scannf = FALSE)
Dp <- count.points(locs, kasc)[litab$index]

## Example of use with Dp = Reference
gn <- gnesfa(pc, Reference = Dp, scannf=FALSE)

## One main axis:
barplot(gn$eig)

## The correlation with variables indicate that
## the elevation, the proximity to grass and to
## deciduous forests:
s.arrow(gn$cor)

## The factorial map of the niche...
scatterniche(gn$li, Dp, pts = TRUE)

## The chamois is rather located at high elevation,
## in the grass, far from deciduous forests




##########################################################
##########################################################
##
##
##       Some interesting properties of the GNESFA
##
##
##########################################################
##########################################################



################################
################################
##
## Interesting properties of the
## choice: Dp as Reference
## identical to the MADIFA
## (Calenge et al. in revision,
## See the help page of the function madifa
## for other properties)

gn <- gnesfa(pc, Reference = Dp, scannf=FALSE,
             nfFirst = 7)
gn

## This is the same as the MADIFA:
mad <- madifa(pc, Dp, scannf=FALSE)

## Indeed:
plot(gn$li[,1], mad$li[,1])
cor(gn$li[,1], mad$li[,1])


## And consequently the sum of the squared scores,
## On the axes of the GNESFA...
su <- apply(gn$l1,1,function(x) sum(x^2))


## ... is equal to the Mahalanobis distances between
## the points and the centroid of the niche
## (Clark et al. 1993, see the help page of mahasuhab)

su2 <- mahasuhab(kasc, locs)[litab$index]


## Indeed:
all(su - su2 < 1e-7)
plot(su, su2)




################################
################################
##
## Centering twice is identical to
## the ENFA (Hirzel et al. 2002, see the help
## page of the function enfa)...


#######
##
## ... If Dp is the Reference:

gn <- gnesfa(pc, Reference = Dp, center = "twice", scannf = FALSE)
gn

enf <- enfa(pc, Dp, scannf = FALSE)
plot(enf$li[,2], gn$li[,1])
cor(enf$li[,2], gn$li[,1])

## The first specialization axis of the ENFA
## is the first axis of the GNESFA!


#######
##
## ... If Dp is the Focus:

gn <- gnesfa(pc, Focus = Dp, center = "twice",
             scannf = FALSE, nfFirst = 6)
plot(enf$li[,2], gn$li[,6])
cor(enf$li[,2], gn$li[,6])

## The first specialization axis of the ENFA
## is the last axis of the GNESFA!


#######
##
## Whatever the distribution chosen as Reference,
## projecting the cloud of points on the hyperplane
## orthogonal to the marginality axis, and performing
## a GNESFA in this subspace is identical to an ENFA!


## The marginality axis of the ENFA is identical
## to the component "projmar" of the GNESFA

plot(enf$li[,1],gn$projmar)
cor(enf$li[,1],gn$projmar)




################################
################################
##
## Interesting properties of the
## case: Dp as Focus, one categorical
## variable. Relationships with the selection
## ratios of Manly et al. (1972, see the
## help page of wi)


## For example, take the Elevation, and
## define a factor with 4 levels
elev <- data.frame(el = cut(litab$tab$Elevation, 4))

## Now, compute the complete disjonctive table
dis <- acm.disjonctif(elev)
head(dis)

## Now perform the GNESFA with Dp as Focus:
pc <- dudi.pca(dis, scannf = FALSE)
gn <- gnesfa(pc, Dp, scannf = FALSE, nfFirst = 3)


#######
##
## This analysis is closely related to the concept of
## selection ratios

## Compute the percentage of use of each level:
us <- apply(dis, 2, function(x) sum(x*Dp)/sum(Dp))
av <- apply(dis, 2, function(x) sum(x)/length(x))

## The selection ratios
wi <- widesI(us, av)$wi

## Compute the sum of the eigenvalue
sum(gn$eig)

## Compute the sum of the selection ratios - 1
sum(wi) - 1

## In other words, when the GNESFA (Dp as Focus) is
## applied on only one categorical variable, this
## analysis finds a set of axes which partition the
## sum of the selection ratios so that it is maximum
## on the first axes!!


}

}
\keyword{multivariate}
