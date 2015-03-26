\name{enfa}
\alias{enfa}
\alias{hist.enfa}
\alias{print.enfa}
\alias{data2enfa}
\alias{print.dataenfa}
\title{Ecological Niche Factor Analysis}
\description{
  \code{enfa} performs an Ecological Niche Factor Analysis.
  \code{hist.enfa} draws histograms of the row scores or of the initial
  variables of the ENFA.\cr
  \code{data2enfa} prepares data (kasc and localizations) to be analyzed
  by the ENFA.
}
\usage{
enfa(dudi, pr, scannf = TRUE, nf = 1)
\method{hist}{enfa}(x, scores = TRUE, type = c("h", "l"), adjust = 1, Acol, Ucol, 
          Aborder, Uborder, Alwd = 1, Ulwd = 1, \dots)
data2enfa(kasc, pts)
}
\arguments{
  \item{dudi}{a duality diagram, object of class \code{dudi} (see Details)}
  \item{pr}{a vector giving the utilization weights associated to each unit}
  \item{scannf}{logical.  Whether the eigenvalues barplot should be displayed}
  \item{nf}{an integer indicating the number of kept
    specialization axes }
  \item{x}{an object of class \code{enfa}}
  \item{scores}{logical.  If \code{TRUE}, the histograms display
    the row scores of the ENFA.  If \code{FALSE}, they display the
    niche on the environmental variables (in this case, this is equivalent to
    \code{histniche})}
  \item{type}{what type of plot should be drawn.  Possible types are:\cr
     \code{"h"} for histograms,\cr
     \code{"l"} for kernel density estimates (see \code{?density}).\cr
    By default, \code{type = "h"} is used.  If \code{type = "l"} is used,
    the position of the mean of each distribution is indicated by dotted
    lines}
  \item{adjust}{if \code{type = "l"}, a parameter used to control the
    bandwidth of the density estimates (see \code{?density})}
  \item{Acol}{if \code{type = "h"}, a color to be used to fill the histogram
    of the available pixels. if \code{type = "l"}, a color to be used for the
    kernel density estimates of the available pixels}
  \item{Ucol}{if \code{type = "h"}, a color to be used to fill the histogram
    of the used pixels. if \code{type = "l"}, a color to be used for the
    kernel density estimates of the used pixels}
  \item{Aborder}{color for the border of the histograms of the available
    pixels}
  \item{Uborder}{color for the border of the histograms of the used
    pixels}
  \item{Alwd}{if \code{type = "l"}, the line width of the kernel density 
    estimates of the available pixels}
  \item{Ulwd}{if \code{type = "l"}, the line width of the kernel density 
    estimates of the used pixels}
  \item{kasc}{a raster map of class \code{kasc}}
  \item{pts}{a data frame with two columns, giving the coordinates of
    the species locations}
  \item{\dots}{further arguments passed to or from other methods }
}
\details{
  The niche concept, as defined by Hutchinson (1957), considers the
  ecological niche of a species as an hypervolume in the multidimensional
  space defined by environmental variables, within
  which the populations of a species can persist.  The Ecological Niche
  Factor Analysis (ENFA) has been developped by Hirzel et al. (2002) to
  analyse the position of the niche in the ecological space.  Nicolas
  Perrin (1984) described the position of the niche in the n-dimensional
  space using two measures: the M-specialization (hereafter termed marginality)
  and the S-specialization (hereafter termed specialization).  The
  marginality represents the squared distance of the niche barycentre from
  the mean available habitat.  A large specialization corresponds to a
  narrow niche relative to the habitat conditions available to the species.\cr

  The ENFA first extracts an axis of marginality (vector from the average
  of available habitat conditions to the average used habitat conditions).
  Then the analysis extracts successives orthogonal axes (i.e. uncorrelated),
  which maximises the specialization of the species.  The calculations
  used in the function are described in Hirzel et al. (2002).\cr

  The function \code{enfa} can be used on both quantitative variables and
  qualitative variables (though the interpretation of the results of the
  ENFA for qualitative variables is still under research), provided that
  the table containing the values of habitat variables (columns) for
  each resource unit (rows) is correctly transformed
  (e.g. column-centered and standardised for tables containing only
  quantitative variables), and that appropriate column weights are given
  (e.g. the sum of the weights for
  the levels of a factor should be the same as the weight of one
  quantitative variable). Therefore, the function \code{enfa} requires
  that a preliminary multivariate analysis is performed on the table
  (using analysis of the family of duality diagram, e.g. principal
  component analysis or Hill and Smith analysis). The object returned
  by this preliminary analysis contains the appropriate weights and
  transformation of the original data frame. For example, 
  the function \code{dudi.mix} can be used first on the data.frame
  containing the value of both quantitative (e.g. slope, elevation) and
  qualitative habitat variables (e.g. vegetation) for each pixel of a
  raster map. The result of this analysis can then be passed as argument
  to the function \code{enfa} (see examples below).\cr
}
\value{
  \code{enfa} returns a list of class \code{enfa} containing the
  following components: 
  \item{call}{original call.}
  \item{tab}{a data frame with n rows and p columns.}
  \item{pr}{a vector of length n containing the number of points in each
    pixel of the map.}
  \item{nf}{the number of kept specialization axes.}
  \item{m}{the marginality (squared length of the marginality vector).}
  \item{s}{a vector with all the eigenvalues of the analysis.}
  \item{lw}{row weights, a vector with n components.}
  \item{li}{row coordinates, data frame with n rows and nf columns.}
  \item{cw}{column weights, a vector with p components.}
  \item{co}{column coordinates, data frame with p rows and nf columns.}
  \item{mar}{coordinates of the marginality vector.}
  \code{data2enfa} returns a list of class \code{dataenfa} containing the
  following components: 
  \item{tab}{a data frame with n rows and p columns.}
  \item{pr}{a vector of length n containing the number of points in each
    pixel of the map.}
  \item{index}{an integer vector giving the position of the rows of
    \code{tab} in the  initial object of class \code{kasc}.}
  \item{attr}{an object of class \code{mapattr} with the attributes of the
    initial \code{kasc}.} 
}
\references{
  Hutchinson, G.E. (1957) Concluding Remarks. \emph{Cold Spring Harbor
  Symposium on Quantitative Biology}, \bold{22}: 415--427. 
  
  Perrin, N. (1984) Contribution a l'ecologie du genre Cepaea
  (Gastropoda) : Approche descriptive et experimentale de l'habitat et
  de la niche ecologique. These de Doctorat. Universite de Lausanne,
  Lausanne. 
  
  Hirzel, A.H., Hausser, J., Chessel, D. and Perrin, N. (2002)
  Ecological-niche factor analysis: How to compute habitat-suitability
  maps without absence data? \emph{Ecology}, \bold{83}, 2027--2036.

  Basille, M., Calenge, C., Marboutin, E., Andersen, R. and Gaillard,
  J.M. (2008) Assessing habitat selection using multivariate statistics:
  Some refinements of the ecological-niche factor
  analysis. \emph{Ecological Modelling}, \bold{211}, 233--240.

}
\author{Mathieu Basille \email{basille@ase-research.org} }

\seealso{
  \code{\link[ade4]{niche}}, \code{\link{kselect}} for other types of
  analysis of the niche, when several species are under studies,
  \code{\link{niche.test}} to perform a test of the marginality and the
  tolerance of the niche, and \code{\link{scatter.enfa}} to have a
  graphical display of objects of class
  \code{enfa}.  See \code{\link{madifa}} for another factorial analysis of
  the ecological niche.
}
\examples{

data(lynxjura)

map <- lynxjura$map

## We keep only "wild" indices.
tmp <- lynxjura$locs[,4]!="D"
locs <- lynxjura$locs[tmp, c("X","Y")]
hist(map, type = "l")
## The variable artif is far from symetric

## We perform a square root transformation
## of this variable
## We therefore normalize the variable 'artif'
map[,4] <- sqrt(map[,4])
hist(map, type = "l")

## We prepare the data for the ENFA
(dataenfa1 <- data2enfa(map, locs))

## We then perform the PCA before the ENFA
pc <- dudi.pca(dataenfa1$tab, scannf = FALSE)

## The object 'pc' contains the transformed table (i.e.
## centered so that all columns have a mean of 0
## and scaled so that all columns have a variance of 1
## 'pc' also contains the weights of the habitat variables,
## and the weights of the pixels in the analysis

(enfa1 <- enfa(pc, dataenfa1$pr,
               scannf = FALSE))
hist(enfa1)
hist(enfa1, scores = FALSE, type = "l")

## randomization test and scatterplot
\dontrun{
(renfa <- randtest(enfa1))
plot(renfa)
scatter(enfa1)
}
}
\keyword{multivariate}
