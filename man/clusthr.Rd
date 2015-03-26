\name{clusthr}
\alias{clusthr}
\alias{print.clusthr}
\alias{plot.clusthr}
\alias{clusthr.area}
\alias{getverticesclusthr}
\title{Estimation of Home Range by Clustering}
\description{
  \code{clusthr} allows the estimation of the home range by clustering
  (see details).\cr
  \code{plot.clusthr} plots the results.\cr
  \code{clusthr.area} computes the home-range size for given percents of
  relocations included in the home range, and plots the results.\cr
  \code{getverticesclusthr} computes the coordinates of the vertices of
  the home ranges (\code{kver.rast} and \code{kver2shapefile} can be
  used to export the results toward GIS, see the help page of these
  functions).\cr
}
\usage{
clusthr(xy, id = NULL)
\method{print}{clusthr}(x, \dots)
\method{plot}{clusthr}(x, whi = names(x), pch = 21,
             bgpts = "white", colpts = "black", cex = 0.7,
             plotit = TRUE, colpol = "grey", \dots)
clusthr.area(x, percent = seq(20, 100, by = 5),
             unin = c("m", "km"), unout = c("ha", "km2", "m2"),
             plotit = TRUE)
getverticesclusthr(x, whi = names(x), lev=95)
}
\arguments{
  \item{xy}{a data frame with two columns containing the coordinates of
    the relocation of the monitored animals}
  \item{id}{a factor giving the identity of the animal for each
    relocation}
  \item{x}{an object of class \code{clusthr} returned by the function
    \code{clusthr}}
  \item{whi}{a vector of character indicating the animals to be plotted}
  \item{pch}{either an integer specifying a symbol or a single character
    to be used as the default in plotting points.  See \code{points}
    for possible values and their interpretation.}
  \item{bgpts}{background ("fill") color for the open plot symbols given by
    \code{pch=21:25}}
  \item{colpts}{character. The color of the points}
  \item{cex}{The size of the points (see \code{help(par)})}
  \item{plotit}{logical.  Whether the plot should be drawn.}
  \item{colpol}{a character string indicating the colors to be used. Can
    be set to \code{"grey"}, \code{"terrain.colors"},
    \code{"heat.colors"}, \code{"cm.colors"}, or \code{NA}}
  \item{percent,lev}{100 minus the proportion of outliers to be excluded
    from the home range }
  \item{unin}{the units of the relocations coordinates.  Either
    \code{"m"} (default) for
    meters or \code{"km"} for kilometers}
  \item{unout}{the units of the output areas.  Either \code{"m2"} for
    square meters, \code{"km2"} for square kilometers or \code{"ha"} for
    hectares (default)}
  \item{\dots}{additional arguments to be passed to the functions
    \code{plot} and \code{print}.}
}
\details{
  This method estimates home range using a modification of
  single-linkage cluster analysis developped by Kenward et
  al. (2001). The clustering process is described hereafter: the three
  locations with the minimum mean of nearest-neighbour joining distances
  (NNJD) form the first cluster. At each step, two distances are
  computed: (i) the minimum mean NNJD between three locations (which
  corresponds to the next potential cluster) and (ii) the minimum of the
  NNJD between a cluster "c" and the closest location.  If (i) is smaller
  that (ii), another cluster is defined with these three locations.
  If (ii) is smaller than (i), the cluster "c" gains a new location.  If
  this new location belong to another cluster, the two cluster
  fuses. The process stop when all relocations are assigned to the same
  cluster.\cr

  At each step of the clustering process, the proportion of all
  relocations which are assigned to a cluster is computed (so that the
  home range can be defined to enclose a given proportion of the
  relocations at hand, i.e. to an uncomplete process).
  At a given step, the home range is defined as
  the set of minimum convex polygon enclosing the relocations in the
  clusters.
}
\value{
  \code{clusthr} returns a list of class \code{clusthr}. This list has one
  component per animal (named as the levels of argument \code{id}). Each
  component is itself a list, with the following sub-components: 
  \item{xy}{the animals' relocations}
  \item{results}{a data.frame with three columns: \code{step}
    indicates the step number of the algorithm, \code{clust} corresponds
    to the cluster assigned to some relocations, and \code{reloc}
    indicates the relocation(s) which is (are) assigned to the cluster
    "clust" at step "step"}
  
  \code{plot.clusthr} returns an invisible list (see \code{invisible}),
  with one component per animal.  Each component is itself a list with
  one component per step of the clustering process. At each step, an
  object of class \code{area} describes the home range of the animal

  \code{clusthr.area} returns a data.frame of class \code{hrsize}, which
  can be plotted using the generic function \code{plot}.
  
  \code{getverticesclusthr} returns a list of class \code{kver}, where each
  element is an object of class \code{area}.
}
\references{
  Kenwward R.E., Clarke R.T., Hodder K.H. and Walls S.S. (2001) Density
  and linkage estimators of homre range: nearest neighbor clustering
  defines multinuclear cores. \emph{Ecology}, \bold{82}, 1905-1920.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{kver}} for further information on objects
  of class \code{kver}.}
\examples{

data(puechabon)
lo<-puechabon$locs[,c("X","Y")]

## Home Range Estimation
res <- clusthr(lo, puechabon$locs$Name)

## Displays the home range
plot(res)

## Computes the home range size
clusthr.area(res)


}
\keyword{spatial}
\keyword{hplot}
