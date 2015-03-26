\name{as.traj}
\alias{as.traj}
\alias{traj}
\alias{print.traj}
\alias{summary.traj}
\alias{plot.traj}
\alias{getburst}
\alias{traj2df}
\alias{df2traj}
\title{Working with Trajects in 2D Space}
\description{
  The class \code{traj} is intended to explore trajects of
  animals monitored using radio-tracking.\cr
  \code{as.traj} creates an object of this class.\cr
  \code{summary.traj} returns the number of relocations for each "burst"
  of relocations and each animal.\cr
  \code{plot.traj} allows various graphical displays of the
  trajects.\cr
  \code{getburst} returns an object of class \code{traj} satisfying the
  specified criteria (selection of one focus animal, of a period of interest,
  of special "bursts" (see details)).\cr
  \code{traj2df}, and the reciprocal function \code{df2traj}
  respectively converts an object of class \code{traj} to an object of
  class \code{data.frame}, and conversely.
}
\usage{
as.traj(id, xy, date, burst = id, ...)
\method{print}{traj}(x, ...)
\method{summary}{traj}(object, id = levels(object$id), date = NULL, \dots)
\method{plot}{traj}(x, id = levels(x$id), burst = levels(x$burst), date = NULL,
          asc = NULL, area = NULL,
          xlim = range(x$x), ylim = range(x$y),
          colasc = gray((256:1)/256), colpol = "green",
          addpoints = TRUE, addlines = TRUE,
          perani = TRUE, final = TRUE, \dots)
getburst(x, burst = levels(x$burst),
         id = levels(x$id), date = NULL)
traj2df(x)
df2traj(df)
}
\arguments{
  \item{id}{a factor or a character vector giving for each relocation the
    identity of the individual monitored in \code{as.traj}. \cr
    a character vector containing the identity of the
    individuals of interest in other functions}
  \item{xy}{a data frame with 2 columns containing the x and y
    coordinates of the relocations}
  \item{date}{a vector of class \code{POSIXct}
    giving the date for each relocation in \code{as.traj}.\cr
    a vector of class \code{POSIXct} with length 2, indicating the
    beginning and the end of the period of interest in other
    functions}
  \item{burst}{a factor or a character vector giving the identity of each
    "burst" of relocations in \code{as.traj} (e.g. the circuit id, see
    details).  The burst level needs to be unique (two animals cannot
    have the same burst levels).  \cr
    a character vector containing the burst levels of interest in
    \code{plot.traj} and \code{getburst}} 
  \item{x}{an object of class \code{traj}}
  \item{object}{an object of class \code{traj}}
  \item{asc}{an object of class \code{asc}}
  \item{area}{an object of class \code{area} (see \code{help(area)})}
  \item{xlim}{the ranges to be encompassed by the x axis}
  \item{ylim}{the ranges to be encompassed by the y axis}
  \item{colasc}{a character vector giving the colors of the
    map of class \code{asc}}
  \item{colpol}{a character vector giving the colors of the
    polygon contour map, when \code{area} is not \code{NULL}}
  \item{addlines}{logical.  If \code{TRUE}, lines joining consecutive
    relocations are drawn}
  \item{addpoints}{logical.  If \code{TRUE}, points corresponding to
    each relocation are drawn}
  \item{perani}{logical.  If \code{TRUE}, one plot is drawn for each
    level of the factor \code{id}, and for a given animal, the several
    bursts are superposed on the same plot.  If \code{FALSE}, one plot
    is drawn for each level of the factor \code{burst}}
  \item{final}{logical. If \code{TRUE}, the initial and final
    relocations of each burst are indicated in blue and red, respectively}
  \item{df}{a data frame to be converted to the class \code{traj}} 
  \item{\dots}{other optional vectors containing some variables measured at
    each relocation (e.g. temperature, wind, 
    elevation, etc.) in \code{as.traj}.  \cr
    For other functions, arguments to be passed to
    the generic functions \code{plot}, \code{summary} and \code{print}}
}
\details{
  For a given individual, trajects are often sampled as "bursts"
  of relocations (Dunn and Gipson, 1977).  For example,
  when an animal is monitored using radio-tracking, the data may consist
  of several circuits of activity (two successive relocations on one
  circuit are often highly autocorrelated, but the data from two
  circuits may be sampled at long intervals in time).  These bursts are
  indicated by the factor \code{burst}.
}
\value{
  An object of class \code{traj} is a data frame with one column
  named \code{id}, one column named \code{x}, one column named \code{y},
  one column named \code{date} and one column named \code{burst}.  This
  class therefore inherits from the class \code{data.frame}.
}
\references{
  Dunn, J.E. and Gipson, P.S. (1977) Analysis of radio telemetry data in
  studies of home range. \emph{Biometrics}. \bold{59}, 794--800.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{ \code{\link{as.POSIXct}} and \code{\link{strptime}} for
  additional information of the class \code{POSIX}. }
\examples{

data(puechabon)
locs <- puechabon$locs
locs[1:4,]

### Conversion of the date to the format POSIX
da <- as.character(locs$Date)
da <- as.POSIXct(strptime(as.character(locs$Date),
               "\%y\%m\%d"))


### Creation of the object of class "traj"
(tr <- as.traj(id = locs$Name, xy = locs[,c("X", "Y")],
               date = da))
summary(tr)
plot(tr)

### Displays on maps of the study area
k <- puechabon$kasc
ele <- getkasc(k, "Elevation")
plot(tr, asc = ele)


}
\keyword{spatial}
