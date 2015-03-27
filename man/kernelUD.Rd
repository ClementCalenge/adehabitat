\name{kernelUD}
\alias{kernelUD}
\alias{print.khr}
\alias{image.khr}
\alias{plotLSCV}
\alias{getvolumeUD}
\alias{kernel.area}
\alias{getverticeshr}
\title{Estimation of Kernel Home-Range}
\description{
  \code{kernelUD} is used to estimate the utilization distribution (UD)
  of animals monitored by radio-tracking, with the classical kernel method.\cr
  \code{getvolumeUD} and \code{kernel.area} provide utilities
  for home-range size estimation.\cr
  \code{getverticeshr} stores the
  home range contour as objects of class \code{area} in a list of class
  \code{kver}, with one component per animal.
}
\usage{
kernelUD(xy, id = NULL, h = "href", grid = 40, same4all = FALSE,
         hlim = c(0.1, 1.5), kern = c("bivnorm", "epa"), extent = 0.5)
\method{print}{khr}(x, \dots)
\method{image}{khr}(x, axes = FALSE, mar = c(0,0,2,0),
          addcontour = TRUE, addpoints = TRUE, \dots)
plotLSCV(x)
getvolumeUD(x)
kernel.area(xy, id, h = "href", grid = 40,
            same4all = FALSE, hlim = c(0.1,1.5), kern = "bivnorm",
            levels = seq(20,95, by = 5),
            unin = c("m", "km"),
            unout = c("ha", "km2", "m2"), extent = 0.5)
getverticeshr(x, lev = 95)
}

\arguments{
  \item{xy}{a data frame with two columns (x and y coordinates of the
    animal relocations)}
  \item{id}{an optional factor giving the animals identity associated to
    \code{xy}}
  \item{h}{a character string or a number.  If \code{h} is set to
    \code{"href"}, the ad hoc method is used for the smoothing parameter
    (see details).  If \code{h} is set to 
    \code{"LSCV"}, the least-square cross validation method is
    used. Note that \code{"LSCV"} is not available if \code{kern = "epa"}. 
    Alternatively, \code{h} may be set to any given numeric value}
  \item{grid}{a number giving the size of the grid on
    which the UD should be estimated.  Alternatively, this parameter may
    be an object of class \code{asc}, or a list of objects of class
    \code{asc}, with named elements corresponding to each level of the
    factor id (see examples)}
  \item{same4all}{logical. If \code{TRUE}, the same grid
    is used for all animals.  If \code{FALSE}, one grid per
    animal is used}
  \item{hlim}{a numeric vector of length two.  If \code{h = "LSCV"},
    the function minimizes the cross-validation criterion for values of
    h ranging from \code{hlim[1]*href} to \code{hlim[2]*href}, where
    \code{href} is the smoothing parameter computed with the ad hoc
    method (see below)}
  \item{kern}{a character string.  If \code{"bivnorm"}, a bivariate
    normal kernel is used.  If \code{"epa"}, an Epanechnikov kernel is
    used.}
  \item{extent}{a value indicating the extent of the grid used for the
    estimation (the extent of the grid on the abscissa is equal to
    \code{(min(xy[,1]) + extent * diff(range(xy[,1])))}).}
  \item{x}{an object of class \code{khr} returned by \code{kernelUD}.}
  \item{axes}{logical.  Whether the axes are to be plotted}
  \item{mar}{the margin parameter (see \code{help(par)})}
  \item{addcontour}{logical. If \code{TRUE}, contours are drawn on the
    graphics}
  \item{addpoints}{logical. If \code{TRUE}, the animal relocations are
    drawn on the graphics}
  \item{levels}{a vector of percentage levels for home-range size
    estimation}
  \item{unin}{the units of the relocations coordinates.  Either \code{"m"} for
    meters (default) or \code{"km"} for kilometers}
  \item{unout}{the units of the output areas.  Either \code{"m2"} for
    square meters, \code{"km2"} for square kilometers or \code{"ha"} for
    hectares (default)}
  \item{lev}{the percentage level for home range contour estimation.}
  \item{\dots}{additionnal parameters to be passed to the generic
    functions \code{print} and \code{image}}
}
\details{
  The Utilization Distribution (UD) is the bivariate function giving the
  probability density that an animal is found at a point according to
  its geographical coordinates.  Using this model, one can define 
  the home range as the minimum area in which an animal has some
  specified probability of being located. The
  functions used here correspond to the approach described in Worton
  (1995).\cr

  The kernel method has been recommended by many authors for the
  estimation of the utilization distribution (e.g. Worton, 1989, 1995).
  The default method for the estimation of the
  smoothing parameter is the \emph{ad hoc} method, i.e. for a bivariate
  normal kernel
  \deqn{h = \sigma n^{- \frac{1}{6}}}{h = Sigma*n^(-1/6)}
  where \deqn{\sigma = 0.5 (\sigma(x)+\sigma(y))}{Sigma =
    0.5*(sd(x)+sd(y))} which supposes that the UD is
  bivariate normal.  If an Epanechnikov kernel is used, this value is
  multiplied by 1.77 (Silverman, 1986, p. 86).
  Alternatively, the smoothing parameter h may be
  computed by Least Square Cross Validation (LSCV).  The estimated value
  then minimizes the Mean Integrated Square Error (MISE), i.e. the
  difference in volume between the true UD and the estimated UD.  Note
  that the cross-validation criterion cannot be minimized in some
  cases.  According to Seaman and Powell (1998) \emph{"This is a difficult
  problem that has not been worked out by statistical
  theoreticians, so no definitive response is available at this
  time"} (see Seaman and Powell, 1998 for further details and tricky
  solutions).  \code{plotLSCV} allows to have a diagnostic of the
  success of minimization of the cross validation criterion (i.e. to
  know whether the minimum of the CV criterion occurs within the scanned
  range).  Finally, the UD is then estimated over a grid.\cr

  The default kernel is the bivariate normal kernel, but the
  Epanechnikov kernel, which requires less computer time is also
  available for the estimation of the UD. \cr
  
  The function \code{getvolumeUD} modifies the UD component of the
  object passed as argument, so
  that the contour of the UD displayed by the functions \code{contour}
  and \code{image.khr} corresponds to the different percentage levels of
  home-range estimation (see examples).  In addition, this function is
  used in the function \code{kernel.area}, to compute the home-range
  size.  Note, that the function \code{plot.hrsize} (see the help page
  of this function) can be used to display the home-range size estimated
  at various levels.
}
\value{
  The class \code{khr} is a class grouping three sub-classes,
  \code{khrud}, \code{kbbhrud} and \code{khrudvol}:
  \code{kernelUD} returns a list of the class \code{khrud}.  This list
  has one component per animal (named as the levels of argument
  \code{id}). Each component is itself a list, with the following
  sub-components: 
  \item{UD}{an object of class \code{asc}, with the values of density
    probability in each cell of the grid}
  \item{h}{if \code{LSCV} is not used, the value of the smoothing
    parameter.  if
    \code{LSCV} is used, a list with three components: \describe{
      \item{\code{CV}}{the results of the cross-validation procedure.
	The first 
	column contains the sequence of values tested for the smoothing 
	parameter, and the second column contains the value of the
	cross-validation criterion.}
      \item{\code{convergence}}{TRUE if the LSCV succeeds (i.e. if the
	optimum smoothing parameter have been found by the procedure),
	FALSE otherwise.}
      \item{\code{h}}{the value of the smoothing parameter used in UD
	estimation.}
    }
  }
  \item{locs}{The relocations used in the estimation procedure.}
  \item{hmeth}{The argument \code{h} of the function kernelUD}

  \code{getvolumeUD} returns a list of class \code{khrvol},
  with the same components as lists of class \code{khrud}.  \cr

  \code{kernel.area} returns a data frame of subclass \code{hrsize},
  with one column per animal and one row per level of
  estimation of the home range.\cr

  \code{getverticeshr} returns an object of class \code{kver}.
}
\references{

  Silverman, B.W. (1986)
  \emph{Density estimation for statistics and data
    analysis}. London: Chapman \& Hall.
  
  Worton, B.J. (1989) Kernel methods for estimating the utilization
  dirstibution in home-range studies. \emph{Ecology}, \bold{70}, 164--168.

  Worton, B.J. (1995) Using Monte Carlo simulation to evaluate
  kernel-based home range estimators. \emph{Journal of Wildlife Management},
  \bold{59},794--800.

  Seaman, D.E. and Powell, R.A. (1998) \emph{Kernel home range estimation
    program (kernelhr)}. Documentation of the
  program.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\link{asc}} for additionnal informations on
  objects of class \code{asc},  \code{\link{mcp}} for estimation of home
  ranges using the minimum convex polygon, and for help on the function
  \code{plot.hrsize}. \code{\link{kver}} for information on objects of
  class \code{kver}, \code{\link{kernelbb}} for an alternative approach
  of the kernel estimation for trajectory data.
}

\examples{

data(puechabon)
loc <- puechabon$locs[, c("X", "Y")]
id <- puechabon$locs[, "Name"]

## Estimation of UD for the four animals
(ud <- kernelUD(loc, id))

image(ud) ## Note that the contours
          ## corresponds to values of probability density
udvol <- getvolumeUD(ud)
image(udvol)
## Here, the contour corresponds to the
## home ranges estimated at different probability
## levels (i.e. the contour 90 corresponds to the 90 percent
## kernel home-range)
## udvol describes, for each cell of the grid,
## the smaller home-range to which it belongs 

## Calculation of the 95 percent home range
ver <- getverticeshr(ud, 95)
elev <- getkasc(puechabon$kasc, "Elevation") # Map of the area
image(elev)
plot(ver, add=TRUE)
legend(696500, 3166000, legend = names(ver), fill = rainbow(4))


## Example of estimation using LSCV
udbis <- kernelUD(loc, id, h = "LSCV")
image(udbis)

## Compare the estimation with ad hoc and LSCV method
## for the smoothing parameter
(cuicui1 <- kernel.area(loc, id)) ## ad hoc
plot(cuicui1)
(cuicui2 <- kernel.area(loc, id, h = "LSCV")) ## LSCV
plot(cuicui2)

## Diagnostic of the cross-validation
plotLSCV(udbis)


## Use of the same4all argument: the same grid
## is used for all animals
udbis <- kernelUD(loc, id, same4all = TRUE)
image(udbis)

## Estimation of the UD on a map
## (e.g. for subsequent analyses on habitat selection)
elev <- getkasc(puechabon$kasc, "Elevation")
opar <- par(mfrow = c(2, 2), mar = c(0, 0, 2, 0))
cont <- getcontour(elev)

for (i in 1:length(udbis)) {
   image(elev, main = names(udbis)[i], axes = FALSE)
   points(udbis[[i]]$locs, pch = 21, bg = "white", col = "black")
}


## Measures the UD in each pixel of the map
udbis <- kernelUD(loc, id, grid = elev)
opar <- par(mfrow = c(2, 2), mar = c(0, 0, 2, 0))
for (i in 1:length(udbis)) {
  image(udbis[[i]]$UD, main = names(udbis)[i], axes = FALSE)
  box()
  polygon(cont[, 2:3])
}
par(opar)




## Estimation of the UD with a list of objects of class "asc" passed as
## argument grid (useful for large datasets)

## For example, consider the following limits:
lim <- rbind(c(697901,701061,3160198,3162604),
             c(698936,701089,3159969,3162518),
             c(698461,701928,3157362,3160427),
             c(698265,701369,3157219,3162661))

gro <- lapply(1:4, function(i) {
              subsetmap(elev, xlim = lim[i,1:2], ylim=lim[i,3:4])
})
names(gro) <- levels(id)

## show the data:
opar <- par(mfrow=c(2,2), mar=c(0.1,0.1,2,0.1))
lapply(1:4, function(i) {
  image(gro[[i]], main=names(gro)[i], axes=FALSE)
  points(loc[id==names(gro)[i],])
  box()
})
gro

## The map has been subset to fit the relocations.
## Now, estimate the UD:
ud.one.per.grid <- kernelUD(loc, id, grid = gro)
image(ud.one.per.grid)


## The UD can then be matched to habitat maps

}
\keyword{spatial}

