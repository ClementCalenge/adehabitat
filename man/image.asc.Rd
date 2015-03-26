\name{image.asc}
\alias{image.asc}
\alias{persp.asc}
\alias{contour.asc}
\alias{plot.asc}
\title{Displays a Color Image of a Matrix of Class 'asc' }
\description{
  These functions display a raster matrix of class \code{asc}.
}
\usage{
\method{image}{asc}(x, col = gray((240:1)/256), clfac = NULL, \dots)
\method{contour}{asc}(x, \dots)
\method{persp}{asc}(x, \dots)
\method{plot}{asc}(x, \dots)
}
\arguments{
  \item{x}{a matrix of class \code{asc}}
  \item{col}{for maps of type \code{"numeric"}, the colors to be used
    (see \code{help(par)})} 
  \item{clfac}{for maps of type \code{"factor"}, a character vector giving
    the names of colors for each level of the factor (see
    \code{help(colasc)})}
  \item{\dots}{additional arguments to be passed to the generic
    function \code{image}, \code{persp}, \code{contour}, and
    \code{filled.contour} (but see below)}
}

\section{Note}{
The function \code{plot.asc} uses the function
\code{filled.contour}. The output produced by \code{filled.contour} is
actually a combination of two plots; one is the filled contour and one
is the legend.  Two separate coordinate systems are set up for these
two plots, but they are only used internally - once the function has
returned these coordinate systems are lost.  If you want to annotate
the main contour plot, for example to add points, you can specify
graphics commands in the \code{plot.axes} argument of the function
\code{filled.contour} (this argument is to be passed to the function
\code{plot.asc}). 
An example is given below.
}

\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\link{image}},  \code{\link{contour}},
  \code{\link{persp}}, \code{\link{filled.contour}},
  \code{\link{asc}}
}

\examples{

data(puechabon)

# Case of a continuous variable: the elevation
my.map <- getkasc(puechabon$kasc, "Elevation")
image(my.map, main = "Elevation in Puechabon")
contour(my.map, add = TRUE)

# use of the function plot.asc: adding points on the map
plot(my.map, plot.axes = {points(puechabon$locs[,c("X","Y")])},
     main = "Elevation")

# Case of a factor: the aspect
asp <- getkasc(puechabon$kasc, "Aspect")
cl <- colasc(asp, NorthEast = "blue", SouthEast = "red", 
             SouthWest = "orange", NorthWest = "green")

## graphical display
image(asp, clfac = cl, main = "Aspect", xlab = "Lambert X",
      ylab = "Lambert Y") 
legend(706500, 3162000, legend = levels(asp), fill = cl,
       cex = 0.7)


opar<-par(mar = c(0,0,3,0), bg = "slategray")
persp(my.map, scale = FALSE, box = FALSE, border = NA, shade = 0.75, 
      col = "darkolivegreen3", expand = 2, theta = -60, phi = 30, 
      main = "The topography of Puechabon")
par(opar)


}
\keyword{hplot}
