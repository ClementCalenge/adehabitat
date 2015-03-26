\name{plot.area}
\alias{plot.area}
\title{Graphical Display of Objects of Class "area"}
\description{
  \code{plot.area} allows a graphical display of objects of class
  \code{"area"}. 
}
\usage{
\method{plot}{area}(x, which = levels(x[,1]),
          colpol = rep("green", nlevels(x[, 1])),
          colborder = rep("black", nlevels(x[, 1])),
          lwd = 2, add = FALSE, ...) 
}
\arguments{
  \item{x}{an object of class \code{"area"} }
  \item{which}{a vector of character indicating the polygons to be plotted}
  \item{colpol}{a vector of the color for filling the polygon.  The default,
    \code{NA}, is to leave polygons unfilled }
  \item{colborder}{a vector of the color to draw the border.  The
    default.  Use \code{border = NA} to omit borders}
  \item{lwd}{the border width, a \strong{positive} number}
  \item{add}{logical. if \code{TRUE}, the polygons are added to a
    previous plot}
  \item{\dots}{additional arguments to be passed to the generic function
    \code{plot} }
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link[ade4]{area.plot}} for another way to display
  objects of class \code{"area"} }
\examples{

## Loading the relocations of wild boars
## monitored using radio-tracking
data(puechabon)
locs <-  puechabon$locs
el <- getkasc(puechabon$kasc, "Elevation")

## Estimation of the MCP home ranges of the animals
cp <- mcp(locs[,4:5], locs[,1])

## Use of plot.area to display the results
plot(cp)

## different colors:
color <- c("red", "blue", "green", "yellow")
lev <- levels(cp[,1])
image(el)
plot(cp, colpol = color, add = TRUE)
legend(697198, 3165529, legend = lev, fill = color)

## or:
image(el)
plot(cp, colborder = color, colpol = NA, add = TRUE)
legend(697198, 3165529, legend = lev, fill = color)

## plots one animal
image(el)
plot(cp, which = "Brock", add = TRUE)

}
\keyword{hplot}

