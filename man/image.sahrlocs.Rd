\name{image.sahrlocs}
\alias{image.sahrlocs}
\title{Graphical Display of the Habitat Composition of
  the Home Ranges of Animals Monitored Using Radio-Tracking }
\description{
  \code{image.sahrlocs} allows a gray-level display of the composition
  of home ranges (different colors are used for factors).
  For a given variable, the minimum gray level (default is
  \code{"white"}) and the maximum gray level (default is \code{"black"})
  represents respectively the minimum and the maximum of the variable
  \emph{***on the study area***}.
  
}
\usage{
\method{image}{sahrlocs}(x, ani = names(x$hr), var = names(x$sa),
               mar = c(0, 0, 0, 0), axes = FALSE, dfidxy = NULL, colpts =
               "black", pch = 21, bg = "white", inv = FALSE, cexpts = 0.6, 
               csub = 2, possub = c("bottomleft", "bottomright", "topleft",
               "topright"), \dots)
}
\arguments{
  \item{x}{an object of class \code{sahrlocs}}
  \item{ani}{a character vector giving the names of the
    variables of the \code{"hr"} component (the animals) for which a
    display is wanted}
  \item{var}{a character vector giving the names of the
    variables of the \code{"sa"} component (the habitat variables) for
    which a display is wanted} 
  \item{mar}{the graphical parameter \code{mar} (see
    \code{\link{par}})}
  \item{axes}{logical.  Whether the axes should
    be plotted}
  \item{csub}{the character size for the legend, used with
    \code{par("cex")*csub}}
  \item{possub}{a character string indicating the sub-title position
    (\code{"topleft"}, \code{"topright"}, \code{"bottomleft"},
    \code{"bottomright"})} 
  \item{dfidxy}{an optional data frame with three columns giving the
    identity and the coordinates of the relocations of each animal.  (if
    not \code{NULL}, the relocations of each animal are plotted in its
    home range, see Examples)}
  \item{colpts}{if \code{dfidxy} is not \code{NULL}, the color of the
    points to be used for the plot of the relocations}
  \item{pch}{if \code{dfidxy} is not \code{NULL}, the size of the points
    to be used for the plot of the relocations (see \code{\link{par}})}
  \item{bg}{if \code{dfidxy} is not \code{NULL}, the background color to
    be used for the plot of the relocations (see \code{\link{par}})}
  \item{inv}{by default, lower values of the mapped variables are
    brighter.  If \code{FALSE}, the lower values are darker}
  \item{cexpts}{if \code{dfidxy} is not \code{NULL}, the size of the points}
  \item{\dots}{additionnal parameters to be passed to the generic
    function \code{image}}
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\link{as.sahrlocs}} for additionnal information on
  objects of class \code{sahrlocs} }

\examples{
data(puechabon)
(sahr <- puechabon$sahr)

## Displays all the variables for a given animal
image(sahr, ani = "Chou")

## Displays all the animals for a given variable
image(sahr, var = "Elevation")

## Load and displays the relocations of the animals
locs <-  puechabon$locs[,c(1,4:5)]
image(sahr, var = "Elevation", dfidxy = locs, pch = 21)
}
\keyword{spatial}
\keyword{multivariate}
