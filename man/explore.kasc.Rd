\name{explore.kasc}
\alias{explore.kasc}
\title{Interactive Exploration of Maps of Class kasc (requires the
  package tkrplot)}
\description{
  This interface allows to explore distances, values, etc. on a map of
  class \code{kasc}.
}
\usage{
explore.kasc(ka, coltxt = "blue", hscale = 1, vscale = 1, \dots)
}
\arguments{
  \item{ka}{An object of class \code{kasc} }
  \item{coltxt}{character.  the color of the text to be printed}
  \item{hscale}{passed to tkrplot}
  \item{vscale}{passed to tkrplot}
  \item{\dots}{not yet implemented}
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{ \code{\link{kasc}} for further information about objects of
class \code{kasc} and \code{\link{image.kasc}} for graphical display}
\examples{
\dontrun{
if (require(tkrplot)) {
 data(puechabon)
 explore.kasc(puechabon$kasc)
}
}
}
\keyword{hplot}
\keyword{spatial}
