\name{plotltr}
\alias{plotltr}
\title{Changes in Traject Parameters Over Time}
\description{
  This function allows a graphical examination of the changes in
  descriptive parameters in objects of class \code{ltraj}
}
\usage{
plotltr(x, which = "dist", ...)
}
\arguments{
  \item{x}{An object of class \code{ltraj}}
  \item{which}{a character string giving any syntactically correct R
    expression implying the descriptive elements in \code{x}.}
  \item{\dots}{additional parameters to be passed to the generic
    function \code{plot}}
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{ \code{\link{ltraj}} for additional information about
  objects of class \code{ltraj}, and \code{\link{sliwinltr}} for a
  sliding window smoothing}
\examples{
data(puechcirc)

plotltr(puechcirc, "cos(rel.angle)")
plotltr(puechcirc, "dist")
plotltr(puechcirc, "dx")

}
\keyword{hplot}
