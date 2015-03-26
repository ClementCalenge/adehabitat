\name{storemapattr}
\alias{storemapattr}
\alias{mapattr}
\title{Store Attributes of Maps of Class asc and kasc}
\description{
  \code{storemapattr} stores attributes of maps of class \code{asc} and
  \code{kasc} in an object of class \code{mapattr}. 
}
\usage{
storemapattr(x)
}
\arguments{
  \item{x}{an object of class \code{asc} or \code{kasc} }
}
\details{
  This function is essentially used by programmers in functions dealing with
  maps of class \code{asc} or \code{kasc}.  The function
  \code{getXYcoords} is an example of function using \code{storemapattr}.
}
\value{
  Returns an object of class \code{mapattr}.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{ \code{\link{kasc}} and \code{\link{asc}} for
  additional information on the classes \code{kasc} and \code{asc}. }
\examples{
data(puechabon)
(kasc <- puechabon$kasc)

(toto <- storemapattr(kasc))

}
\keyword{programming}
