\name{labcon}
\alias{labcon}
\title{Labelling Connected Features }
\description{
  This function attributes unique labels to pixels belonging to connected
  features on a map of class \code{asc}.
}
\usage{
labcon(x)
}
\arguments{
  \item{x}{an object of class \code{asc}}
}
\value{
  Returns a matrix of class \code{asc}, of type \code{"factor"}, with a
  number of levels equals to the number of connected components
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{asc}} for further information on the class \code{asc}}

\examples{
data(puechabon)
hr <- getsahrlocs(puechabon$sahr, "hr")
u <- getkasc(hr, "Jean")
image(u)

## numbering of the connected components
p <- labcon(u)
nlevels(p)
image(p)

##  stores the first component
c1 <- p
c1[c1 != 1] <- NA
image(c1)

##  stores the second component
c2 <- p
c2[c2 != 2] <- NA
image(c2)
}
\keyword{spatial}

