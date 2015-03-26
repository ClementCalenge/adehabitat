\name{bauges}
\alias{bauges}
\docType{data}
\title{Census of chamois (Rupicapra rupicapra) in the Bauges mountains}
\description{
  This data set contains the relocations of 198 chamois groups in the
  Bauges mountains, as well as maps of 7 environmental variables on the
  study area.
}
\usage{data(bauges)}
\details{
  This dataset contains a subsample of the data collected by volunteers
  and professionals working in various French wildlife and forest
  management, from 1994 to 2004 in the wildlife reserve of "les Bauges"
  (French Alps). Note that both the maps and the relocations have been
  slightly destroyed to preserve copyright on the data.
}
\source{
  Daniel Maillard, Office National de la chasse et de la faune sauvage,
  95 rue Pierre Flourens, 34000 Montpellier, France
}
\examples{
data(bauges)

image(bauges$kasc)
image(getkasc(bauges$kasc, "Elevation"))
points(bauges$locs, pch=16)

}
\keyword{datasets}
