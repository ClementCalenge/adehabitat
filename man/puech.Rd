\name{puech}
\alias{puech}
\docType{data}
\title{Radio-Tracking Data of Wild Boar (2)}
\description{
  This data set stores the results of the monitoring of 6 wild boar at
  Puechabon (Mediterranean habitat, South of France). These data have
  been collected by Daniel Maillard (Office national de la chasse et de
  la faune sauvage).
}
\usage{data(puech)}
\details{
  The list \code{puech} has two components:\cr
  \code{puech$relocations} is a data frame containing the
  relocations of the wild boar resting sites in summer.  It contains the
  coordinates of the relocations and the name of the corresponding wild
  boar.\cr
  \code{puech$map} is a list of objects of class \code{asc} (see
  \code{help(asc)}) that describe nine environmental variables on the
  study area (the elevation, the tree cover, the shrub
  cover, the distance to recreational trails, the distance to crops,
  the distance to water points, the grass cover, the slope and the
  sunshine).\cr 
  Note that both the maps and the relocations have been slightly
  "damaged" to preserve the copyright on the data.
}
\references{
  Maillard, D. (1996). \emph{Occupation et utilisation de la garrigue et
    du vignoble mediterraneens par le Sanglier}.
  Universite d'Aix-Marseille III: PhD thesis.
}
\keyword{datasets}
