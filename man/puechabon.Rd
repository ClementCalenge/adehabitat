\name{puechabon}
\docType{data}
\alias{puechabon}
\title{ Radio-Tracking Data of Wild Boar }
\description{
This data set stores the results of the monitoring of 4 wild boars at
Puechabon (Mediterranean habitat, South of France).
These data have been collected by Daniel Maillard (Office national de la
chasse et de la faune sauvage).
}
\usage{
data(puechabon)
}
\details{
  The list \code{puechabon} has three components:\cr
  \code{puechabon$kasc}
  is an  object of class \code{kasc} (see \code{help(kasc)}) that
  describes several variables on the study area.\cr
  \code{locs} is a data frame of the relocations of the wild boar
  resting sites in summer. Information on
  wild boars is provided by factors \code{Name}, \code{Sex}, \code{Age}.\cr
  \code{sahr} is the associated object of class
  \code{sahrlocs}.  The home ranges were 
  estimated with buffers including all pixels within 500 m of a boar
  relocation.
}
\references{
Maillard, D. (1996). \emph{Occupation et utilisation de la garrigue et
  du vignoble mediterraneens par le Sanglier}.
Universite d'Aix-Marseille III: PhD thesis.
}

\keyword{datasets}
