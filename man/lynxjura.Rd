\name{lynxjura}
\alias{lynxjura}
\docType{data}
\title{Monitoring of Lynx}
\description{
  This data set stores the results of the monitoring of lynx in the
  French Jura between 1980 and 1999.  These data have been collected by
  the Lynx Network of the french wildlife management office (Office
  national de la chasse et de la faune sauvage).
}
\usage{data(lynxjura)}
\format{
  The list \code{lynxjura} has two components:
  \itemize{
    \item{\code{map}}{an object of class \code{kasc} (see
      \code{help(kasc)}) that describes several variables on the study
      area:  \code{forets} is the density of forests, \code{hydro} is the
      density of rivers, \code{routes} is the density of roads and
      \code{artif} is the distance from urbanized areas.}
    \item{\code{locs}}{a data frame containing the locations of presence
      indices of the lynx.  \code{X} and \code{Y} are the x and y
      coordinates, \code{Date} is the date of the collection of the 
      indice and \code{Type} represents the type of data (\code{C}: alive
      lynx captured, \code{D}: attacks on livestock, \code{E}: prints or
      tracks, \code{F}: feces, \code{J}: hairs, \code{L}: corpse of lynx,
      \code{O}: sightings and \code{P}: attacks on wild prey).}
  }
}
\source{
Vandel, J.M. (2001) \emph{Repartition du Lynx (Lynx lynx) en France (Massif
Alpin, Jurassien et Vosgien). Methodologie d'etude et statut
actuel}. Ecole Pratique des Haute Etudes de Montpellier II:
Dissertation.
}
\keyword{datasets}
