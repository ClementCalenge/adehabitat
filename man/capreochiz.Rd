\name{capreochiz}
\alias{capreochiz}
\docType{data}
\title{GPS Monitoring of one Roe Deer in Chize (France)}
\description{
  This dataset contains the relocations of a roe deer collected using
  GPS collars in the Chize reserve (Deux-Sevres, France) by the ONCFS 
  (Office national de la chasse et de la faune sauvage).
}
\usage{data(capreochiz)}
\format{
  This dataset is a list containing the relocations in an object of
  class \code{ltraj} names \code{locs} and information on these
  relocations in an object named \code{"info"} (DOP, status, etc.).
}
\source{
  Sonia Said, Office national de la chasse et de la faune sauvage,
  CNERA-CS, 1 place Exelmans, 55000 Bar-le-Duc (France).
}
\examples{
data(capreochiz)

plot(capreochiz$locs)

head(capreochiz$info)

}
\keyword{datasets}
