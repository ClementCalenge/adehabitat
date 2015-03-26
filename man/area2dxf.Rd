\name{area2dxf}
\alias{area2dxf}
\title{Exportation of Areas}
\description{
\code{area2dxf} exports a data frame of class \code{area} in a DXF file.  Such
files can be read by virtually all Geographic Information Systems.
}
\usage{
area2dxf(x, file, lay = 1:nlevels(factor(x[, 1])))
}
\arguments{
  \item{x}{an object of class \code{area}}
  \item{file}{a character string giving the name of the DXF file to be
    created}
  \item{lay}{an optional vector with a length equal to the number of
    polygons in \code{x}. This vector is then stored in the field
    "Layer" of the DXF file (see examples), for GIS mapping}
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\link{area}}, \code{\link[ade4]{area.plot}}}
\examples{

\dontrun{
  ## Loads the dataset elec88 from the package ade4
  data(elec88)
  ar <- as.area(elec88$area)
  area.plot(ar, values=elec88$tab[,1])
  
  ## exportation of the election results toward a dxf file
  area2dxf(ar, file = "Dept", lay = elec88$tab[,1])
  
  ## Removes the file
  file.remove("Dept.dxf")
}



}
\keyword{hplot}

