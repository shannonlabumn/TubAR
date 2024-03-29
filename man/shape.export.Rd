% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/export.R
\name{shape.export}
\alias{shape.export}
\title{Turn the output from shape.all into a matrix of median trait values}
\usage{
shape.export(shape_dat, remove_tag = TRUE, roundness_thresh = 0.7)
}
\arguments{
\item{shape_dat}{Output from shape.all function}

\item{remove_tag}{If TRUE, objects below the `roundness_thresh` will be removed from the data before taking median values}

\item{roundness_thresh}{The minimum roundness for an object to not be removed from the output, default is 0.7}
}
\value{
A matrix of the median values of BB width, BB length, perimeter, convex perimeter, area, convex hull area, compactness, roundness, and maximum length values for each image.
}
\description{
Turn the output from shape.all into a matrix of median trait values
}
\examples{
  load(system.file("example_data", "shape_exmp.rds", package = "TubAR"))
  shape.export(shape_exmp,remove_tag=TRUE, roundness_thresh=0.7)
}
