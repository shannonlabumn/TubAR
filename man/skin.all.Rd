% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skin_func.R
\name{skin.all}
\alias{skin.all}
\title{Perform `find.skin` for all JPEG files in the working directory.}
\usage{
skin.all(
  n.core = 1,
  display = T,
  mode = "debug",
  write.clean = F,
  pix.min = 4000,
  scaledown = 4,
  colorcard = "bottomright",
  color.correct = T
)
}
\arguments{
\item{n.core}{The number of processor cores to use in processing. Default is 1.}

\item{display}{If TRUE, a black and white image of each object will be created for each image}

\item{mode}{If set to “debug” the display will show numbered objects and skinned area on the objects.}

\item{write.clean}{If TRUE, a JPEG image with backgroung removed will be created for each input image}

\item{pix.min}{The minimum pixel size required for objects to not be removed with the background,}

\item{scaledown}{The amount image sizes will be reduced in order to aid processing speed}

\item{colorcard}{Declares which corner the color card is in. "bottomright" by default. Can also be set to "bottomleft", "topright", and "topleft"}

\item{color.correct}{Determines if pixel values will be color corrected based on the color card. TRUE by default.}
}
\value{
A nested list of median values for each image of rednees, skinning, and lightness and the values forthese traits for every object in every image.
}
\description{
Perform `find.skin` for all JPEG files in the working directory.
}
\examples{
  setwd(system.file("images", package = "TubAR"))
  skin.all(n.core=1, display=T, mode="debug", write.clean=F, pix.min=4e3, scaledown=4, colorcard="bottomright")
}
