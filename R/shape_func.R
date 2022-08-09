# Cari Schmitz Carley, Shannon Lab, Jan. 2019
# *Tub*er *A*nalysis in *R*

#' Perform `find.shape` for all JPEG files in the working directory.
#'
#' @param n.core The number of processor cores to use in processing. Default is 1.
#' @param pix.min is number of pixels expected for the smallest object to be designated a tuber
#' @param scaledown by which image is divided for faster computing (a reduction in image size)
#' @param colorcard = NULL/"bottomright" to remove a color card if used
#' @param background.color The color of the background of the image. "black", "blue", and "white" are currently supported. Default is white.
#'
#' @return A dataframe containing the image,BB width, BB length, perimeter, convex perimeter, area, convex hull area, compactness, roundness, and maximum length values for each object in all images in the working directory.
#' @examples
#'   setwd(system.file("images", package = "TubAR"))
#'   shape.all(n.core=1)
#' @import future
#' @import future.apply
#' @import boot
#' @export
shape.all <- function(n.core=1, pix.min=4e3, scaledown=4, colorcard="bottomright", background.color="white"){
  files <- list.files( pattern=".*\\.(jpg|JPG)$")
  names <- gsub(".jpg", "", files, ignore.case = T)

	if(n.core > 1){
		plan(multiprocess, workers=n.core)

		results <- future_lapply(files, function(x) find.shape(x, pix.min=pix.min, scaledown=scaledown, colorcard=colorcard, background.color=background.color))

	}else{
			results <- lapply(files, function(x) find.shape(x, pix.min=pix.min, scaledown=scaledown, colorcard=colorcard, background.color=background.color))

		}

	names(results) <- names
	return(unpack(results))

}

#' Unpack shape.all into a dataframe in melted format.
#' @param results Input from shape.all
#' @return A dataframe in melted format
 unpack <- function(results){ # unpack shape.all into a dataframe in melted format
	 t <- unlist(lapply(results, function(x) length(x[[1]])))
	 traits <- names(results[[1]])
	 plots <- names(results)

	 dat <- matrix(NA, nrow=sum(unlist(t)), ncol=length(traits)+1)
	 colnames(dat) <- c("image", traits)

	 for(i in 1:length(traits)){
		dat[, traits[i]] <- unlist(lapply(results, function(x) x[[traits[i]]]))
	 }
	dat[,"image"] <- unlist(sapply(1:length(t), function(x) rep(names(t[x]), t[x])))
	return(dat)
 }

 #' Returns the outline of an image object.
 #' @param Morig A binarized image matrix containing a single object
 #' @return Binarized image matrix of the outline of an object
get.trace <- function(Morig){
	fac <- setdiff(unique(c(Morig)), 0)
	M <- Morig/fac

	horiz <- rep(0, dim(M)[2])
	vert <- rep(0, dim(M)[1])

	W <- cbind(M[ , -c(1)], vert)
	NW <- rbind(W[-c(1), ], horiz)
	N <- rbind(M[-c(1), ], horiz)
	NE <- cbind(vert, N[,-c(dim(N)[2])] )
	E <- cbind(vert, M[ , -c(dim(M)[2])])
	SE <- rbind(horiz, E[-c(dim(E)[1]), ])
	S <- rbind(horiz, M[-c(dim(M)[1]), ])
	SW <- cbind(S[ ,-c(dim(S)[2])], vert)

	M2 <- unname(((W+NW+N+NE+E+SE+S+SW) >0)*1 - M)
	return( (M2 !=0)*fac)

}


#' Calculates the perimeter length from the coordinates of an object outline
#' @param xy An integer vector of points on the convex hull of an object
#' @return The perimeter length of the convex hull in pixels
get.perim <- function(xy){
	perim <- rep(NA, dim(xy)[1])
	perim[1] <- sqrt(abs(xy[length(perim), 1]-xy[1, 1])^2 + abs(xy[length(perim), 2]-xy[1,2])^2)
	for(i in 2:length(perim)){
		perim[i] <- sqrt(abs(xy[i-1, 1]-xy[i, 1])^2 + abs(xy[i-1, 2]-xy[i,2])^2)
	}
	return(sum(perim))
}

#' Gives the length and width of the bounding box that would contain an object.
#' @param xy A matrix of object outline coordinates
#' @return A list of the bounding box endpoint coordinates, width, height, and angle
############################################################################
# Title: Convex hull, (minimum) bounding box, and minimum enclosing circle
# Author: Daniel Wollschlaeger
# Date: 2019
# Source: http://dwoll.de/rexrepos/posts/diagBounding.html
###########################################################################
getMinBBox <- function(xy) {
    stopifnot(is.matrix(xy), is.numeric(xy), nrow(xy) >= 2, ncol(xy) == 2)

    ## rotating calipers algorithm using the convex hull
    H    <- chull(xy)      ## hull indices, vertices ordered clockwise
    n    <- length(H)      ## number of hull vertices
    hull <- xy[H, ]        ## hull vertices

    ## unit basis vectors for all subspaces spanned by the hull edges
    hDir  <- diff(rbind(hull, hull[1, ])) ## hull vertices are circular
    hLens <- sqrt(rowSums(hDir^2))        ## length of basis vectors
    huDir <- diag(1/hLens) %*% hDir       ## scaled to unit length

    ## unit basis vectors for the orthogonal subspaces
    ## rotation by 90 deg -> y' = x, x' = -y
    ouDir <- cbind(-huDir[ , 2], huDir[ , 1])

    ## project hull vertices on the subspaces spanned by the hull edges, and on
    ## the subspaces spanned by their orthogonal complements - in subspace coords
    projMat <- rbind(huDir, ouDir) %*% t(hull)

    ## range of projections and corresponding width/height of bounding rectangle
    rangeH  <- matrix(numeric(n*2), ncol=2)  ## hull edge
    rangeO  <- matrix(numeric(n*2), ncol=2)  ## orthogonal subspace
    widths  <- numeric(n)
    heights <- numeric(n)

    for(i in seq(along=numeric(n))) {
        rangeH[i, ] <- range(projMat[  i, ])

        ## the orthogonal subspace is in the 2nd half of the matrix
        rangeO[i, ] <- range(projMat[n+i, ])
        widths[i]   <- abs(diff(rangeH[i, ]))
        heights[i]  <- abs(diff(rangeO[i, ]))
    }

    ## extreme projections for min-area rect in subspace coordinates
    ## hull edge leading to minimum-area
    eMin  <- which.min(widths*heights)
    hProj <- rbind(   rangeH[eMin, ], 0)
    oProj <- rbind(0, rangeO[eMin, ])

    ## move projections to rectangle corners
    hPts <- sweep(hProj, 1, oProj[ , 1], "+")
    oPts <- sweep(hProj, 1, oProj[ , 2], "+")

    ## corners in standard coordinates, rows = x,y, columns = corners
    ## in combined (4x2)-matrix: reverse point order to be usable in polygon()
    ## basis formed by hull edge and orthogonal subspace
    basis <- cbind(huDir[eMin, ], ouDir[eMin, ])
    hCorn <- basis %*% hPts
    oCorn <- basis %*% oPts
    pts   <- t(cbind(hCorn, oCorn[ , c(2, 1)]))

    ## angle of longer edge pointing up
    dPts <- diff(pts)
    e    <- dPts[which.max(rowSums(dPts^2)), ] ## one of the longer edges
    eUp  <- e * sign(e[2])       ## rotate upwards 180 deg if necessary
    deg  <- atan2(eUp[2], eUp[1])*180 / pi     ## angle in degrees

    return(list(pts=pts, width=widths[eMin], height=heights[eMin], angle=deg))
}


#' Measures the BB width, BB length, perimeter, convex perimeter, area, convex hull area, compactness, roundness, and maximum length values for each object in the image.
#'
#' @param image The image to be processed, must be a JPEG
#' @param pix.min is number of pixels expected for the smallest object to be designated a tuber
#' @param scaledown by which image is divided for faster computing (a reduction in image size)
#' @param colorcard = NULL/"bottomright" to remove a color card if used
#' @param background.color The color of the background of the image. "black", "blue", and "white" are currently supported. Default is white.
#'
#' @return A nested list with the BB width, BB length, perimeter, convex perimeter, area, convex hull area, compactness, roundness, and maximum length values for each object.
#' @examples
#'   find.shape(system.file("images", "2020fy2_317.jpg", package = "TubAR"))
#'   find.shape(system.file("images", "2020fy2_317.jpg", package = "TubAR"), pix.min=4e4, scaledown=5, colorcard="bottomright")
#' @import EBImage
#' @import sp
#' @import future
#' @import future.apply
#' @import boot
#' @export
find.shape <- function(image, pix.min=4e3, scaledown=4, colorcard="bottomright", background.color="white"){

	#read in the image
	im <- readImage(image)

	# reduce image size
	im2 <- resize(im, h=dim(im)[2]/scaledown)
	if(scaledown>1){ # adjust expected tuber size to scaling of image
			pix.min <- pix.min/scaledown
	}

	# grayscale the image
	if(background.color=="black"){
	  gr <- im2@.Data[,,1] # for the black background red seems to allow potatoes of different types to stand out the most
	  bi <- gr < 0.5
	}else if(background.color=="blue"){
	  gr <- im2@.Data[,,3]#blue(B) is obviously the best spectrum to look at here
	  bi <- gr > 0.55 #hard to get a good threshold here with the example photos I had
	}else{
	  gr <- im2@.Data[,,3] # for red potatoes the best separation btwn bkgrnd & tuber is in the B spectrum
	  #gr <- channel(im2, "gray") # an alternative, but for our images more sensitive to reflections on tubers
	  bi <- gr > 0.75
	}

	# fill any holes in the objects
	bifil <- fillHull(1-bi)

	# number the objects and remove those that are smaller than pix.min (soil, shadows, plot numbering, etc)
	lab <- bwlabel(bifil)
	tab <- data.frame(table(lab)); tab$lab <- as.numeric(as.character(tab$lab))
	stab <- tab$lab[which(tab$Freq < pix.min)]

	labels <- lab
	labels[which(lab@.Data %in% stab, arr.ind=T)] <- 0

	# find and remove the color card
	dist <- X <- Y  <-rep(NA, max(labels))
	for(i in 1:max(labels)){
		ix <- which(labels==i, arr.ind=T)
		X[i]=mean(ix[,1])
		Y[i]=mean(ix[,2])

		dist[i] <- sqrt(X[i]^2 + Y[i]^2)
	}

	labels2<-labels

	if(colorcard=="bottomright"){
		labels2[which(labels==which.max(dist), arr.ind=T)] <- 0
	}

	if(colorcard=="topleft"){
		labels2[which(labels==which.min(dist), arr.ind=T)] <- 0
	}

	if(colorcard=="topright" | colorcard=="bottomleft"){
		X2 <- X-dim(labels)[1] +1
		#Y2 <- Y-dim(labels)[2]+1
		dist2 <- sqrt(X2^2 + Y^2)

		if(colorcard=="topright"){
			labels2[which(labels==which.min(dist2), arr.ind=T)] <- 0
		}

		if(colorcard=="bottomleft"){
			labels2[which(labels==which.max(dist2), arr.ind=T)] <- 0
		}
	}

	p <- unique(c(labels2)); p <- setdiff(p, c(0))


	bbox.width<- bbox.height <- perimeter <- convex.perim <- area <-chull.area <-
	roundness <-compactness <- max.length<- ends <-vector("list", length(p))

	labelsX <- medianFilter(labels2, 50/scaledown) # smooth out hyperlocal effects such as peeling skin

  for(i in 1:length(p)){
	  labelsx <- labelsX
	  labelsx[which(labels2!=p[i])] <- 0
	  outline <- get.trace(labelsx)

	# tuber perimeter (in pixel units, sum of units)
		perimeter[[i]] <- sum(outline >0)
		area[[i]] <- sum(labelsx)

	# length - max euclidian distance between coordinates - usually crosses the axis of the potato
		outline2 <- which(outline >0, arr.ind=T)
		dist <- as.matrix(dist(outline2))
		max.length[[i]] <- max(dist)
		#ends[[i]] <- outline2[which(dist==max(dist), arr.ind=T)[1, ], ]

	# compactness
		compactness[[i]] <- (4*pi*area[[i]]) /(perimeter[[i]]^2)

	# circularity (excludes local irregularities)
		ic <- chull(outline2)
		convexhull <- outline2[c(ic, ic[1]), ]
		convex.perim[[i]] <- get.perim(convexhull)

		chull.poly <- Polygon(convexhull, hole=F)
		chull.area[[i]] <- chull.poly@area

		roundness[[i]] <- (4*pi*chull.area[[i]]) /(convex.perim[[i]]^2)

	# bounding box from which to calculate LW ratio
		box <- getMinBBox(outline2)
		bbox.width[[i]] <- min(box$width, box$height)
		bbox.height[[i]] <- max(box$width, box$height)

	}
	return(list(bbox.width=bbox.width, bbox.height=bbox.height, perim=perimeter, convex.perim=convex.perim,
	area=area, chull.area=chull.area, roundness=roundness, compactness=compactness, max.length=max.length)) #, box.ends=ends))

}


