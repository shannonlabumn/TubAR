# Cari Schmitz Carley, Shannon Lab, Jan. 2019
# *Tub*er *A*nalysis in *R*

#' get color card, segment and return the 24 colors as RGB.
#'
#' @param image Input JPEG image
#' @param colorcard Declares which corner the color card is in. "bottomright" by default. Can also be set to "bottomleft", "topright", and "topleft"
#' @param scaledown The amount image sizes will be reduced in order to aid processing speed
#' @param pix.min The minimum pixel size required for objects to not be removed with the background,
#'
#' @return An array of the 24 colors as RGB
#' @import Morpho

# get color card, segment and return the 24 colors as RGB
grabcard <- function(image, colorcard="bottomright", scaledown=8, pix.min=4e4){
		#read in the image
	im <- readImage(image)

	# reduce image size
	im2c <- resize(im, h=dim(im)[2]/scaledown)
	if(scaledown>1){ # adjust expected tuber size to scaling of image
			pix.min <- pix.min/scaledown
	}

	# grayscale the image
	gr <- im2c@.Data[,,3] # for red potatoes the best separation btwn bkgrnd & tuber is in the B spectrum
	# gr <- channel(im2c, "gray") # an alternative, but for our images more sensitive to reflections on tubers
	bi <- gr > 0.75

	# fill any holes in the objects
	bifil <- fillHull(1-bi)

	# number the objects and remove those that are smaller than pix.min (soil, shadows, plot numbering, etc)
	lab <- bwlabel(bifil)
	tab <- data.frame(table(lab)); tab$lab <- as.numeric(as.character(tab$lab))
	stab <- tab$lab[which(tab$Freq < pix.min)]

	labels <- lab
	labels[which(lab@.Data %in% stab, arr.ind=T)] <- 0

	# find and grab the color card
	dist <- X <- Y  <-rep(NA, max(labels))
	for(i in 1:max(labels)){
		ix <- which(labels==i, arr.ind=T)
		X[i]=mean(ix[,1])
		Y[i]=mean(ix[,2])

		dist[i] <- sqrt(X[i]^2 + Y[i]^2)
	}
	labels2<-labels

	if(colorcard=="bottomright"){
		labels2[which(labels!=which.max(dist), arr.ind=T)] <- 0
	}

	if(colorcard=="topleft"){
		labels2[which(labels!=which.min(dist), arr.ind=T)] <- 0
	}

	if(colorcard=="topright" | colorcard=="bottomleft"){
		X2 <- X-dim(labels)[1] +1
		#Y2 <- Y-dim(labels)[2]+1
		dist2 <- sqrt(X2^2 + Y^2)

		if(colorcard=="topright"){
			labels2[which(labels!=which.min(dist2), arr.ind=T)] <- 0
		}

		if(colorcard=="bottomleft"){
			labels2[which(labels!=which.max(dist2), arr.ind=T)] <- 0
		}
	}



	im3 <- im2c
	im3@.Data[,,1][which(labels2==0)] <- 0
	im3@.Data[,,2][which(labels2==0)] <- 0
	im3@.Data[,,3][which(labels2==0)] <- 0



	# blur <- medianFilter(im3, 5)
	# cluster the grayscale version
	# gr <- channel(im3, "gray")

	# get top 24 colors... cluster dendrogram of RGB, cut tree at 24 groups, get means?
	# border of CC is same as black square? text same as white square?

	dat <- cbind(cbind(im3[,,1][which(labels2!=0, arr.ind=T)],
			 im3[,,2][which(labels2!=0, arr.ind=T)]),
			 im3[,,3][which(labels2!=0, arr.ind=T)]
			)

	# ix <- which(apply(dat, 1, sum)==0)
	# dat2 <- dat[-ix, ]

	card <- matrix(c(116,81,67,199,147,129,91,122,156,90,108,64,130,128,176,92,190,172,224,124,47,68,91,170,198,82,97,94,58,106,159,189,63,230,162,39,34,63,147,67,149,74,180,49,57,238,198,32,193,84,151,12,136,170,243,238,243,200,202,202,161,162,161,120,121,120,82,83,83,49,48,51), nrow = 24, ncol = 3, byrow = T)
	card2 <- as.matrix(card/255)

	# coordinates of the center of each color chip on the card
	firstD<-c(0.20,0.41,0.63,0.85)
	secondD<-c(0.11,0.26,0.40,0.55,0.69,0.84)

	# crop image to only include color card
	im5<-im3[min(which(im3[,,]>0, arr.ind=T)[,1]):max(which(im3[,,]>0, arr.ind=T)[,1]),min(which(im3[,,]>0, arr.ind=T)[,2]):max(which(im3[,,]>0, arr.ind=T)[,2]),]

	# get color of center pixel on each color chip
	chip.masks=c(1:24)
	firstD.count=c(rep.int(4,6), rep.int(3,6),rep.int(2,6),rep.int(1,6))
	secondD.count=rep.int(1:6,4)
	obs.land<-matrix(data=NA, nrow=24, ncol=3)
	for(z in chip.masks){
	  obs.land[z,]<-rbind(im5@.Data[round(length(im5[,1,1])*firstD[firstD.count[z]]),round(length(im5[1,,1])*secondD[secondD.count[z]]),1],im5@.Data[round(length(im5[,1,1])*firstD[firstD.count[z]]),round(length(im5[1,,1])*secondD[secondD.count[z]]),2],im5@.Data[round(length(im5[,1,1])*firstD[firstD.count[z]]),round(length(im5[1,,1])*secondD[secondD.count[z]]),3])
	}

  return(obs.land)

}


