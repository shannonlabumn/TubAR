# Cari Schmitz Carley, Shannon Lab, Jan. 2019
# *Tub*er *A*nalysis in *R*

#' Perform `find.skin` for all JPEG files in the working directory.
#'
#' @param n.core The number of processor cores to use in processing. Default is 1.
#' @param display If TRUE, a black and white image of each object will be created for each image
#' @param mode If set to “debug” the display will show numbered objects and skinned area on the objects.
#' @param write.clean If TRUE, a JPEG image with backgroung removed will be created for each input image
#' @param pix.min The minimum pixel size required for objects to not be removed with the background,
#' @param scaledown The amount image sizes will be reduced in order to aid processing speed
#' @param colorcard Declares which corner the color card is in. "bottomright" by default. Can also be set to "bottomleft", "topright", and "topleft"
#' @param color.correct Determines if pixel values will be color corrected based on the color card. TRUE by default.
#' @param background.color The color of the background of the image. Black, blue, and white are currently supported. Default is white.
#' @param color.values Color correction will change the color card RGB values to match this value. "Default" uses RGB values from the manufacturer CameraTrax. "Revised" uses values that result in more "natural" colors. Otherwise, this should be a matrix of the RGB values for the 24-color card starting with the top right color going down, then left.
#' @param color.center Two item list of the coordinates of the center pixel of each color chip in terms of the proportion of the card to the left/above the pixel. The list should contain 2 vectors, the first going left to right and the second going up to down.
#'
#' @return A nested list of median values for each image of rednees, skinning, and lightness and the values forthese traits for every object in every image.
#' @examples
#'   setwd(system.file("images", package = "TubAR"))
#'   skin.all(n.core=1, display=T, mode="debug", write.clean=F, pix.min=4e3, scaledown=4, colorcard="bottomright")
#' @export

# wrapper -- before running, the wd needs to be set to the directory containing the images
skin.all <- function(n.core=1, display=T, mode="debug", write.clean=F, pix.min=4e3, scaledown=4, colorcard="bottomright", color.correct=T, background.color="white",  color.values= "revised", color.center = "default"){
  files <- list.files( pattern=".*\\.(jpg|JPG)$")
  names <- gsub(".jpg", "", files, ignore.case = T)

	results <- lapply(files, function(x) find.skin(x, display=display, mode=mode, write.clean=write.clean, pix.min=pix.min, scaledown=scaledown, colorcard=colorcard, n.core=n.core, color.correct=color.correct, background.color=background.color, color.values= color.values, color.center = color.center))

	names(results) <- names

	table <- matrix(NA, length(names), 3); rownames(table) <- names; colnames(table)<- c("%skinned", "redness", "lightness")
	for(t in 1:dim(table)[1]){
		table[t, ] <- unlist(lapply(results[[t]], function(x) median(x)))
	}

	return(list(median.values=table , by.tuber =results))
}



#' Find skinning threshold across the image and return skinning percent by tuber
#' assumes all potatoes in the image are a sampling of the same clone ( will have same difference between skin/flesh colors)
#'
#' @param image The image to be processed, must be a JPEG
#' @param display =T plots to quartz
#' @param mode ="debug" plots skinned area, numbers each tuber
#' @param write.clean =T saves a 'clean_image.jpg' with background whited out (for presentations, etc)
#' @param pix.min is number of pixels expected for the smallest object to be designated a tuber
#' @param scaledown Multiple by which image is divided for faster computing (a reduction in image size)
#' @param colorcard = NULL/"bottomright" to remove a color card if used
#' @param n.core The number of processor cores to use in processing. Default is 1.
#' @param color.correct Determines if pixel values will be color corrected based on the color card. TRUE by default.
#' @param background.color The color of the background of the image. "black", "blue", and "white" are currently supported. Default is white.
#' @param color.values Color correction will assume the color card is the one described in the Vignette, noted by "default". Otherwise, this should be a matrix of the RGB values for the 24-color card starting with the top right color going down, then left.
#' @param color.center Two item list of the coordinates of the center pixel of each color chip in terms of the proportion of the card to the left/above the pixel. The list should contain 2 vectors, the first going left to right and the second going up to down.
#'
#' @return A list containing the redness, skinning, and lightness of each object in the image
#' @examples
#'   find.skin(system.file("images", "2020fy2_317.jpg", package = "TubAR"))
#'   find.skin(system.file("images", "2020fy2_317.jpg", package = "TubAR"),display=T, mode="debug", write.clean=F, pix.min=4e4, scaledown=7, colorcard="bottomright", n.core=2, color.correct=T)
#' @import EBImage
#' @import future
#' @import future.apply
#' @import boot
#' @import minpack.lm
#' @export

find.skin <- function(image, display=T, mode="debug", write.clean=F, pix.min=4e4, scaledown=8, colorcard="bottomright", n.core=1, color.correct=T, background.color="white",  color.values= "revised", color.center = "default"){

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

	# RGB values for each tuber
	rgb <- sapply(p, function(x) cbind(im2@.Data[,,1][which(labels2==x)],
									   im2@.Data[,,2][which(labels2==x)],
									   im2@.Data[,,3][which(labels2==x)]))

	# Color correction
	if(color.correct==T){
	  obs.land<-grabcard(image, colorcard=colorcard, scaledown=scaledown, pix.min=pix.min, color.center = color.center)
	  if (length(obs.land)==72){
	    if(color.values == "default"){
	      #Matrix of each chip 255 scale RGB starting in top right corner
	      #RGB values given by color card manufacturer
	      card <- matrix(c(116,81,67,199,147,129,91,122,156,90,108,64,130,128,176,92,190,172,224,124,47,68,91, +
	                         170,198,82,97,94,58,106,159,189,63,230,162,39,34,63,147,67,149,74,180,49,57,238,198, +
	                         32,193,84,151,12,136,170,243,238,243,200,202,202,161,162,161,120,121,120,82,83,83,49,48,51), nrow = 24, ncol = 3, byrow = T)
	      card2 <- as.matrix(card/255)
	    }else if(color.values == "revised"){
	      #RGB of FY2_126.jpg 2019 FY2 sample. Skinning function was developed around similar lighting.
	      card2 <- matrix(c(0.4823529,0.4156863,0.4784314,0.9019608,0.7058824,0.7490196,0.5058824,0.6627451, +
	                          0.9058824,0.3490196,0.5960784,0.5137255,0.6823529,0.7019608,0.9607843,0.6431373, +
	                          0.9372549,0.9568627,0.9254902,0.5803922,0.4000000,0.3960784,0.4901961,0.9058824, +
	                          0.8352941,0.3803922,0.5568627,0.3607843,0.2666667,0.6117647,0.7254902,0.9176471, +
	                          0.5607843,0.9647059,0.7647059,0.4470588,0.2549020,0.3058824,0.7058824,0.3411765, +
	                          0.7764706,0.6078431,0.7098039,0.2313725,0.3215686,0.9803922,0.9137255,0.4666667, +
	                          0.8705882,0.4588235,0.8549020,0.4235294,0.7215686,0.9568627,0.9764706,0.9529412, +
	                          0.9607843,0.8862745,0.9137255,0.9372549,0.7568627,0.8196078,0.8823529,0.5490196, +
	                          0.6627451,0.7960784,0.3333333,0.4509804,0.5843137,0.2784314,0.3333333,0.4470588), nrow = 24, ncol = 3, byrow = T)
	    }else{
	      card <- color.values
	      card2 <- as.matrix(card/255)
	    }
	    rgb<-lapply(rgb, function(x)tps3d(x,obs.land,card2))
	  }else{
	    warning("Error: color correction failure. Is card crooked?")
	  }
	}

	# to Lab for just the tuber pixels
	if(n.core >1){
		# go parallel with future and furture.apply packages
  	  plan(multiprocess, workers=n.core)

 	   	Lab <- future_lapply(rgb, function(x) t(apply(x, 1, function(y) convertColor(y, from="sRGB", to="Lab"))))
  } else {
 	   	Lab <- lapply(rgb, function(x) t(apply(x, 1, function(y) convertColor(y, from="sRGB", to="Lab"))))

    	}

	thresh <- vector("list", length(p)); names(thresh)<-p # thresholds along the range of blue to yellow (b)
	skinper <- vector("list", length(p)); names(thresh)<-p # the percent skinning along the values of b

	for(i in 1:length(p)){
		thresh[[i]] <- c(min(Lab[[i]][,3], na.rm=T):max(Lab[[i]][,3], na.rm=T)) # along the range of b values
		skinper[[i]] <- sapply(thresh[[i]], function(x) sum(Lab[[i]][, 3] >x)/ (dim(Lab[[i]])[1]) )
	}

	# fit a sigmoid curve across the potatoes
	curve <- data.frame(x=unlist(thresh), y=unlist(skinper))
	sigm <- nlsLM(y~ a/(1+exp(-b * (x-c)))+d, data=curve, start=list(a=max(curve$y), b=1, c=median(curve$x), d=min(curve$y)))
	c=summary(sigm)$parameters["c", "Estimate"]
	thr <- c*1.5
	#thr <- c*0.78

	# get skinning percent at the threshold for each potato
	skinpot <- sapply(1:length(p), function(i) sum(Lab[[i]][,3] > thr)/( dim(Lab[[i]])[1]) )
	#skinpot <- sapply(1:length(p), function(i) sum(Lab[[i]][,3] < thr)/( dim(Lab[[i]])[1]) )

	# get median red color intensity within potato, outside of skinned area
	red.intens <- sapply(1:length(p), function(i) median(Lab[[i]][,2][which(Lab[[i]][,3] < thr)]))
	#red.intens <- sapply(1:length(p), function(i) median(Lab[[i]][,2][which(Lab[[i]][,3] > thr)]))

	# get median lightness within potato, outside of the skinned area
	lightness <- sapply(1:length(p), function(i) median(Lab[[i]][, 1][which(Lab[[i]][,3] < thr)]))
	#lightness <- sapply(1:length(p), function(i) median(Lab[[i]][, 1][which(Lab[[i]][,3] > thr)]))

	# VISUALIZATION
		labels3 <- labels2
		labels3[which(labels3 >0, arr.ind=T)] <- 1
		template <- abs (labels3 -1)


	if(display==T){

		for(i in 1:length(p)){
		  template[which(labels2==p[i])[which(Lab[[i]][,3]>thr)]] <- .5
			#template[which(labels2==p[i])[which(Lab[[i]][,3]<thr)]] <- .5
		}

		display(template, "raster")
    if(mode=="debug"){
		  L <- setdiff(unique(c(labels2)),0)
		 	for(i in L){
		 		 ix <- which(labels==i, arr.ind=T)
				 text(x=mean(ix[,1]), y=mean(ix[,2]),labels=which(L %in% i), cex=2, col="lightblue")

		 	}
    }

	}

	if(write.clean==T){
		test <- resize(template, h=dim(im2)[2]*scaledown)
		clean.im <- im

		R <- clean.im@.Data[,,1]
		R[which(test@.Data==1, arr.ind=T)] <- 1

		G <- clean.im@.Data[,,2]
		G[which(test@.Data==1, arr.ind=T)] <- 1

		B <- clean.im@.Data[,,3]
		B[which(test@.Data==1, arr.ind=T)] <- 1

		clean.im@.Data[,,1] <- R
		clean.im@.Data[,,2] <- G
		clean.im@.Data[,,3] <- B

		writeImage(clean.im, paste0("clean_", image))

	}


	result <- list(skinning=round(skinpot,2), redness=round(red.intens, 1), lightness=round(lightness, 1))
	return(result)

}





