# Michael Miller, Shannon Lab, Jan. 2021
# *Tub*er *A*nalysis in *R*

##############
##COLOR DATA##
##############

#' Turn the output from skin.all into a matrix of median trait values
#'
#' @param color_list Output from skin.all function
#' @param remove_tag If TRUE, objects below `green_thresh` or above `red_thresh` will be removed from the data before taking median values
#' @param green_thresh Minimum CIELAB green-red value for object data to not be removed
#' @param red_thresh Maximum CIELAB green-red value for object data to not be removed
#'
#' @return A matrix of the median values of redness, skinning, and lightness for each image
#' @examples
#'   load(system.file("example_data", "skin_exmp.rds", package = "TubAR"))
#'   skin.export(skin_exmp, remove_tag=TRUE, green_thresh=-5, red_thresh=50)
#' @importFrom Biobase subListExtract
#' @export
skin.export<-function(color_list,remove_tag=TRUE, green_thresh=-5, red_thresh=50){
  # seperate data for each metric
  TBskin<-subListExtract(color_list$by.tuber,"skinning")
  TBred<-subListExtract(color_list$by.tuber,"redness")
  TBlight<-subListExtract(color_list$by.tuber,"lightness")

  # turn the data lists into arrays
  TBskint=array()
  for(item in TBskin){

    TBskint<-suppressWarnings(cbind(TBskint,as.array(item)))

  }
  TBskint<-TBskint[,2:length(colnames(TBskint))]
  colnames(TBskint)<-names(TBskin)

  TBredt=array()
  for(item in TBred){

    TBredt<-suppressWarnings(cbind(TBredt,as.array(item)))

  }

  TBredt<-TBredt[,2:length(colnames(TBredt))]
  colnames(TBredt)<-names(TBred)

  TBlightt=array()
  for(item in TBlight){

    TBlightt<-suppressWarnings(cbind(TBlightt,as.array(item)))

  }

  TBlightt<-TBlightt[,2:length(colnames(TBlightt))]
  colnames(TBlightt)<-names(TBlight)

  # remove red and green tag data and replace it with the median value of the line
  if(remove_tag==TRUE){
    redInd<-arrayInd(which(TBredt>red_thresh),dim(TBredt))
    greenInd<-arrayInd(which(TBredt<green_thresh),dim(TBredt))

    tag_replace_red<-vector(length=length(redInd[,2]))
    tag_replace_green<-vector(length=length(greenInd[,2]))

    skin_replace_red<-vector(length=length(redInd[,2]))
    skin_replace_green<-vector(length=length(greenInd[,2]))

    light_replace_red<-vector(length=length(redInd[,2]))
    light_replace_green<-vector(length=length(greenInd[,2]))
    if(length(redInd)!=0){
      c=0
      for(i in redInd[,2]){
        c<-c+1
        tag_replace_red[c]<-median(TBredt[,i])
        skin_replace_red[c]<-median(TBskint[,i])
        light_replace_red[c]<-median(TBlightt[,i])
      }
    }
    if(length(greenInd)!=0){
      c=0
      for(i in greenInd[,2]){
        c<-c+1
        tag_replace_green[c]<-median(TBredt[,i])
        skin_replace_green[c]<-median(TBskint[,i])
        light_replace_green[c]<-median(TBlightt[,i])
      }
    }
    if(length(redInd)!=0){
      for(i in c(1:length(tag_replace_red))){
        TBredt[redInd[i,1],redInd[i,2]]<-tag_replace_red[i]
        TBskint[redInd[i,1],redInd[i,2]]<-skin_replace_red[i]
        TBlightt[redInd[i,1],redInd[i,2]]<-light_replace_red[i]
      }
    }
    if(length(greenInd)!=0){
      for(i in c(1:length(tag_replace_green))){
        TBredt[greenInd[i,1],greenInd[i,2]]<-tag_replace_green[i]
        TBskint[greenInd[i,1],greenInd[i,2]]<-skin_replace_green[i]
        TBlightt[greenInd[i,1],greenInd[i,2]]<-light_replace_green[i]
      }
    }
  }

  # flip arrays so that lines are rows
  TBredt<-t(TBredt)
  TBskint<-t(TBskint)
  TBlightt<-t(TBlightt)

  # find medians for rows
  TBred_med<-apply(TBredt,1,median)
  TBskin_med<-apply(TBskint,1,median)
  TBlight_med<-apply(TBlightt,1,median)

  # Make export matrix
  color_matrix<-matrix(c(TBred_med,TBskin_med,TBlight_med),nrow=length(TBred_med), ncol=3)
  rownames(color_matrix)<-names(TBred_med)
  colnames(color_matrix)<-c("redness","skinning","lightness")

  return(color_matrix)
}


##############
##SHAPE DATA##
##############

#' Turn the output from shape.all into a matrix of median trait values
#'
#' @param shape_dat Output from shape.all function
#' @param remove_tag If TRUE, objects below the `roundness_thresh` will be removed from the data before taking median values
#' @param roundness_thresh The minimum roundness for an object to not be removed from the output, default is 0.7
#'
#' @return A matrix of the median values of BB width, BB length, perimeter, convex perimeter, area, convex hull area, compactness, roundness, and maximum length values for each image.
#' @examples
#'   load(system.file("example_data", "shape_exmp.rds", package = "TubAR"))
#'   shape.export(shape_exmp,remove_tag=TRUE, roundness_thresh=0.7)
#' @export
shape.export<-function(shape_dat, remove_tag=TRUE, roundness_thresh=0.7){
  if(remove_tag == TRUE){
    # Remove tags from shape data by removing tubers below a roundness threshold
    shape_dat=shape_dat[which(shape_dat[,8]>roundness_thresh),]#0.7 is used as a threshold that appeared to be higher than the roundness of any tags, this may need to be changed if using different tags/unusual potatoes
  }

  files <- list.files( pattern="*.jpg")
  shape_row_names <- gsub(".jpg", "", files)

  # create vectors to hold shape medians for each line
  bbox.width<-vector(length = length(shape_row_names))
  bbox.height<-vector(length = length(shape_row_names))
  perim<-vector(length = length(shape_row_names))
  convex.perim<-vector(length = length(shape_row_names))
  area<-vector(length = length(shape_row_names))
  chull.area<-vector(length = length(shape_row_names))
  roundness<-vector(length = length(shape_row_names))
  compactness<-vector(length = length(shape_row_names))
  max.length<-vector(length = length(shape_row_names))

  # turn the first row into the rownames
  rownames(shape_dat)<-shape_dat[,1]

  # find line medians and put medians into vectors
  count<-c(1:length(shape_row_names))
  for(i in count){
    bbox.width[i]<-median(as.numeric(shape_dat[startsWith(rownames(shape_dat),shape_row_names[i]),2]))
    bbox.height[i]<-median(as.numeric(shape_dat[startsWith(rownames(shape_dat),shape_row_names[i]),3]))
    perim[i]<-median(as.numeric(shape_dat[startsWith(rownames(shape_dat),shape_row_names[i]),4]))
    convex.perim[i]<-median(as.numeric(shape_dat[startsWith(rownames(shape_dat),shape_row_names[i]),5]))
    area[i]<-median(as.numeric(shape_dat[startsWith(rownames(shape_dat),shape_row_names[i]),6]))
    chull.area[i]<-median(as.numeric(shape_dat[startsWith(rownames(shape_dat),shape_row_names[i]),7]))
    roundness[i]<-median(as.numeric(shape_dat[startsWith(rownames(shape_dat),shape_row_names[i]),8]))
    compactness[i]<-median(as.numeric(shape_dat[startsWith(rownames(shape_dat),shape_row_names[i]),9]))
    max.length[i]<-median(as.numeric(shape_dat[startsWith(rownames(shape_dat),shape_row_names[i]),10]))
  }
  colheads<-matrix(c("bbox.width","bbox.height","perim","convex.perim","area","chull.area","roundness","compactness","max.length"),nrow=9,ncol=1)

  # Make export matrix
  shape_matrix<-matrix(c(bbox.width,bbox.height,perim,convex.perim,area,chull.area,roundness,compactness,max.length),nrow=length(count),ncol=9, byrow=FALSE)
  dimnames(shape_matrix)<-list(shape_row_names,colheads[1:9])

  return(shape_matrix)
}




