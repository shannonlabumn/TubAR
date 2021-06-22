####Failed color correction masking methods####
install.packages("Morpho")
library(Morpho)

tps3d()
computeTransform()
applyTransform()
help(Morpho)
Morpho:nose
data(nose)


data(nose)
## define some landmarks
refind <- c(1:3,4,19:20)
## use a subset of shortnose.lm as anchor points for a TPS-deformation
reflm <- shortnose.lm[refind,]
tarlm <- reflm
##replace the landmark at the tip of the nose with that of longnose.lm
tarlm[4,] <- longnose.lm[4,]
## deform a set of semilandmarks by applying a TPS-deformation
## based on 5 reference points
deformed <- tps3d(shortnose.lm, reflm, tarlm,threads=1)
## Not run:
##visualize results by applying a deformation grid
deformGrid3d(shortnose.lm,deformed,ngrid = 5)

obs.land2<-matrix(c(im4@.Data[144,75,],im4@.Data[144,79,],im4@.Data[144,83,],im4@.Data[144,86,],im4@.Data[144,91,], im4@.Data[144,95,],im4@.Data[140,75,], im4@.Data[140,79,], im4@.Data[140,83,],im4@.Data[140,87,], im4@.Data[140,91,],im4@.Data[140,95,],im4@.Data[136,75,], im4@.Data[136,79,], im4@.Data[136,83,],im4@.Data[136,87,], im4@.Data[136,91,],im4@.Data[136,95,],im4@.Data[132,75,], im4@.Data[132,80,], im4@.Data[132,83,],im4@.Data[132,87,], im4@.Data[132,91,],im4@.Data[131,95,]), ncol=3,byrow=TRUE)
obs.land<-grabcard(image)
colnames(obs.land)<-colnames(exp.land)
obs.land
exp.land<-card2
exp.land
rgb[1]
display(im4)
rgb.obs<-as.matrix(rgb[1])[[1]]
colnames(rgb.obs)<-colnames(exp.land)

255*tps3d(rgb.obs, obs.land3, exp.land)

image<-images[1]
obs.land<-grabcard(image)
obs.land<-newcenters
obs.land <-255*obs.land

obs.land<-obs.land+obs.diff
obs.land[which(obs.land<0)]<-0
obs.land[which(obs.land>255)]<-255

obs.colors<-rgb(round(obs.land[,1]), green=round(obs.land[,2]), blue=round(obs.land[,3]), maxColorValue=255)
barplot(rep(1, 24), axes = FALSE, space = 0, col = obs.colors)


obs.land2<-obs.land2*255
obs.colors2<-rgb(round(obs.land2[,1]), green=round(obs.land2[,2]), blue=round(obs.land2[,3]), maxColorValue=255)
barplot(rep(1, 24), axes = FALSE, space = 0, col = obs.colors2)


exp.land<-card2
exp.land <-255*exp.land
exp.colors<-rgb(round(exp.land[,1]), green=round(exp.land[,2]), blue=round(exp.land[,3]), maxColorValue=255)
barplot(rep(1, 24), axes = FALSE, space = 0, col =exp.colors)


obs.diff<-obs.land2-obs.land
obs.diff
obs.land-exp.land


cor.pota<-tps3d(rgb.obs, obs.land2, exp.land)[200:250,]
cor.pota<-255*cor.pota
pota.colors<-rgb(round(cor.pota[,1]), green=round(cor.pota[,2]), blue=round(cor.pota[,3]), maxColorValue=255)
barplot(rep(1, 50), axes = FALSE, space = 0, col =pota.colors)


display(im3)
im5<-resize(im3, 100)
im5@.Data<-im3@.Data[which(im3@.Data >0, arr.ind=TRUE)]

im3r<-resize(im3, 500)
im5dims<-as.vector(apply(which(im3r@.Data >0, arr.ind=TRUE),2,max)-apply(which(im3r@.Data >0, arr.ind=TRUE),2,min)+1)
im5dims[2]<-87
im5@.Data<-array(im3r@.Data[which(im3r@.Data >0)],dim = im5dims)
display(im5)
display(im3r)
obs.crop<-matrix(c(im5@.Data[17,4,],im5@.Data[17,8,],im5@.Data[17,12,],im5@.Data[17,16,],im5@.Data[17,20,], im5@.Data[17,24,],im5@.Data[12,4,], im5@.Data[12,8,], im5@.Data[12,12,],im5@.Data[12,16,], im5@.Data[12,20,],im5@.Data[12,24,],im5@.Data[9,4,], im5@.Data[9,8,], im5@.Data[9,12,],im5@.Data[9,16,], im5@.Data[9,20,],im5@.Data[9,24,],im5@.Data[5,4,], im5@.Data[5,8,], im5@.Data[5,12,],im5@.Data[5,16,], im5@.Data[5,20,],im5@.Data[5,24,]), ncol=3,byrow=TRUE)
obs.crop<-obs.crop*255
obs.crop2<-rgb(round(obs.crop[,1]), green=round(obs.crop[,2]), blue=round(obs.crop[,3]), maxColorValue=255)
barplot(rep(1, 24), axes = FALSE, space = 0, col = obs.crop2)

browseVignettes(package= "EBImage")

im5<-im3
display(translate(im5, c(7,11)))
display(im5)

im5dims<-as.vector(apply(which(im5 !=0, arr.ind=TRUE),2,max)-apply(which(im5 !=0, arr.ind=TRUE),2,min)+1)
ind0<-which(im5 != 0, arr.ind=TRUE)
ind1r<-ind0[,1][which(ind0[,3]==1)]
ind2r<-ind0[,2][which(ind0[,3]==1)]
ind1g<-ind0[,1][which(ind0[,3]==2)]
ind2g<-ind0[,2][which(ind0[,3]==2)]
ind1b<-ind0[,1][which(ind0[,3]==3)]
ind2b<-ind0[,2][which(ind0[,3]==3)]





i=0
redvec<-vector(length=length(ind1r))
grnvec<-vector(length=length(ind1r))
bluvec<-vector(length=length(ind1r))
while(i<5272){
  i=i+1
  redvec[i]<-im5[ind1r[i],ind2r[i],1]
  grnvec[i]<-im5[ind1g[i],ind2g[i],2]
  bluvec[i]<-im5[ind1b[i],ind2b[i],3]
}
print_factors(44250)
im5dims
allvec<-rgb(redvec, green=grnvec, blue=bluvec)
177*250
print_factors(15816)
barplot(rep(1, 204), axes = FALSE, space = 0, col = unique(allvec)[colsort$ix[1:204]])
allmat<-matrix(allvec,,)
display(Image(allmat))
max(allvec)
length(unique(allvec))

unique(allvec)%in% allvec
colsort<-sort(tabulate(match(allvec,unique(allvec))), decreasing = T, index.return=T)

unique(allvec)[colsort$ix[1:50]]

data.frame(table(bwlabel(im5)))
mask()
display(colorLabels(propagate(seeds = im5, x = im5, lambda = 100)))
display(thresh(im5,w=200,h=100))


unique(allvec)[1139]

filter2(im5,makeBrush(9)/sum(makeBrush(9)))
im5t<-thresh(im5, w=10, h=10, offset=0.01)
im5t<-opening(im5t, makeBrush(5, shape='box'))
im5t<-fillHull(im5t)
im5t<-bwlabel(im5t)
display(im5t)

table(im5t)

max(im5t)

im5lab<-bwlabel(im5)
im5[which(im5lab==1, arr.ind = TRUE)]

display(translate(im5, c(30,30)))
image=translate(im5, c(30,30))


w = makeBrush(size = 31, shape = 'gaussian', sigma = 4)
im3<-filter2(im3, w)

display(im3)

nmask = thresh(im3, w=5, h=5, offset=0.05)
nmask = opening(nmask, makeBrush(5, shape='box'))
nmask = fillHull(nmask)
nmask = bwlabel(nmask)

display(nmask)

im3@.Data[which(im3@.Data<0.5, arr.ind=T)]<- 0
display(im3[1283:1455,721:982,])
display(which(im3[,,]>0, arr.ind=T))

orderX[1:6]
Kxy$centers
Kxy$centers[25,2]
image<-"FY3Rd19_9_26.jpg"
grabcard("2020fy2_317.jpg")
# my plan was to use the warped spline method to correct image colors
# see matlab script in Appendix to Menesatti et al 2012 Sensors 12:7063-7079
orderX<-orderX[which(order(Kxy$centers, decreasing=TRUE)<25)]
order(Kxy$centers, decreasing=TRUE)
image<-"FY3Ru12_8_7.jpg"
image<-"FY3Ru10_2_7.jpg"
w = makeBrush(size = 31, shape = 'gaussian', sigma = 5)
im3<-filter2(im3, w)

####kmeans for color correction####

# 3D k-means clustering, give it the color card points as centroids
# Khw <- kmeans(dat, centers=card[, c(6:8)]/100, algorithm="Hartigan-Wong", iter.max=100)
# Kl <- kmeans(dat, centers=card[, c(6:8)]/100, algorithm="Lloyd", iter.max=100) # == Forgy
# Km <- kmeans(dat, centers=card[, c(6:8)]/100, algorithm="MacQueen", iter.max=100)

# 5D (XY RGB) k-means -- much closer but not quite right due to edges
ixy <- which(labels2!=0, arr.ind=T)
datxy <- cbind(ixy, dat)

Kxy <- kmeans(datxy, centers=24, iter.max=50)


# explore alternatives to kmeans? doesn't seem necessary

# assign new centers to color card numbering
orderX <- order(Kxy$centers, decreasing=TRUE)

col4 <- names(sort(Kxy$centers[orderX[1:6],2], decreasing=F))
col3 <- names(sort(Kxy$centers[orderX[7:12],2], decreasing=F))
col2 <- names(sort(Kxy$centers[orderX[13:18],2], decreasing=F))
col1 <- names(sort(Kxy$centers[orderX[19:24],2], decreasing=F))

cardorder <- as.numeric(c(col4,col3,col2,col1))





