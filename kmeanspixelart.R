# Restaurando pixelart con R y k-means
# www.datosimagensonido.com

library(tiff)

# RGB 3D-plot of matrix M
PlotScatterRGB=function(M, sub=NULL, radius=0.02) {
    library(rgl)
    open3d()
    for (i in 1:nrow(M)) {
        spheres3d(M[i,1], M[i,2], M[i,3],
                  color=rgb(M[i,1], M[i,2], M[i,3]), radius=radius)
    }
    bg3d(color="lightgray")
    axes3d()
    title3d(sub=sub, xlab='R', ylab='G', zlab='B')
}


# INPUT IMAGES PARAMETERS
pixelart=c("bobafett", "pokemon", "octocat", "sonic")  # TIFF filename
dimx=c(38, 28, 45, 23)  # width
dimy=c(43, 24, 41, 33)  # height
SAFE=c(2, 10, 1, 2)  # border pixels to drop
K=c(15, 4, 8, 14)  # k-means k


# LOOP THROUGH PIXELART
for (n in 1:length(pixelart)) {

    img=readTIFF(paste0(pixelart[n], ".tif"), native=F, convert=F)
    
    # MEDIAN AVERAGING
    DIMX=dimx[n]
    DIMY=dimy[n]
    print(paste0("Rescaling '", pixelart[n],
        "' to ", DIMX, "x", DIMY, " pixels..."))
    img_median=array(0, c(DIMY,DIMX,3))
    
    NX=(dim(img)[2]/DIMX)
    NY=(dim(img)[1]/DIMY)
    img_used=img*0
    for (i in 1:DIMY) {
        for (j in 1:DIMX) {
            for (k in 1:3) {
                xini=round((j-1)*NX+1+SAFE[n])
                xfin=round(j*NX-SAFE[n])
                yini=round((i-1)*NY+1+SAFE[n])
                yfin=round(i*NY-SAFE[n])
                
                img_median[i,j,k]=median(img[yini:yfin, xini:xfin, k])
                img_used[yini:yfin, xini:xfin, k]=img[yini:yfin, xini:xfin, k]
            }
        }
    }
    
    writeTIFF(img_median, paste0(pixelart[n],"_median_",DIMX,"x",DIMY,".tif"),
        bits.per.sample=8, compression="LZW")
    writeTIFF(img_used, paste0(pixelart[n],"_used.tif"),
        bits.per.sample=8, compression="LZW")
    
    
    # 3D-PLOT OF RGB VALUES
    # Visual cluster estimation
    R=img_median[,,1]
    G=img_median[,,2]
    B=img_median[,,3]
    dim(R)=length(R)
    dim(G)=length(G)
    dim(B)=length(B)
    
    M=array(0,c(length(R),3))  # Matriz RGB
    colnames(M)=c("R","G","B")
    M[,1]=R
    M[,2]=G
    M[,3]=B
    PlotScatterRGB(M, sub=paste0("'",pixelart[n],"' median values"),
        radius=0.015)
    
    
    # K-MEANS CLUSTERING
    NCOLOURS=K[n]  # k clusters
    set.seed(0)  # reproducible segmentation
    kmeansfit=kmeans(subset(M, select=c("R","G","B")), centers=NCOLOURS,
        nstart=90000, iter.max=100000)   # high nstart can prevent from
    clustering=kmeansfit$cluster         # missing the tiniest clusters
    
    # Clustering histogram
    png(paste0(pixelart[n],"_histogram.png"), width=512, height=400)
    breaks=seq(0, NCOLOURS, length.out=NCOLOURS+1)
    hist(clustering, breaks=breaks, col='gray',
         main=paste0("'",pixelart[n],"' cluster histogram (k=",
        K[n], ")"), axes=F)
    axis(1, at=breaks, labels=T)
    dev.off()
    
    imgcentersR=array(0, DIMY*DIMX)
    imgcentersG=imgcentersR
    imgcentersB=imgcentersR
    centers=kmeansfit$centers
    for (i in 1:NCOLOURS) {
        indices=which(clustering==i)
        imgcentersR[indices]=centers[i,1]
        imgcentersG[indices]=centers[i,2]
        imgcentersB[indices]=centers[i,3]
    }
    dim(imgcentersR)=c(DIMY,DIMX)
    dim(imgcentersG)=c(DIMY,DIMX)
    dim(imgcentersB)=c(DIMY,DIMX)
    
    img_restored=array(0, c(DIMY,DIMX,3))
    img_restored[,,1]=imgcentersR
    img_restored[,,2]=imgcentersG
    img_restored[,,3]=imgcentersB
    writeTIFF(img_restored, paste0(pixelart[n],"_restored_",
        DIMX,"x",DIMY,"_",NCOLOURS,".tif"),
        bits.per.sample=8, compression="LZW")
    
    # 3D-plot of centroids
    PlotScatterRGB(centers,
        sub=paste0("'",pixelart[n],"' centroids (k=",K[n],")"))

}
