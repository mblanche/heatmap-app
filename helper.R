library(edgeR)
library(RColorBrewer)

## Reading the DGE data from disk
saved.data <- readRDS("data/dge.rds")
gene2name <- readRDS("data/gene2name.rds")
hc.cols <- brewer.pal(9,"Set1")

## This function takes a matrix of fold change
## with the column names as sample names
## It also takes a gene wise cluster to sort the matrix and draw a dendrogram
## Finaly, it takes a zlim (a min and max range) to adjust the heatmap contrast
plotHeatMap <- function(d,hc,zlim,height,noMarker=FALSE){
    ## Picking up a color scheme
    ## because the zlim have different value, I need to center them on black
    ## So, i'll have a high and a low palette
    cols <- colorRampPalette(c("blue","black","yellow"))(256)
    
    ## Let's try to map the continous value in the matrix
    ## to the discrete color in cols, keeping 0 as black (128)
    ## Rescale the d to set the range to zlim
    d[d < zlim[1]] <- zlim[1]
    d[d > zlim[2]] <- zlim[2]
    t <- d
    t[d<=0] <- c(1:128)[cut(d[d <= 0],128,labels=FALSE)]
    t[d>0] <- c(127:256)[cut(d[d > 0],127,labels=FALSE)]
    d <- t

    ## Creating the canvas
    layout(matrix(c(1,2,0,3),ncol=2,byrow=TRUE)
           ,widths=c(30,70),heights=c(85,15))
    ## Plotting the gene dendrogram
    par(mar = c(2,1,2,0))
    plot(hc,
         axes=TRUE,
         yaxt='s',
         yaxs='i',
         xaxt='n',
         xaxs='i',
         horiz=TRUE,
         leaflab='none',
         ylab="Height")
    ## Add a marker where dendrogram is broken
    if(noMarker == FALSE) abline(v=height,lty=2,col='red')
    ## Plotting the heatmap
    par(mar = c(2,0.2,2,1))
    image(t(d[unlist(hc),]),
          col=cols,
          zlim=c(1,256),
          xaxt='n',
          yaxt='n'
          )
    ## Ploting the sample names
    if(ncol(d) == 1){
        at <- 1
    } else {
        at <- seq(0,1,1/(ncol(d)-1))
    }
    axis(side=1,at=at,labels=colnames(d))
    par(mar = c(3,0.2,2,1))
    ## Plotting the color scale
    image(matrix(seq(zlim[1],zlim[2],length=100)),
          col=cols,
          xaxt='n',
          yaxt='n')
    axis(side=1,at=c(0,0.5,1),labels=format(c(zlim[1],0,zlim[2]),digits=2))

}

colbranches <- function(n, col) {
    a <- attributes(n) 
    attr(n, "edgePar") <- c(a$edgePar, list(col=col,lty=1,lwd=1))
    return(n) # Don't forget to return the node!
}


colBranches <- function(hc,height,cols=cols){
    .col.counter <<- 0
    breakAt <- function(hc,height,cols){
        res <- lapply(hc,function(sub.hc){
            if (attributes(sub.hc)$height < height){
                if (.col.counter == length(cols))
                    .col.counter <<- 0
                .col.counter <<- .col.counter+1
                return(dendrapply(sub.hc,colbranches,cols[.col.counter]))
            } else {
                res <- breakAt(sub.hc,height,cols)
                sub.hc <- do.call(merge,c(res,height=attributes(sub.hc)$height))
                attr(sub.hc,'edgePar') <- list(lty=2,col='black')
                return(sub.hc)
            }
        })
    }
    res <- breakAt(hc,height,cols)
    do.call(merge,c(res,height=attributes(hc)$height))
}
