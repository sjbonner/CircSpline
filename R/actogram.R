actogram <- function(data,range=NULL,double=FALSE,main=NULL){

    if(double){
        xlim <- c(-1,49)
        xaxis.at <- seq(0,48,6)
        xaxis.lab <- c(seq(0,24,6),seq(6,24,6))
    }
    else{
        xlim <- c(-1,25)
        xaxis.lab <- xaxis.at <- seq(0,48,3)
    }
    
    if(is.null(range))
        range <- c(1,data$D)

    K <- range[2]-range[1]+1
    
    plot(NA,NA,xlim=xlim,ylim=c(0,K*max(data$Y,na.rm=TRUE)),axes=FALSE,
         xlab="Time of Day (Hours)",ylab="",main=main)
    axis(1,at=xaxis.at+24/(2*data$n),label=xaxis.lab); box()

    for(d in range[1]:range[2]){
        base <- (range[2]-d)*max(data$Y,na.rm=TRUE)
        #abline(h=base)
        lines(range(xaxis.at),rep(base,2))

        if(double)
            tmp <- intersect(which(data$times > 24*(d-1)),
                             which(data$times <= 24*(d+1)))
        else
            tmp <- intersect(which(data$times > 24*(d-1)),
                             which(data$times <= 24*(d+1)))

        tmp1 <- tmp[which(data$Y[tmp]>0)]
        
        for(k in tmp1){
            lines(rep(data$times[k]-24*(d-1),2),
                  base + c(0,data$Y[k]))
        }
        text(xlim[1],base,d,cex=.5,adj=c(0.5,0))
    }

    abline(v=c(0,24,48)+24/(2*data$n),col="black")
    
    polygon(x=c(24,24,48,48,24)+24/(2*data$n),
            y=c(0,max(data$Y,na.rm=TRUE),max(data$Y,na.rm=TRUE),0,0),col="grey")

    text(12+24/(2*data$n),(range[2]+.5)*max(data$Y,na.rm=TRUE),"Day i",cex=.7)
    
    if(double)
        text(36+24/(2*data$n),(range[2]+.5)*max(data$Y,na.rm=TRUE),"Day i+1",cex=.7)
}

