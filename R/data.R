##' Build data object.
##'
##' Build data object.
##' @title Build Data Object
##' @param Y Activity series (length D*n)
##' @param D Number of days in study
##' @param n Observations per day
##' @param J Number of activity periods
##' @param offset Offset of initial time
##' @param treat Length of treatmens in days
##' @param segments Number of activity periods in each spline segment
##' @return Data object
##' @author Simon Bonner
##' @export
##' @importFrom splines spline.des
buildData <- function(Y,D,n,J=NULL,
                      treat=NULL,
                      segments=NULL,
                      skip=NULL,
                      model="Model1",
                      offset=0){

    ## Create basic data list
    data <- list(Y=Y,
                 D=D,
                 J=ifelse(is.null(J),D,J),
                 n=n,            # Observations per day
                 N=D*n,               # Total number of observations
                 delta=24/n,          # Time between observations (hours)
                 times=offset + seq(24/n,24*D,24/n), # Observation times (hours)
                 skip=skip)

    ## Add information on treatments
    if(is.null(treat)){
        data$treat <- D
        data$ntreat <- 1
        data$jtreat <- rep(1,data$J)
    }
    else{
        data$treat <- treat
        data$ntreat <- length(treat)
        data$jtreat <- c(rep(1:data$ntreat,segments),data$ntreat)
    }

    ## Add information on spline segments
    if(is.null(segments)){
        ## Use a single spline
        data$nsegments <- 1
        data$segments <- data$J-1
        data$segstart <- 1
        data$segend <- data$J-1
    }
    else{
        ## Spline is divided into multiple segments
        data$nsegments <- length(segments)
        data$segments <- segments
        data$segstart <- 1 + cumsum(c(0,segments[-data$nsegments]))
        data$segend <- cumsum(segments)
    }

    ## Construct spline design matrices
    Psi <- lapply(data$segments,function(j){
        genSpline(j)
    })

    data$Psi <- lapply(1:data$nsegments,function(k){
        Psi[[k]]$X
    })

    ## Store variance matrices for constructing priors
    data$V0 <- lapply(1:data$nsegments,function(k){
        c(1,rep(0,length(Psi[[k]]$V)))
    })

    data$V1 <- lapply(1:data$nsegments,function(k){
        c(0,Psi[[k]]$V)
    })

    ## Consolidate segments for Model 2
    if(model=="Model2"){
        dims <- t(sapply(data$Psi,dim))

        ## Consolidate Psi
        Psi <- cbind(data$Psi[[1]],
                     matrix(0,dims[1,1],sum(dims[-1,2])))

        if(length(data$Psi)>2){
            for(k in 2:(length(data$Psi)-1)){
                Psi <- rbind(Psi,
                             cbind(matrix(0,dims[k,1],sum(dims[1:(k-1),2])),
                                   data$Psi[[k]],
                                   matrix(0,dims[k,1],sum(dims[-(1:k),2]))))
            }
        }

        if(length(data$Psi)>1){
            k <- length(data$Psi)
            Psi <- rbind(Psi,
                         cbind(matrix(0,dims[k,1],sum(dims[1:(k-1),2])),
                               data$Psi[[k]]))
        }

        data$Psi <- list(Psi)

        ## Redefine segment data
        data$segments2 <- data$segments

        data$nsegments <- 1
        data$segments <- sum(data$segments)
        data$segstart <- 1
        data$segend <- data$segments

        ## Redefine prior parameters
        data$V0 <- list(unlist(data$V0))
        data$V1 <- list(t(bdiag(data$V1)))
    }

    return(data)
}
