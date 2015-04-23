buildDesign <- function(D,n,J=NULL,offset=0,treat=NULL,segments=NULL){
    design <- list(D=D,
                   J=ifelse(is.null(J),D,J),
                   treat=treat,
                   n=n,            # Observations per day
                   N=D*n,               # Total number of observations
                   delta=24/n,          # Time between observations (hours)
                   times=offset + seq(24/n,24*D,24/n)) # Observation times (hours)

    if(is.null(segments)){
        ## Use a single spline
        design$nsegments <- 1
        design$segments <- design$J-1
        design$segstart <- 1
        design$segend <- design$J-1
    }
    else{
        ## Spline is divided into multiple segments
        design$nsegments <- length(segments)
        design$segments <- segments
        design$segstart <- 1 + cumsum(c(0,segments[-design$nsegments]))
        design$segend <- cumsum(segments)
    }
    
    return(design)
}
