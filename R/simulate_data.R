simulate <- function(D,
                     n,
                     skip=NULL,
                     treat=c(6,6,6),
                     params,
                     phi.star=NULL,
                     Psi=NULL,
                     truth=FALSE,
                     debug=FALSE){

    if(debug)
        browser()

    ## Define observation times
    N <- D*n
    delta <- 24/n
    times <- seq(delta,24*D,delta)
    
    ## Activity start times
    if(!is.null(phi.star)){
        s <- rnorm(1,params$mu.s1,params$sigma.s)

        snext <- rnorm(1,s+24+phi.star(s),params$sigma.s)
        
        while(snext < max(times)){
            s <- c(s,snext)
            snext <- rnorm(1,snext+24+phi.star(snext),params$sigma.s)
        }

        J <- length(s)
    }
    else{
        phi <- Psi %*% params$c
        s <- params$mu.s1 + cumsum(rnorm(D,c(0,24 + phi),params$sigma.s))
    }

    ## Activity end times
    jtreat <- apply(outer(s,24*cumsum(treat),">"),1,sum) + 1
    e <- rnorm(J,s + params$mu.e[jtreat],params$sigma.e)

    ## Activity
    A1 <- outer(times,s,">")
    A2 <- outer(times,e+delta,"<=")

    if(!is.null(skip)){
        A1[,skip] <- 0
        A2[,skip] <- 0
    }
    
    A <- apply(A1*A2,1,sum)

    ## Activity counts

    ## 1) Segment specific parameters
    shape <- params$alpha
    rate <- params$alpha/params$mu

    ## 2) Time periods per treatment
    N1 <- treat*n

    lambda <- do.call("rbind",
                      lapply(1:length(treat),function(i){
                          cbind(rgamma(N1[i],shape=shape[i,1],rate=rate[i,1]),
                                rgamma(N1[i],shape=shape[i,2],rate=rate[i,2]))
                      }))
    
    Y <- rpois(N,(1-A)*lambda[,1] + A*lambda[,2])

    segments <- apply(outer(s,24*cumsum(treat),"<") *
                      outer(s,24*c(0,cumsum(treat[-length(treat)])),">"),
                      2,sum)

    segments[length(treat)] <- segments[length(treat)]-1
    
    
    if(truth)
        return(list(J=J,
                    treat=treat,
                    segments=segments,
                    s=s,
                    e=e,
                    sstar=delta*ceiling(s/delta),
                    estar=delta*ceiling(e/delta),
                    A=A,
                    lambda=lambda,
                    Y=Y))
    else
        return(Y)
}



    
