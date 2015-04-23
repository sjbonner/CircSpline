parameterList <-
    setRefClass("parameterList",
                fields=list(lambda="matrix",
                    mu.s1="numeric",
                    mu="matrix",
                    alpha="matrix",
                    s="numeric",
                    e="numeric",
                    A="numeric",
                    mu.e="numeric",
                    sigma.s="numeric",
                    sigma.e="numeric",
                    sigma.c="numeric",
                    c="list",
                    phi="numeric"))

parameterList$methods(initialize=function(data,inits=NULL,method=1,plot=FALSE,verbose=TRUE){

    if(verbose)
        cat("  Generating initial values...\n")

    ## mu.s1
    if(is.null(inits$mu.s1)){
        if(verbose)
            cat("    Setting mu.s1 to default value (6).\n")
        mu.s1 <<- 6
    }
    else{
        mu.s1 <<- inits$mu.s1
    }

    ## s and e
    if(method==1){
        ## Data based
        if(is.null(inits$s) && is.null(inits$e)){
            tmp <- start.and.end(data,6/data$delta,plot)

            s <<- data$times[tmp$start + 6/data$delta] - runif(data$J,0,data$delta)

            e <<- data$times[tmp$end + 6/data$delta] - runif(data$J,0,data$delta)
        }

        else if(is.null(inits$s) || is.null(inits$e))
            stop("You can only specify initial values for both or neither of s and e.\n")

        else{
            s <<- inits$s
            e <<- inits$e
        }
    }
    else if(method==2){
        ## Constant 24 hour period and 12 hour days with added noise
        s <<- 6 + cumsum(c(0,rep(24,data$J-1))) + rnorm(data$J,0,.5)
        e <<- s + 12 + rnorm(data$J,0,.5)
    }
    else if(method==3){
        ## Constant 25 hour period and 12 hour days with added noise
        s <<- 6 + cumsum(c(0,rep(25,data$J-1))) + rnorm(data$J,0,.5)
        e <<- s + 12 + rnorm(data$J,0,.5)
    }
    else{
        stop("There are only three methods available for choosing s and e.\n")
    }

    ## mu.e
    if(verbose)
        cat("    Computing mean activity length...\n")

    mu.e <<- sapply(1:length(data$treat),function(j){
        mean(e[data$jtreat==j]-s[data$jtreat==j])
    })
    
    ## A
    if(verbose)
        cat("    Computing activity indicators...\n")
    A1 <- outer(data$times,s-data$delta,">")
    A2 <- outer(data$times,e+data$delta,"<")

    if(!is.null(data$skip)){
        A1[,data$skip] <- 0
        A2[,data$skip] <- 0
    }

    A <<- apply(A1*A2,1,sum)

    ## mu
    if(is.null(inits[["mu"]]))
        mu <<- cbind(rep(mean(data$Y[which(A==0)]),data$ntreat),
                     rep(mean(data$Y[which(A==1)]),data$ntreat))
    else
        mu <<- inits$mu

    ## alpha
    if(is.null(inits$alpha)){
        V <- cbind(rep(var(data$Y[which(A==0)]),data$ntreat),
                   rep(var(data$Y[which(A==1)]),data$ntreat))

        alpha <<- mu^2/(V-mu)
    }
    else
        alpha <<- inits$alpha

    ## lambda
    if(is.null(inits$lambda)){
        if(verbose)
            cat("    Generating initial values for lambda...\n")

        treats <- rep(1:length(data$treat),data$n*data$treat)

        lambda <<- matrix(NA,nrow=data$N,ncol=2)

        for(t in 1:data$N){
            shape <- data$Y[t]+alpha[treats,A[t]+1]
            rate <- 1 + alpha[treats,A[t]+1]/mu[treats,A[t]+1]

           lambda[t,A[t]+1] <<- rgamma(1,shape=shape,rate=rate)
        }

    }
    else{
        lambda <<- inits$lambda
    }

    ## missing lambda
    if(is.null(inits$lambda)){
        tmp0 <- which(A==0)
        lambda[tmp0,2] <<- rgamma(length(tmp0),shape=alpha[treats[tmp0],2],
                                  rate=alpha[treats[tmp0],2]/mu[treats[tmp0],2])

        tmp1 <- which(A==1)
        lambda[tmp1,1] <<- rgamma(length(tmp1),shape=alpha[treats[tmp1],1],
                                  rate=alpha[treats[tmp1],1]/mu[treats[tmp1],1])
    }

    ## sigma.e
    if(is.null(inits$sigma.e))
        sigma.e <<- sd(e-s-12)
    else
        sigma.e <<- inits$sigma.e

    ## c
    if(is.null(inits$c)){
        c <<- lapply(1:data$nsegments,function(k){
            stmp <- s[data$segstart[k]:(data$segend[k]+1)]
            x <- stmp[-1]-stmp[-length(stmp)]-24
            lm(x ~ data$Psi[[k]]-1)$coeff
        })
    }
    else
        c <<- inits$c

    ## sigma.c
    # sigma.c <<- rep(1,length(data$segments2))
    sigma.c <<- 1
    
    ## phi
    if(is.null(inits$phi)){
        if(verbose)
            cat("    Computing phi...\n")
        phi <<- unlist(lapply(1:data$nsegments,function(k){
            as.numeric(data$Psi[[k]] %*% c[[k]])
        }))
    }
    else
        phi <<- inits$phi

    ## sigma.s
    if(is.null(inits$sigma.s)){
        sigma.s <<- sd(s[-1]-s[-data$J]-24-phi)
    }
    else
        sigma.s <<- inits$sigma.s
})

start.and.end <- function(data,m=NULL,plot=FALSE){
    ## Define width of window
    if(is.null(m))
        m <- 6/data$delta

    ## medY <- median(data$Y)

    ## active.smooth <- sapply((m+1):(data$N-m+1),function(t){
    ##     sum((data$Y[t+(0:(m-1))] > medY))+
    ##         sum((data$Y[t-(m:1)] < medY)) -
    ##             sum((data$Y[t+(0:m-1)] < medY)) -
    ##                 sum((data$Y[t-(m:1)] > medY))
    ## })

    active.smooth <- sapply((m+1):(data$N-m+1),function(t){
        sum((data$Y[t+(0:(m-1))])) -
            sum((data$Y[t-(m:1)])) 
    })

    ## Plot smoothed activity
    if(plot)
        plot((m+(1:length(active.smooth)))*data$delta,active.smooth,type="h")

    ## Identify start times
    pos <- list(c(1,NA))

    k <- 1
    t <- min(which(active.smooth>0))

    while(t < length(active.smooth)){
        if(active.smooth[t]>0){
            t <- t+1
        }
        else{
            pos[[k]][2] <- t

            while(active.smooth[t] <= 0 && t<length(active.smooth)){
                t <- t+1
            }

            k <- k+1

            pos[[k]] <- c(t,NA)
        }
    }
    
    if(is.na(pos[[k]][2])) pos[[k]] <- NULL

    diff <- sapply(pos,function(t){
        sum(active.smooth[t[1]:t[2]])
    })

    diff.ranks <- rank(diff,ties="first")

    start <- rep(NA,data$J-length(data$skip))
    tmp <- pos[[which(diff.ranks==length(diff))]]
    start[1] <- tmp[1] + which.max(active.smooth[tmp[1]:tmp[2]]) - 1

    j <- 2
    l <- 1

    while(j < (data$J+1-length(data$skip))){

        tmp <- pos[[which(diff.ranks==length(diff)-l)]]
        tmp <- tmp[1] + which.max(active.smooth[tmp[1]:tmp[2]]) - 1
        l <- l+1

        if(min(abs(tmp-start),na.rm=TRUE)>12/data$delta){
            start[j] <- tmp
            j <- j+1
        }
    }

    start <- sort(start)


    ## start1 <- sapply(length(diff):(length(diff)-data$J+1),function(i){
    ##     floor(mean(pos[[which(diff.ranks==i)]]))
    ## })
    ## start1 <- sort(start1)

    end <- sapply(1:(data$J-1-length(data$skip)),function(j){
        start[j] + which.min(active.smooth[start[j]:start[j+1]])-1
    })

    end <- c(end,start[data$J-length(data$skip)] +
             which.min(active.smooth[-(1:start[data$J-length(data$skip)])]))

    ## Add start and end times for skipped periods
    for(j in data$skip){
        start <- c(start[1:(j-1)],
                   floor((2*end[j-1] + start[j])/3),
                   start[-(1:(j-1))])
                   
        end <- c(end[1:(j-1)],
                 floor((start[j] + start[j+1])/2),
                 end[-(1:(j-1))])
    }

    ## Plot estimated start and end times
    if(plot){
        points((start+m)*data$delta,active.smooth[start],col="green",pch=16)
        points((end+m)*data$delta,active.smooth[end],col="red",pch=16)
    }

    return(list(start=start,
                end=end))
}
