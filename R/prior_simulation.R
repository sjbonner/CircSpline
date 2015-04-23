priorSim <- function(skip,data,priors){
    ## Prior predictive simulation

    ## 0) Simulate variance parameters
    shape <- priors$sigma.s$nu.0/2
    rate <- priors$sigma.s$nu.0*priors$sigma.s$sigma.0^2/2

    sigma.s <- 1/sqrt(rgamma(1,shape=shape,rate=rate))

    shape <- priors$sigma.e$nu.0/2
    rate <- priors$sigma.e$nu.0*priors$sigma.e$sigma.0^2/2

    sigma.e <- 1/sqrt(rgamma(1,shape=shape,rate=rate))
        
    ## 1) Simulate spline
    ## c <- lapply(1:length(data$Psi),function(k){
    ##     l <- ncol(data$Psi[[k]])
    ##     b <- (-1)^(runif(l)<.5) * rexp(l,priors$c$sigma^2)
    ##     diag(priors$c$d[[k]]) %*% t(priors$c$V[[k]]) %*% b
    ## })
    
    c <- lapply(1:data$nsegments,function(k){
        l <- ncol(data$Psi[[k]])
        (-1)^(runif(l)<.5) * rexp(l,priors$c$sigma^2)
    })
    
    phi <- unlist(lapply(1:length(data$Psi),function(k){
        data$Psi[[k]] %*% c[[k]]
    }))

    ## 2) Generate activity start times
    s <- rep(NA,data$J)
    s[1] <- rnorm(1,6,sigma.s)

    for(j in 2:data$J)
        s[j] <- rnorm(1,s[j-1]+24+phi[j-1],sigma.s)

    ## 3) Generate activity end times
    e <- s + 12 + rnorm(data$J,0,sigma.e)

    ## 4) Define activity indicators
    A1 <- outer(data$times,s-data$delta,">")
    A2 <- outer(data$times,e+data$delta,"<")

        if(!is.null(skip)){
        A1[,skip] <- 0
        A2[,skip] <- 0
    }
    
    A <<- apply(A1*A2,1,sum)

    ## 5) Count parameters
    mu <- 1/rgamma(2,shape=priors$mu$shape,rate=priors$mu$rate)
    alpha <- rgamma(2,shape=priors$alpha$shape,rate=priors$alpha$rate)

    ## 6) Simulate lambda
    lambda <- cbind(rgamma(data$N,shape=alpha[1],rate=alpha[1]/mu[1]),
                    rgamma(data$N,shape=alpha[2],rate=alpha[2]/mu[2]))

    ## 7) Simulate data
    Y <- rpois(data$N,(1-A)*lambda[,1] + A*lambda[,2])

    return(list(sigma.s=sigma.s,
                sigma.e=sigma.e,
                c=c,
                s=s,
                e=e,
                mu=mu,
                alpha=alpha,
                lambda=lambda,
                Y=Y))
}                      
    
