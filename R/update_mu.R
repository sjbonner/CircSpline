updateMu <- function(data,params,priors,proposals,accept){
    ## Update mu conditional on lambda (Gibbs step)

    ## Compute parameters of the inverse gamma full conditional
    shape <- data$treat * data$n * params$alpha + priors$mu$shape

    treats <- rep(1:length(data$treat),data$n*data$treat)

    rate <- t(sapply(1:length(data$treat),function(j){
        params$alpha[j,] * apply(params$lambda[which(treats==j),],2,sum) +
            priors$mu$rate[j,]
    }))
    
    ## Simulate from full conditional
    params$mu <- matrix(1/rgamma(2*length(data$treat),shape=shape,rate=rate),
                        ncol=2)
}

updateMu2 <- function(data,params,priors,proposals,accept){
    ## Update mu conditional on Y (MH step)
    
    ## Generate proposal
    mu.prop <- exp(rnorm(2,log(params$mu),proposals$mu$sigma))

    ## Compute MH ratios
    lalpha <- sapply(1:2,function(i){
        tmp <- which(params$A==i-1)
        
        n <- length(tmp)
        
        num <- sum(dnbinom(data$Y[tmp],size=params$alpha[i],
                           mu=mu.prop[i],log=TRUE))

        num <- num + prior.mu(mu.prop[i],i,priors,log=TRUE)

        #num <- num + log(mu.prop)

        den <- sum(dnbinom(data$Y[tmp],size=params$alpha[i],
                           mu=params$mu[i],log=TRUE))
        
        den <- den + prior.mu(params$mu[i],i,priors,log=TRUE)

        #den <- den + log(params$alpha[i])

        num - den
    })

    if(any(!is.finite(lalpha)))
        browser()
    
    ## Accept or reject
    for(k in 1:2){
        if(log(runif(1)) < lalpha[k]){
            params$mu[k] <- mu.prop[k]
            accept$setmu(k,1)
        }
        else{
            accept$setmu(k,0)
        }
    }
}
