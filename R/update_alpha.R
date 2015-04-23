updateAlpha <- function(data,params,priors,proposals,accept){
    ## Update alpha conditional on lambda
    
    ## Generate proposals for alpha
    alpha.prop <- matrix(exp(rnorm(2*length(data$treat),
                                   log(params$alpha),proposals$alpha$sigma)),
                         ncol=2)
    
    ## Compute MH ratios
    treats <- rep(1:length(data$treat),data$n*data$treat)

    for(j in 1:length(data$treat)){
        index <- which(treats==j)
        
        lalpha <- sapply(1:2,function(k){
            num <- sum(dgamma(params$lambda[index,k],shape=alpha.prop[j,k],
                              rate=alpha.prop[j,k]/params$mu[j,k],log=TRUE))
            
            num <- num + prior.alpha(alpha.prop[j,k],j,k,priors,log=TRUE) +
                log(alpha.prop[j,k])
            
            
            den <- sum(dgamma(params$lambda[index,k],shape=params$alpha[j,k],
                              rate=params$alpha[j,k]/params$mu[j,k],log=TRUE))
        
            den <- den + prior.alpha(params$alpha[j,k],j,k,priors,log=TRUE) +
                log(params$alpha[j,k])
            
            
            num - den
        })
        
        ## Accept or reject
        for(k in 1:2){
            
            if(log(runif(1)) < lalpha[k]){
                ## if(j==2 && k==2)
                ##     browser()
                
                params$alpha[j,k] <- alpha.prop[j,k] 
                accept$setalpha(j,k,1)
            }
            else{
                accept$setalpha(j,k,0)
            }
        }
    }
}

updateAlpha2 <- function(data,params,priors,proposals,accept){
    ## Update alpha conditional on y directly

    ## Generate proposal
    alpha.prop <- exp(rnorm(2,log(params$alpha),proposals$alpha$sigma))

    ## Compute MH ratios
    lalpha <- sapply(1:2,function(i){
        tmp <- which(params$A==i-1)
        
        num <- sum(dnbinom(data$Y[tmp],size=alpha.prop[i],
                           mu=params$mu[i],log=TRUE))

        num <- num + prior.alpha(alpha.prop[i],i,priors,log=TRUE)

        #num <- num + log(alpha.prop)


        den <- sum(dnbinom(data$Y[tmp],size=params$alpha[i],
                           mu=params$mu[i],log=TRUE))
        
        den <- den + prior.alpha(params$alpha[i],i,priors,log=TRUE)

        #den <- den + log(params$alpha[i])

        num - den
    })

    if(any(!is.finite(lalpha)))
        browser()

    ## Accept or reject
    for(k in 1:2){
        if(log(runif(1)) < lalpha[k]){
            params$alpha[k] <- alpha.prop[k]
            accept$setalpha(k,1)
        }
        else{
            accept$setalpha(k,0)
        }
    }
}
