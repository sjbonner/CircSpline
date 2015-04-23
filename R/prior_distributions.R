prior.sigma.c <- function(sigma,priors,log=TRUE,lower.bound=1e-6){
    ## ## Improper prior (1/sigma)
    ## if(any(sigma < lower.bound))
    ##     return(-Inf)
    ## else
    ##     tmp <- -sum(log(sigma))


    ## Exponential prior
    if(any(sigma < lower.bound))
        return(-Inf)
    else
        tmp <- sum(dexp(sigma,1,log=log))

    ## ## Uniform prior
    ## if(sigma<0)
    ##     return(-Inf)
    ## if(sigma>10)
    ##     return(-Inf)
    ## else
    ##     return(0)

    if(log)
        return(tmp)
    else
        return(exp(tmp))

}

prior.c.dexp <- function(c,k,sigma.c,priors,log=TRUE){
    ## Prior distribution for c[[K]]
    ## tmp <- sum(-abs(priors$c$V[[i]] %*% diag(priors$c$d) %*% c)/priors$c$sigma^2)

    tmp <- -sum(abs(c/sigma.c^2))
    
    if(log)
        return(tmp)
    else
        return(exp(tmp))
}

prior.c.norm <- function(c,k,sigma.c,priors,log=TRUE){
    prior.sd <- sqrt(sigma.c^2 * priors$c$V[[k]])
        
    tmp <- sum(dnorm(c,0,prior.sd,log=TRUE))
    
    if(log)
        return(tmp)
    else
        return(exp(tmp))
}

prior.c.first.deriv <- function(c,k,sigma.c,priors,log=TRUE,debug=FALSE){
    if(debug)
        browser()
    
    prior.sd <- sqrt(as.numeric(priors$c$v0 * priors$c$V0[[k]] +
                                sigma.c^2 * priors$c$V1[[k]]))
    
    tmp <- sum(dnorm(c,0,prior.sd,log=TRUE)) 
    
    if(log)
        return(tmp)
    else
        return(exp(tmp))
}

prior.c <- function(c,k,sigma.c,priors,log=TRUE,dist="first.deriv"){
    ## Prior density for c
    if(dist=="norm")
        tmp <- prior.c.norm(c,k,sigma.c,priors)
    else if(dist=="first.deriv")
        tmp <- prior.c.first.deriv(c,k,sigma.c,priors)
    else if(dist=="dexp"){
        tmp <- prior.c.dexp(c[[k]],k,sigma.c,priors)
    }   
    else{
        stop("Unrecognized distribution prior distribution for c.\n\n")
    }
    
    if(log)
        return(tmp)
    else
        return(exp(tmp))
}

prior.c.all <- function(c,sigma.c,priors,log=TRUE,dist="first.deriv",debug=FALSE){
    if(debug)
        browser()
    
    ## Prior density for c
    if(dist=="norm"){
        tmp <- sum(sapply(1:length(c),function(k){
            prior.c.norm(c[[k]],k,sigma.c,priors)
        }))
    }
    else if(dist=="first.deriv"){
        tmp <- sum(sapply(1:length(c),function(k){
            prior.c.first.deriv(c[[k]],k,sigma.c,priors)
        }))
    }
    else if(dist=="dexp"){
        tmp <- sum(sapply(1:length(c),function(k){
            prior.c.dexp(c[[k]],k,sigma.c,priors)
        }))
    }   
    else{
        stop("Unrecognized distribution prior distribution for c.\n\n")
    }
    
    if(log)
        return(tmp)
    else
        return(exp(tmp))
}
               
prior.sigma.s <- function(sigma,priors,log=TRUE){
    ## Prior distribution for sigma_s
    
    shape <- priors$sigma.s$nu.0/2
    rate <- priors$sigma.s$nu.0/2 * priors$sigma.s$sigma.0^2

    dgamma(1/sigma^2,shape=shape,rate=rate,log=log) * 2/sigma^3
}

prior.sigma.e <- function(sigma,priors,log=TRUE){
    ## Prior distribution for sigma_e
    
    shape <- priors$sigma.e$nu.0/2
    rate <- priors$sigma.e$nu.0/2 * priors$sigma.e$sigma.0^2

    dgamma(1/sigma^2,shape=shape,rate=rate,log=log) * 2/sigma^3
}

prior.mu <- function(mu,treat,k,priors,log=TRUE){
    ## Prior density for mu[k]

    dgamma(1/mu,shape=priors$mu$shape[treat,k],
           rate=priors$mu$rate[treat,k],log=log) *
               1/mu^2
}

prior.mu.all <- function(mu,k,priors,log=TRUE){
    tmp <- sum(dgamma(1/mu,shape=priors$mu$shape,
                      rate=priors$mu$rate,log=TRUE) * 1/mu^2) 

    if(log)
        return(tmp)
    else
        return(exp(tmp))
}

prior.alpha <- function(alpha,treat,i,priors,log=TRUE){
    ## Prior density for alpha

    dgamma(alpha,shape=priors$alpha$shape[treat,i],
           rate=priors$alpha$rate[treat,i],log=log)
}

prior.alpha.all <- function(alpha,priors,log=TRUE){
    tmp <- sum(dgamma(alpha,shape=priors$alpha$shape,rate=priors$alpha$rate,log=TRUE))

    if(log)
        return(tmp)
    else
        return(exp(tmp))
}


prior.all <- function(params,priors,log=TRUE){
    ## Prior density for all parameters
    tmp <- prior.c.all(params$c,params$sigma.c,priors=priors) +
        prior.sigma.s(params$sigma.s,priors=priors) +
            prior.sigma.e(params$sigma.e,priors=priors) +
                prior.mu.all(params$mu,priors=priors) +
                    prior.alpha.all(params$alpha,priors=priors)

    if(log)
        return(tmp)
    else
        return(exp(tmp))
}
