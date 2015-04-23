lhd.s <- function(data,params,log=TRUE){
    ## Log density of s given parameters

    ## Compute deviances from spline model
    d <- c(params$s[1]-params$mu.s1,params$s[-1]-params$s[-data$J]-24-params$phi)

    ## Log-likelihood contribution
    tmp <- sum(dnorm(d,0,sd=params$sigma.s,log=TRUE))

    if(log)
        return(tmp)
    else
        return(exp(tmp))
}

lhd.e <- function(data,params,log=TRUE){
    ## Log density of e given s and parameters

    ## Compute deviances
    d <- params$e - params$s - 12

    ## Log-likeilhood contribution
    tmp <- sum(dnorm(d,0,sd=params$sigma.e,log=TRUE))

    if(log)
        return(tmp)
    else
        return(exp(tmp))
}

lhd.lambda <- function(data,params,log=TRUE){
    ## Log density of lambda given parameters

    tmp <- sum(dgamma(params$lambda[,1],shape=params$alpha[1],
                      rate=params$alpha[1]/params$mu[1],log=TRUE)) +
               sum(dgamma(params$lambda[,2],shape=params$alpha[2],
                          rate=params$alpha[2]/params$mu[2],log=TRUE))

    if(log)
        return(tmp)
    else
        return(exp(tmp))
}

lhd.Y <- function(data,params,log=TRUE){
    ## Log density of Y given parameters

    ## Compute A
    A1 <- outer(data$times,params$s-data$delta,">")
    A2 <- outer(data$times,params$e+data$delta,"<")
    A <- apply(A1*A2,1,sum)

    ## Log-likelihood contribution
    tmp0 <- which(A==0)
    tmp1 <- which(A==1)

    tmp <- sum(dpois(data$Y[tmp0],params$lambda[tmp0,1],log=TRUE)) +
        sum(dpois(data$Y[tmp1],params$lambda[tmp1,2],log=TRUE))

    if(log)
        return(tmp)
    else
        return(exp(tmp))
}

lhd.Y2 <- function(data,params,log=TRUE){
    ## Log density of Y given parameters
    ## Does not condition on lambda

    ## Compute A
    A1 <- outer(data$times,params$s-data$delta,">")
    A2 <- outer(data$times,params$e+data$delta,"<")
    A <- apply(A1*A2,1,sum)

    ## Log-likelihood contribution
    tmp1 <- which(A==1)
    n <- length(tmp1)

    tmp <- dnbinom(sum(data$Y[-tmp1]),size=(data$N-n)*params$alpha[1],
                   mu=(data$N-n)*params$mu[1],log=TRUE) +
                       dnbinom(sum(data$Y[tmp1]),size=n*params$alpha[1],
                               mu=n*params$mu[1],log=TRUE)
    
    if(log)
        return(tmp)
    else
        return(exp(tmp))
}

lhd <- function(data,params,log=TRUE){
    ## Compute likelihood

    tmp <- lhd.s(data,params) +
        lhd.e(data,params) +
            lhd.lambda(data,params) +
                lhd.Y(data,params)

    if(log)
        return(tmp)
    else
        return(exp(tmp))
}

lhd2 <- function(data,params,log=TRUE){
    ## Compute likelihood
    
    tmp <- lhd.s(data,params) +
        lhd.e(data,params) +
            lhd.Y2(data,params)

    if(log)
        return(tmp)
    else
        return(exp(tmp))
}
