updateSigma.s <- function(data,params,priors,proposals,model){
    if(model=="Model1")
        updateSigma.s.model1(data,params,priors,proposals)
    else if(model=="Model2")
        updateSigma.s.model2(data,params,priors,proposals)
    else
        stop(paste0("Unknown model: ",model,".\n\n"))
}

updateSigma.s.model1 <- function(data,params,priors,proposals){
    ## Gibbs step
    shape <- (priors$sigma.s$nu.0 + data$J)/2

    mu.s <- c(6,params$s[-data$J] + 24 + params$phi)
    
    SSE <- sum((params$s - mu.s)^2)

    rate <- (priors$sigma.s$nu.0 * priors$sigma.s$sigma.0^2 + SSE)/2

    params$sigma.s <- 1/sqrt(rgamma(1,shape,rate=rate))
}

updateSigma.s.model2 <- function(data,params,priors,proposals){
    ## Gibbs step
    shape <- (priors$sigma.s$nu.0 + data$J)/2

    mu.s <- cumsum(c(6,24 + params$phi))
    
    SSE <- sum((params$s - mu.s)^2)

    rate <- (priors$sigma.s$nu.0 * priors$sigma.s$sigma.0^2 + SSE)/2

    params$sigma.s <- 1/sqrt(rgamma(1,shape,rate=rate))
}
