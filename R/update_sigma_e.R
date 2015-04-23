updateSigma.e <- function(data,params,priors,proposals){
    ## Gibbs step
    shape <- (priors$sigma.e$nu.0 + data$J)/2

    SSE <- sum((params$e-params$s-params$mu.e[data$jtreat])^2)

    rate <- (priors$sigma.e$nu.0 * priors$sigma.e$sigma.0^2 + SSE)/2

    params$sigma.e <- 1/sqrt(rgamma(1,shape,rate=rate))
}

    

    
    
    
