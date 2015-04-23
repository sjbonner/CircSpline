posterior <- function(data,params,priors,log=TRUE){
    ## Compute posterior density
    tmp <- prior.all(params,priors) + lhd(data,params)

    if(log)
        return(tmp)
    else
        return(exp(tmp))
}
