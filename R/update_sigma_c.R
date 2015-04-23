updateSigma.c <- function (data, params, priors, proposals) 
{
    ## MH Step to update sigma.c

    #browser()
    
    ## Generate proposal
    prop <- rlnorm(1, log(params$sigma.c), proposals$sigma.c$sigma)

    ## Compute log MH ratio
    num <- prior.c.all(params$c, prop, priors) +
        prior.sigma.c(prop,priors) + log(prop)
    
    den <- prior.c.all(params$c, params$sigma.c, priors) +
        prior.sigma.c(params$sigma.c,priors) + log(params$sigma.c)
    
    lalpha <- num - den

    if(!is.finite(lalpha))
        browser()
    
    ## Accept or reject
    if (log(runif(1)) < lalpha) 
        params$sigma.c <- prop
}
