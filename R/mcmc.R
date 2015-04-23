##' MCMC wrapper for spline Model 1. (Version 1)
##'
##' This version of the MCMC algorithm augments the data with lambda
##' (the individual time mean counts).
##' @title MCMC for Spline MOdel
##' @param data Data list
##' @param niter Number of iterations
##' @param inits Initial value list
##' @param priors List of parameters of the prior distributions
##' @param proposals List of parameters of the proposal distributions
##' @param updates List of parameters to update (if NULL then update all)
##' @param nadapt Length of adapting phase
##' @param ladapt Iterations between adapting for multivariate parameters
##' @param temper Tempering value
##' @return List of output
##' @author Simon Bonner
mcmc <- function(data,niter,inits=NULL,priors=NULL,proposals=NULL,updates=NULL,
                 nadapt=1000,ladapt=500,radapt=1.01,temper=1,model,debug=FALSE,verbose=TRUE){

    if(debug)
        browser()

    if(verbose)
        cat("  Running MCMC sampling (It\'s time for a coffee!)...\n")

    ## Define prior parameters
    if(is.null(priors))
        priors <- priorParameters(data)

    ## Define proposal paramters
    if(is.null(proposals))
        proposals <- propParameters(data,priors)

    ## Set parameters to update

    if(is.null(updates))
        updates <- c("c","mu","alpha",#"mu_s1",
                     "mu_e","sigma_s","sigma_e",
                     "sigma_c","lambda","s","e")

    ## Initialize parameter list
    if(is.null(inits))
        inits <- parameterList$new(data,method=1)

    params <- inits
    
    ## Initialize storage
    store <- traceList$new(niter,data)

    ## Store initial values
    store$storeAll(0,params)

    ## Acceptance monitors
    accept <- acceptList$new(niter,data)

    ## Run chain
    flush(stdout())
    cat("\n")
    cat("Starting MCMC sampling:\n")
    pb <- txtProgressBar(0,niter,style=3)

    for(i in 1:niter){
        setTxtProgressBar(pb,i)
        flush(stdout())

        ##### MCMC Updates #####
        ## Update mu_s1
        ## if("mu_s1" %in% updates){
        ##     updateMu.s1(data,params,priors,proposals,model=model)
        ##     store$storeMu.s1(i,params)
        ## }
           
        ## Update sigma_s
        if("sigma_s" %in% updates){
            updateSigma.s(data,params,priors,proposals,model=model)
            store$storeSigma.s(i,params)
        }

        ## Update mu_e
        if("mu_e" %in% updates){
            updateMu.e(data,params,priors,proposals)
            store$storeMu.e(i,params)
        }

        ## Update sigma_e
        if("sigma_e" %in% updates){
            updateSigma.e(data,params,priors,proposals)
            store$storeSigma.e(i,params)
        }

        ## Update sigma.c
        if("sigma_c" %in% updates){
            updateSigma.c(data,params,priors,proposals)
            store$storeSigma.c(i,params)
        }
            
        ## Update c
        if("c" %in% updates){
            updateC(data,params,priors,proposals,accept,model)
            store$storeC(i,params)
        }

        ## Update mu
        if("mu" %in% updates){
            updateMu(data,params,priors,proposals,accept)
            store$storeMu(i,params)
        }

        ## Update alpha
        if("alpha" %in% updates){
            updateAlpha(data,params,priors,proposals,accept)
            store$storeAlpha(i,params)
        }

        ## Update lambda
        if("lambda" %in% updates){
            updateLambda(data,params,priors,proposals,temper=1/temper)
        }

        ## Update s
        if("s" %in% updates){
            updateS(data,params,priors,proposals,temper=1/temper,model=model)
            store$storeS(i,params)
            #store$computeSstar(i,data$delta,params)
        }

        ## Update e
        if("e" %in% updates){
            updateE(data,params,priors,proposals,temper=1/temper,accept)
            store$storeE(i,params)
            #store$computeEstar(i,data$delta,params)
        }

        ##### Tuning #####
        if(i < nadapt){
            
            ## c
            ## for(k in 1:data$nsegments){
            ##     if(accept$getc(k=k)==1)
            ##         proposals$c$Sigma[[k]] <- radapt * proposals$c$Sigma[[k]]
            ##     else
            ##         proposals$c$Sigma[[k]] <- 1/radapt * proposals$c$Sigma[[k]]
            ## }
            if(i %% ladapt==0){
                for(k in 1:data$nsegments){
                    proposals$c$Sigma[[k]] <- 2.38^2/ncol(data$Psi[[k]]) *
                        var(store$c[[k]][(i-ladapt+1):i,])
                }
            }
            
            ## alpha
            for(j in 1:data$ntreat){
                for(k in 1:2){
                    if(accept$getalpha(j,k)==1)
                        proposals$alpha$sigma[j,k] <- radapt * proposals$alpha$sigma[j,k]
                    else
                        proposals$alpha$sigma[j,k] <- 1/radapt * proposals$alpha$sigma[j,k]
                }
            }
            
            ## e
            proposals$e$sigma <- proposals$e$sigma *
                (radapt * accept$e[i,] + 1/radapt * (1-accept$e[i,]))
        }
        
        ## ##### Compute likelihood and posterior density #####
        ## lhd <- lhd2(data,params)
        ## prior <- prior.all(params,priors)
        ## posterior <- lhd + prior

        ## store$storeLhd(i,lhd)
        ## store$storePrior(i,prior)
        ## store$storePosterior(i,posterior)
    }

    ## Convert stored output to coda objects
    trace <- list()

    if(length(store$c)==1)
        trace$c <- list(coda::as.mcmc(store$c[[1]]))
    else
        trace$c <- sapply(store$c,coda::as.mcmc)
    
    trace$phi <- lapply(1:data$nsegments,function(k){
        coda::as.mcmc(trace$c[[k]] %*% t(data$Psi[[k]]))
    })

    trace$mu <- coda::as.mcmc(store$mu)

    trace$alpha <- coda::as.mcmc(store$alpha)

    trace$cv <- sqrt(1 + trace$mu/trace$alpha)

    trace$sigma.s <- coda::as.mcmc(store$sigma.s)

    trace$sigma.e <- coda::as.mcmc(store$sigma.e)

    trace$s <- coda::as.mcmc(store$s)
    trace$sstar <- data$delta*ceiling(trace$s/data$delta)

    trace$e <- coda::as.mcmc(store$e)
    trace$estar <- data$delta*ceiling(trace$e/data$delta)

    cat("\n\n")
    
    return(list(inits=inits,
                curr=params,
                trace=trace,
                accept=accept,
                proposals=proposals))
}
