updateC.mh.model1 <- function(data,params,priors,proposals,accept){
    ## Update coefficients for each spline segment separately

    for(k in 1:data$nsegments){
        ## Generate proposal
        c.prop <- params$c[[k]] +
            MASS::mvrnorm(1,rep(0,length(params$c[[k]])),proposals$c$Sigma[[k]])

        phi.prop <- as.numeric(data$Psi[[k]] %*% c.prop)
        
        ## Compute MH ratio
        ## if(k < data$nsegments)
            x <- params$s[(data$segstart[k]:data$segend[k]) + 1] -
                params$s[data$segstart[k]:data$segend[k]] - 24
        ## else
            ## x <- params$s[(data$segstart[k]+1):data$segend[k]] -
            ##     params$s[data$segstart[k]:(data$segend[k] - 1)] - 24
            
        
        ## 1) Proposal contributions
        llik.prop <- sum(dnorm(x,phi.prop,params$sigma.s,log=TRUE))
        
        lprior.prop <- prior.c(c.prop,k,params$sigma.c,priors)

        ## 2) Current contributions
        llik.curr <- sum(dnorm(x,params$phi[data$segstart[k]:data$segend[k]],
                               params$sigma.s,log=TRUE))

        lprior.curr <- prior.c(params$c[[k]],k,params$sigma.c,priors)
        
        lalpha <- (llik.prop + lprior.prop) - (llik.curr + lprior.curr)
        
        if(!is.finite(lalpha))
            browser()
        
        ## Accept or reject proposal
        if(log(runif(1)) < lalpha){
            params$c[[k]] <- c.prop
            params$phi[data$segstart[k]:data$segend[k]] <- phi.prop
            accept$setc(k,1)
        }
        else{
            accept$setc(k,0)
        }
    }
}

updateC.gibbs.model1 <- function(data,params,priors,proposals,accept){
    
    ## Update coefficients for each spline segment separately
    for(k in 1:data$nsegments){
        ## Compute differences between start times
        x <- params$s[(data$segstart[k]:data$segend[k]) + 1] -
            params$s[data$segstart[k]:data$segend[k]] - 24
        
        V <- solve(params$sigma.s^2 *
                   diag(1/(priors$c$v0 * priors $c$V0[[k]] +
                           params$sigma.c^2 * priors$c$V1[[k]])) +
                   t(data$Psi[[k]]) %*% data$Psi[[k]])
        
        mu <-  V %*% t(data$Psi[[k]]) %*% x
        
        ## Generate new value
        params$c[[k]] <- MASS::mvrnorm(1,mu=mu,Sigma=params$sigma.s^2 * V)
        
        ## Update phi
        params$phi[data$segstart[k]:data$segend[k]] <-
            as.vector(data$Psi[[k]] %*% params$c[[k]])
    }
}

updateC.gibbs.model2 <- function(data,params,priors,proposals,accept){
    
    ## Compute differences between start times
    x <- params$s[-1] - 6 - 24 * 1:(data$J-1)

    L <- 1* (outer(1:(data$J-1),1:(data$J-1),"-") >= 0)

    M <- L %*% data$Psi[[1]]

    sigma <- (priors$c$v0 * priors $c$V0[[1]] + params$sigma.c^2 %*% priors$c$V1[[1]])

    V <- solve(params$sigma.s^2 * diag(1/sigma) + t(M) %*% M)
               
    mu <-  V %*% t(M) %*% x
        
    ## Generate new value
    params$c[[1]] <- MASS::mvrnorm(1,mu=mu,Sigma=params$sigma.s^2 * V)
        
    ## Update phi
    params$phi <- as.vector(data$Psi[[1]] %*% params$c[[1]])
}

updateC.gibbs <- function(data,params,priors,proposals,accept,model){
    if(model=="Model1")
        updateC.gibbs.model1(data,params,priors,proposals,accept)
    else if(model=="Model2")
        updateC.gibbs.model2(data,params,priors,proposals,accept)
    else
        stop(paste0("Unknown model: ",model,".\n\n"))
}

updateC <- function(data,params,priors,proposals,accept,model=model,gibbs=TRUE){
    if(gibbs)
        updateC.gibbs(data,params,priors,proposals,accept,model)
    else
        updateC.mh.model1(data,params,priors,proposals,accept)
}
        
