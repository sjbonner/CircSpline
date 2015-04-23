updateE <- function(data,params,priors,proposals,temper,accept){
    ## Update e conditional on lambda
    
    ## Generate proposal
    e.prop <- rnorm(data$J,params$e,proposals$e$sigma)

    end <- pmax(1,pmin(ceiling(params$e/data$delta),data$N))
    end.prop <- pmax(1,pmin(ceiling(e.prop/data$delta),data$N))

    ## Identify changes in activity
    ## 1) Times which became active under proposal
    ## active <- sapply(1:data$J,function(d){
    ##     intersect(which(data$times>params$e[d]),which(data$times<e.prop[d]))
    ## })

    active <- lapply(1:data$J,function(d){
        if(!d %in% data$skip && end.prop[d] > end[d])
            (end[d]+1):end.prop[d]
    })

    ## 2) which became inactive under proposal
    ## inactive <- lapply(1:data$J,function(d){
    ##     intersect(which(data$times<params$e[d]),which(data$times>e.prop[d]))
    ## })

    inactive <- lapply(1:data$J,function(d){
        if(!d %in% data$skip && end[d] > end.prop[d])
            (end.prop[d]+1):end[d]
    })

    ## Compute MH ratios
    ## 1) Contribution from end times
    num <- sum(dnorm(e.prop,params$s + params$mu.e[data$jtreat],
                     params$sigma.e,log=TRUE))
    den <- sum(dnorm(params$e,params$s + params$mu.e[data$jtreat],
                     params$sigma.e,log=TRUE))
    lalpha1 <- num - den
        
    ## 2) Contribution from changes in activity
    lalpha2 <- sapply(1:data$J,function(d){
        if(d==data$J)
            num <- ifelse(e.prop[d]< params$s[d],-Inf,0)
        else
            num <- ifelse((e.prop[d] < params$s[d]) + (e.prop[d] > params$s[d+1]),-Inf,0)

        den <- 0
        
        if(!is.null(active[[d]])){
            ## e.prop increases activity
            num <- num +
                temper * sum(dpois(data$Y[active[[d]]],
                                   params$lambda[active[[d]],2],log=TRUE))
            
            den <- den +
                temper * sum(dpois(data$Y[active[[d]]],
                                   params$lambda[active[[d]],1],log=TRUE))
        }
        if(!is.null(inactive[[d]])){
            ## e.prop decreases activity
            num <- num +
                temper * sum(dpois(data$Y[inactive[[d]]],
                                   params$lambda[inactive[[d]],1],log=TRUE))
            
            den <- den +
                temper * sum(dpois(data$Y[inactive[[d]]],
                                   params$lambda[inactive[[d]],2],log=TRUE))
        }

        num - den
    })

    lalpha <- lalpha1 + lalpha2

    ## Accept or reject proposals
    for(j in 1:data$J){
        if(log(runif(1)) < lalpha[j]){
            params$e[j] <- e.prop[j]
            accept$sete(j,1)
            
            params$A[inactive[[j]]] <- 0
            params$A[active[[j]]] <- 1
        }
        else{
            accept$sete(j,0)
        }
    }
}
                     
updateE2 <- function(data,params,priors,proposals,temper,accept){
    ## Update E conditional on Y directly
    
    ## Generate proposal
    #e.prop <- rnorm(data$J,params$s + 12,params$sigma.e)
    e.prop <- rnorm(data$J,params$e,proposals$e$sigma)

    end <- pmax(1,pmin(ceiling(params$e/data$delta),data$N))
    end.prop <- pmax(1,pmin(ceiling(e.prop/data$delta),data$N))

    ## Identify changes in activity
    ## 1) Times which became active under proposal
    ## active <- sapply(1:data$J,function(d){
    ##     intersect(which(data$times>params$e[d]),which(data$times<e.prop[d]))
    ## })

    active <- lapply(1:data$J,function(d){
        if(end.prop[d] > end[d])
            (end[d]+1):end.prop[d]
    })

    ## 2) which became inactive under proposal
    ## inactive <- lapply(1:data$J,function(d){
    ##     intersect(which(data$times<params$e[d]),which(data$times>e.prop[d]))
    ## })

    inactive <- lapply(1:data$J,function(d){
        if(end[d] > end.prop[d])
            (end.prop[d]+1):end[d]
    })

    ## Compute MH ratios
    ## 1) Contribution from end times
    num <- sum(dnorm(e.prop,params$s + 12,params$sigma.e,log=TRUE))
    den <- sum(dnorm(params$e,params$s + 12,params$sigma.e,log=TRUE))
    lalpha1 <- num - den
        
    ## 2) Contribution from changes in activity
    lalpha2 <- sapply(1:data$J,function(d){
        if(d==data$J)
            num <- ifelse(e.prop[d]< params$s[d],-Inf,0)
        else
            num <- ifelse((e.prop[d] < params$s[d]) + (e.prop[d] > params$s[d+1]),-Inf,0)

        den <- 0
        
        if(!is.null(active[[d]])){
            ## e.prop increases activity
            num <- num +
                temper * sum(dnbinom(data$Y[active[[d]]],
                                     size=params$alpha[2],mu=params$mu[2],log=TRUE))
            
            den <- den +
                temper * sum(dnbinom(data$Y[active[[d]]],
                                     size=params$alpha[1],mu=params$mu[1],log=TRUE))
        }
        if(!is.null(inactive[[d]])){
            ## e.prop decreases activity
            num <- num +
                temper * sum(dnbinom(data$Y[inactive[[d]]],
                                     size=params$alpha[1],mu=params$mu[1],log=TRUE))
            
            den <- den +
                temper * sum(dnbinom(data$Y[inactive[[d]]],
                                     size=params$alpha[2],mu=params$mu[2],log=TRUE))
        }

        num - den
    })

    lalpha <- lalpha1 + lalpha2
    
    ## Accept or reject proposals
    for(j in 1:data$J){
        if(log(runif(1)) < lalpha[j]){
            params$e[j] <- e.prop[j]
            accept$sete(j,1)
            
            params$A[inactive[[j]]] <- 0
            params$A[active[[j]]] <- 1
        }
        else{
            accept$sete(j,0)
        }
    }
}
                     

