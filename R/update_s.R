updateS <- function(data,params,priors,proposals,temper,model){
    if(model=="Model1")
        updateS.model1(data,params,priors,proposals,temper)
    else if(model=="Model2")
        updateS.model2(data,params,priors,proposals,temper)
    else
        stop(paste0("Unknown model: ",model,".\n\n"))
}

updateS.model1 <- function(data,params,priors,proposals,temper){
    ## Update s conditional on lambda

    ## Model 1: E(s_k|s_{k-1})=s_{k-1} + 24 + \phi_{k-1}

    for(d in 1:data$J){
        ## Generate proposal for continuous start time
        if(d==1){
            mu.s <- .5 * (params$mu.s1 + params$s[2]-24-params$phi[1])
            sigma.s <- params$sigma.s/sqrt(2)
        }
        else if(d==data$J){
            mu.s <- params$s[data$J-1] + 24 + params$phi[data$J-1]
            sigma.s <- params$sigma.s
        }
        else{
            mu.s <- .5 * ((params$s[d-1] + params$phi[d-1]) +
                         (params$s[d+1] - params$phi[d]))
            sigma.s <- params$sigma.s/sqrt(2)
        }

        s.prop <- rnorm(1,mu.s,sigma.s)

        ## Identify current and proposed discrete start time
        start <- pmax(1,pmin(ceiling(params$s[d]/data$delta),data$N))
        start.prop <- pmax(1,pmin(ceiling(s.prop/data$delta),data$N))

        ## Identify changes in activity
        if(! d %in% data$skip){
            ## 1) Periods which become active under proposal
            active <- ifelse(start.prop < start,start.prop:(start-1),NA)
            
            ## 2) Periods which become inactive under proposal
            inactive <- ifelse(start.prop > start,start:(start.prop-1),NA)
        }
        
        ## Compute MH ratio
        if(s.prop > params$e[d] || (d>1 && s.prop < params$e[d-1])){
            lalpha1 <- -Inf
        }
        else if(!is.na(active)){
            num <- temper * sum(dpois(data$Y[active],params$lambda[active,2],log=TRUE))
            
            den <- temper * sum(dpois(data$Y[active],params$lambda[active,1],log=TRUE))

            lalpha1 <- num - den

        }
        else if(!is.na(inactive)){
            num <- temper * sum(dpois(data$Y[inactive], params$lambda[inactive,1],log=TRUE))
            
            den <- temper * sum(dpois(data$Y[inactive],params$lambda[inactive,2],log=TRUE))

            lalpha1 <- num - den

        }
        else{
            lalpha1 <- 0
        }

        ## 2) Contribution of end times
        num <- dnorm(params$e[d],s.prop + params$mu.e[data$jtreat[d]],
                     params$sigma.e,log=TRUE)
        den <- dnorm(params$e[d],params$s[d] + params$mu.e[data$jtreat[d]],
                     params$sigma.e,log=TRUE)

        lalpha2 <- num - den

        lalpha <- lalpha1 + lalpha2

        ## Accept or reject proposal
        if(log(runif(1)) < lalpha){
            params$s[d] <- s.prop
            
            params$A[inactive] <- 0
            params$A[active] <- 1
        }
    }
}

updateS.model2 <- function(data,params,priors,proposals,temper){
    ## Update s conditional on lambda

    ## Model 2: E(s_k|E(s_{k-1}))=E(s_{k-1}) + 24 + \phi_{k-1}

    ## Compute 
    mu.s <- cumsum(c(params$mu.s1,24+params$phi))

    for(d in 1:data$J){
        ## Generate proposal for continuous start time
        s.prop <- rnorm(1,mu.s[d],params$sigma.s)

        ## Identify current and proposed discrete start time
        start <- pmax(1,pmin(ceiling(params$s[d]/data$delta),data$N))
        start.prop <- pmax(1,pmin(ceiling(s.prop/data$delta),data$N))

        ## Identify changes in activity
        if(! d %in% data$skip){
            ## 1) Periods which become active under proposal
            active <- ifelse(start.prop < start,start.prop:(start-1),NA)
            
            ## 2) Periods which become inactive under proposal
            inactive <- ifelse(start.prop > start,start:(start.prop-1),NA)
        }
        
        ## Compute MH ratio
        if(s.prop > params$e[d] || (d>1 && s.prop < params$e[d-1])){
            lalpha1 <- -Inf
        }
        else if(!is.na(active)){
            num <- temper *
                sum(dpois(data$Y[active],params$lambda[active,2],log=TRUE))
            
            den <- temper *
                sum(dpois(data$Y[active],params$lambda[active,1],log=TRUE))

            lalpha1 <- num - den

        }
        else if(!is.na(inactive)){
            num <- temper *
                sum(dpois(data$Y[inactive], params$lambda[inactive,1],log=TRUE))
            
            den <- temper *
                sum(dpois(data$Y[inactive],params$lambda[inactive,2],log=TRUE))

            lalpha1 <- num - den

        }
        else{
            lalpha1 <- 0
        }

        ## 2) Contribution of end times
        num <- dnorm(params$e[d],s.prop + params$mu.e[data$jtreat[d]],
                     params$sigma.e,log=TRUE)
        den <- dnorm(params$e[d],params$s[d] + params$mu.e[data$jtreat[d]],
                     params$sigma.e,log=TRUE)

        lalpha2 <- num - den

        ## 3) Contribution of proposal
        lalpha3 <- -dnorm(s.prop,mu.s[d],params$sigma.s,log=TRUE) +
            dnorm(params$s[d],mu.s[d],params$sigma.s,log=TRUE)

        # browser()
        
        lalpha <- lalpha1 + lalpha2 + lalpha3

        ## Accept or reject proposal
        if(log(runif(1)) < lalpha){
            params$s[d] <- s.prop
            
            params$A[inactive] <- 0
            params$A[active] <- 1
        }
    }
}

updateS2.model1 <- function(data,params,priors,proposals,temper){
    ## Update s conditional on Y directly

    for(d in 1:data$J){
        ## Generate proposal for continuous start time
        if(d==1){
            mu.s <- .5 * (params$mu.s1 + params$s[2]-24-params$phi[1])
            sigma.s <- params$sigma.s/sqrt(2)
        }
        else if(d==data$J){
            mu.s <- params$s[data$J-1] + 24 + params$phi[data$J-1]
            sigma.s <- params$sigma.s
        }
        else{
            mu.s <- .5 * ((params$s[d-1] + params$phi[d-1]) +
                         (params$s[d+1] - params$phi[d]))
            sigma.s <- params$sigma.s/sqrt(2)
        }

        s.prop <- rnorm(1,mu.s,sigma.s)

        ## Identify current and proposed discrete start time
        start <- pmax(1,pmin(ceiling(params$s[d]/data$delta),data$N))
        start.prop <- pmax(1,pmin(ceiling(s.prop/data$delta),data$N))

        ## Identify changes in activity
        ## 1) Periods which become active under proposal
        active <- ifelse(start.prop < start,start.prop:(start-1),NA)
        
        ## 2) Periods which become inactive under proposal
        inactive <- ifelse(start.prop > start,start:(start.prop-1),NA)

        ## Compute MH ratio
        if(s.prop > params$e[d] || (d>1 && s.prop < params$e[d-1])){
            lalpha1 <- -Inf
        }
        else if(!is.na(active)){
            num <- temper * sum(dnbinom(data$Y[active],
                                        size=params$alpha[2],mu=params$mu[2],log=TRUE))
            
            den <- temper * sum(dnbinom(data$Y[active],
                                        size=params$alpha[1],mu=params$mu[1],log=TRUE))

            lalpha1 <- num - den

        }
        else if(!is.na(inactive)){
            num <- temper * sum(dnbinom(data$Y[inactive],
                                        size=params$alpha[1],mu=params$mu[1],log=TRUE))
            
            den <- temper * sum(dnbinom(data$Y[inactive],
                                        size=params$alpha[2],mu=params$mu[2],log=TRUE))

            lalpha1 <- num - den

        }
        else{
            lalpha1 <- 0
        }


        ## 2) Contribution of end times
        num <- dnorm(params$e[d],s.prop + params$mu.e[data$jtreat[d]],
                     params$sigma.e,log=TRUE)
        den <- dnorm(params$e[d],params$s[d] + params$mu.e[data$jtreat[d]],
                     params$sigma.e,log=TRUE)

        lalpha2 <- num - den

        lalpha <- lalpha1 + lalpha2

        ## Accept or reject proposal
        if(log(runif(1)) < lalpha){
            params$s[d] <- s.prop
            
            params$A[inactive] <- 0
            params$A[active] <- 1
        }
    }
}

updateS3.model1 <- function(data,params,priors,proposals,temper){
    ## Update s conditional on Y directly
    ## Propose sstar and then propose s|sstar

    for(d in 1:data$J){
        ## Identify current discrete start time
        start <- pmax(1,pmin(ceiling(params$s[d]/data$delta),data$N))
        
        ## Generate proposal for discrete start time
        start.prop <- start + data$delta*sample(-proposals$s$r:proposals$s$r,1)

        ## Generate continuous start time
        if(d==1){
            mu.s <- .5 * (params$mu.s1 + params$s[2]-24-params$phi[1])
            sigma.s <- params$sigma.s/sqrt(2)
        }
        else if(d==data$J){
            mu.s <- params$s[data$J-1] + 24 + params$phi[data$J-1]
            sigma.s <- params$sigma.s
        }
        else{
            mu.s <- .5 * ((params$s[d-1] + params$phi[d-1]) +
                         (params$s[d+1] - params$phi[d]))
            sigma.s <- params$sigma.s/sqrt(2)
        }

        s.prop <- mu.s +
            sigma.s*qnorm(runif(1,pnorm(start.prop*data$delta-data$delta,mu.s,sigma.s),
                                pnorm(start.prop*data$delta,mu.s,sigma.s)))

        ## Identify changes in activity
        ## 1) Periods which become active under proposal
        active <- ifelse(start.prop < start,start.prop:(start-1),NA)
        
        ## 2) Periods which become inactive under proposal
        inactive <- ifelse(start.prop > start,start:(start.prop-1),NA)

        ## Compute MH ratio
        ## 1) Contribution from start time
        if(s.prop > params$e[d] || (d>1 && s.prop < params$e[d-1])){
            lalpha1 <- -Inf
        }
        else{
            num <- log(pnorm(start.prop*data$delta,mu.s,sigma.s) -
                       pnorm((start.prop-1)*data$delta,mu.s,sigma.s))

            den <- log(pnorm(start*data$delta,mu.s,sigma.s) -
                       pnorm((start-1)*data$delta,mu.s,sigma.s))

            lalpha1 <- num-den
        }
        
        ## 2) Contribution from changes in activity
        if(!is.na(active)){
            num <- temper * sum(dnbinom(data$Y[active],
                                        size=params$alpha[2],mu=params$mu[2],log=TRUE))
            
            den <- temper * sum(dnbinom(data$Y[active],
                                        size=params$alpha[1],mu=params$mu[1],log=TRUE))

            lalpha2 <- num - den

        }
        else if(!is.na(inactive)){
            num <- temper * sum(dnbinom(data$Y[inactive],
                                        size=params$alpha[1],mu=params$mu[1],log=TRUE))
            
            den <- temper * sum(dnbinom(data$Y[inactive],
                                        size=params$alpha[2],mu=params$mu[2],log=TRUE))

            lalpha2 <- num - den

        }
        else{
            lalpha2 <- 0
        }

        ## 3) Contribution of end times
        num <- dnorm(params$e[d],s.prop + params$mu.e[data$jtreat[d]],
                     params$sigma.e,log=TRUE)
        den <- dnorm(params$e[d],params$s[d] + params$mu.e[data$jtreat[d]],
                     params$sigma.e,log=TRUE)

        lalpha3 <- num - den

        lalpha <- lalpha1 + lalpha2 + lalpha3

        ## Accept or reject proposal
        if(log(runif(1)) < lalpha){
            params$s[d] <- s.prop
            
            params$A[inactive] <- 0
            params$A[active] <- 1
        }
    }
}

