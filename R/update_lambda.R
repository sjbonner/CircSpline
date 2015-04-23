updateLambda <- function(data,params,priors,proposals,temper){
    ## Compute parameters of the gamma full conditional distribution

    treats <- rep(1:length(data$treat),data$n*data$treat)
    
    ## 1) Inactive
    shape1 <- (1-params$A) * temper * data$Y + params$alpha[treats,1]
    rate1 <- (1-params$A) * temper + params$alpha[treats,1]/params$mu[treats,1]

    ## 1) Active
    shape2 <- params$A * temper * data$Y + params$alpha[treats,2]
    rate2 <- params$A * temper + params$alpha[treats,2]/params$mu[treats,2]


    ## Simulate from full conditional
    params$lambda <- cbind(rgamma(data$N,shape=shape1,rate=rate1),
                           rgamma(data$N,shape=shape2,rate=rate2))
}
