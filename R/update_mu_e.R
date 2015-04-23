updateMu.e <- function(data,params,priors,proposals){
    ## Update the mean activity length

    params$mu.e <- sapply(1:length(data$treat),function(j){
        tmp <- which(data$jtreat==j)

        K <- length(tmp)
        
        mean <- (sum(params$e[tmp]-params$s[tmp]) +
                 priors$mu.e$mu * params$sigma.e^2/priors$mu.e$sigma^2) /
                     (K + params$sigma.e^2/priors$mu.e$sigma^2)

        var <- params$sigma.e^2 /
            (K + params$sigma.e^2/priors$mu.e$sigma^2)

        rnorm(1,mean,sqrt(var))
    })
}
   
