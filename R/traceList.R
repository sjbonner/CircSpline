traceList <-
    setRefClass("traceList",
                fields=list(
                    sigma.s="numeric",
                    mu.e="matrix",
                    sigma.c="numeric",
                    sigma.e="numeric",
                    c="list",
                    mu="matrix",
                    alpha="matrix",
                    s="matrix",
                    e="matrix",
                    lhd="numeric",
                    prior="numeric",
                    posterior="numeric"
                    ))

traceList$methods(initialize=function(niter,data){
    sigma.s <<- rep(0,niter+1)

    mu.e <<- matrix(0,niter+1,data$ntreat)

    sigma.c <<- rep(0,niter+1)

    sigma.e <<- rep(0,niter+1)

    c <<- lapply(sapply(data$Psi,ncol),function(col){
        matrix(0,niter+1,col)
    })

    mu <<- matrix(0,niter+1,2*data$ntreat)

    alpha <<- matrix(0,niter+1,2*data$ntreat)

    s <<- matrix(0,niter + 1,data$J)

    e <<- matrix(0,niter + 1,data$J)

    lhd <<- numeric(niter+1)

    prior <<- numeric(niter+1)

    posterior <<- numeric(niter+1)
})

traceList$methods(storeSigma.s=function(i,params){
    sigma.s[i+1] <<- params$sigma.s
})

traceList$methods(storeMu.e=function(i,params){
    mu.e[i+1,] <<- params$mu.e
})

traceList$methods(storeSigma.c=function(i,params){
    sigma.c[i+1] <<- params$sigma.c
})

traceList$methods(storeSigma.e=function(i,params){
    sigma.e[i+1] <<- params$sigma.e
})

traceList$methods(storeC=function(i,params){
    for(k in 1:length(params$c)){
        c[[k]][i,] <<- params$c[[k]]
    }
})

traceList$methods(storeMu=function(i,params){
    mu[i+1,] <<- as.vector(params$mu)
})

traceList$methods(storeAlpha=function(i,params){
    alpha[i+1,] <<- as.vector(params$alpha)
})

traceList$methods(storeS=function(i,params){
    s[i+1,] <<- as.vector(params$s)
})

traceList$methods(storeE=function(i,params){
    e[i+1,] <<- as.vector(params$e)
})

traceList$methods(storeLhd=function(i,lhdin){
    lhd[i+1] <<- lhdin
})

traceList$methods(storePrior=function(i,priorin){
    prior[i+1] <<- priorin
})

traceList$methods(storePosterior=function(i,posteriorin){
    posterior[i+1] <<- posteriorin
})

traceList$methods(storeAll=function(i,params){
    storeSigma.s(i,params)
    storeMu.e(i,params)
    storeSigma.c(i,params)
    storeSigma.e(i,params)
    storeC(i,params)
    storeMu(i,params)
    storeAlpha(i,params)
    storeS(i,params)
    storeE(i,params)
})




    
    
                    
                    
                
