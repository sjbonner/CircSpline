acceptList <-
    setRefClass("acceptList",
                fields=list(c="matrix",
                    citer="numeric",
                    alpha="array",
                    alphaiter="matrix",
                    mu="array",
                    muiter="matrix",
                    e="matrix",
                    eiter="numeric"))

acceptList$methods(initialize=function(niter,data){
    ## c
    c <<- matrix(NA,nrow=niter,ncol=data$nsegments)
    citer <<- rep(1,data$nsegments)

    ## alpha
    alpha <<- array(NA,dim=c(niter,length(data$treat),2))
    alphaiter <<- matrix(1,nrow=length(data$treat),ncol=2)

    ## mu
    mu <<- array(NA,dim=c(niter,length(data$treat),2))
    muiter <<- matrix(1,nrow=length(data$treat),ncol=2)

    ## e
    e <<- matrix(NA,nrow=niter,ncol=data$J)
    eiter <<- rep(1,data$J)
})

acceptList$methods(setc=function(k,a){
    c[citer[k],k] <<- a
    citer[k] <<- citer[k] + 1
})

acceptList$methods(getc=function(k,i=NULL){
    if(is.null(i))
        i <- citer[k]-1
    return(c[i,k])
})

acceptList$methods(setalpha=function(treat,k,a){
    alpha[alphaiter[treat,k],treat,k] <<- a
    alphaiter[treat,k] <<- alphaiter[treat,k] + 1
})

acceptList$methods(getalpha=function(treat,k,i=NULL){
    if(is.null(i))
        i <- alphaiter[treat,k]-1
    return(alpha[i,treat,k])
})

acceptList$methods(setmu=function(treat,k,a){
    mu[muiter[k],treat,k] <<- a
    muiter[treat,k] <<- muiter[treat,k] + 1
})

acceptList$methods(getmu=function(treat,k,i=NULL){
    if(is.null(i))
        i <- muiter[treat,k]-1
    return(mu[i,treat,k])
})

acceptList$methods(sete=function(k,a){
    e[eiter[k],k] <<- a
    eiter[k] <<- eiter[k] + 1
})

acceptList$methods(gete=function(k,i=NULL){
    if(is.null(i))
        i <- eiter[k]-1
    return(e[i,k])
})
