genSpline <- function(j,norder=3,nbasis=NULL,tol=1e-6){

    ## Define basis
    if(is.null(nbasis))
        nbasis <- norder + max(0,floor((j-norder)/(norder-1)))

    if(nbasis==norder){
        basis <- fda::create.monomial.basis(c(0,1),nbasis=nbasis)
    }
    else{
        basis <- fda::create.bspline.basis(c(0,1),nbasis=nbasis,norder=norder)
    }

    ## Evaluate basis
    X <- fda::eval.basis((0:(j-1))/(j-1),basis=basis)

    ## Compute penalty matrix
    P <- fda::eval.penalty(basis, Lfdobj=1, rng=basis$rangeval)

    ## Perform eigen decomposition
    eigenP <- eigen(P)

    ## Truncate eigen values to identify penalized terms
    d <- eigenP$values

    ## Reparametrize spline design matrix
    X1 <- cbind(1/sqrt(nbasis),
                X %*% eigenP$vectors[,-nbasis]) #%*% diag(1/sqrt(d[-nbasis])))
    
    #X1 <- X %*% eigenP$vectors[,-nbasis] %*% diag(ifelse(d==0,1,1/sqrt(d)))

    ## Compute prior variance
                #V <- rep(1,nbasis-1)
    V <- d[-nbasis]
                

    return(list(X=X1,V=V))
}


