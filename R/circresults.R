circresults <- function(output,
                    data,
                    niter,
                    nburnin=NULL,
                    nthin=NULL,
                    plotdir="Plots",
                    verbose=TRUE,
                    debug=FALSE){

    if(debug)
        browser()
    
    if(verbose)
        cat("  Summarizing output:\n")

    if(is.null(nburnin))
        nburnin <- ceiling(niter/10)            # Iterations removed before summarizing

    if(is.null(nthin))
        nthin <- max(floor(niter/1000),1)              # For plots only

    ## Create directories
    if(!file.exists(plotdir))
        dir.create(plotdir)
    else
        cat("    Warning: all files in", plotdir,"are about to be destroyed.\n\n")

    ## Initialize summaries list
    summlist <- list()

    ## Plot traces and compute summary statistics for each set of parameters
    ## 1) Coefficients of the spline functions (c)
    cat("    Spline Coefficients...\n")
    x11()
    
    par(mfrow=c(ceiling(data$nsegments/2),2))

    for(k in 1:data$nsegments)
        matplot(seq(1,niter+1,nthin),
                window(output$trace$c[[k]],thin=nthin),type="l",
                main=paste0("Segment ",k,": Trace of c"),
                xlab="Iteration",ylab="c")

    dev.copy2pdf(file=file.path(plotdir,paste0("trace_c_",inter,".pdf")))

    summlist$c <- lapply(1:data$nsegments,function(k){
        summary(window(output$trace$c[[k]],start=nburnin+1))
    })

    ## 2) Phase shift parameters (phi)
    cat("    Phase shift parameters...\n")
    
    x11()
    par(mfrow=c(ceiling(data$nsegments/2),2))

    for(k in 1:data$nsegments)
        matplot(seq(1,niter+1,nthin),type="l",
                window(output$trace$phi[[k]],thin=nthin),
                xlab="Iteration",ylab="Phi")

    dev.copy2pdf(file=file.path(plotdir,paste0("trace_phi_",inter,".pdf")))

    summlist$phi <- lapply(1:data$nsegments,function(k){
        summary(window(output$trace$phi[[k]],start=nburnin+1))
    })

    summlist$phi.all <- list(do.call("rbind",lapply(summlist$phi,function(summ) summ[[1]])),
                             do.call("rbind",lapply(summlist$phi,function(summ) summ[[2]])))

    ## 2) Mean activity levels (mu)
    cat("    Mean activity levels...\n")

    x11()
    par(mfrow=c(2,1))
    matplot(seq(1,niter+1,nthin),
            window(output$trace$mu[,1:4],thin=nthin),
            type="l",ylim=c(0,max(output$trace$mu[,1:4])),
            xlab="Iteration",ylab="mu (Inactive)")

    matplot(seq(1,niter+1,nthin),
            window(output$trace$mu[,5:8],thin=nthin),
            type="l",ylim=c(0,max(output$trace$mu[,5:8])),
            xlab="Iteration",ylab="mu(Active)")

    dev.copy2pdf(file=file.path(plotdir,paste0("trace_mu_",inter,".pdf")))

    summlist$mu <- summary(window(output$trace$mu,start=nburnin+1))

    ## 3) Variance inflation factors for activity levels (alpha)
    cat("    Activity level variance inflation factors...\n")

    x11()
    par(mfrow=c(2,1))
    matplot(seq(1,niter+1,nthin),
            window(output$trace$alpha[,1:4],thin=nthin),
            type="l",ylim=c(0,max(output$trace$alpha[,1:4])),
            xlab="Iteration",ylab="alpha (Inactive)")
    
    matplot(seq(1,niter+1,nthin),
            window(output$trace$alpha[,5:8],thin=nthin),
            type="l",ylim=c(0,max(output$trace$alpha[,5:8])),
            xlab="Iteration",ylab="alpha (Active)")

    dev.copy2pdf(file=file.path(plotdir,paste0("trace_alpha_",inter,".pdf")))

    summlist$alpha <- summary(window(output$trace$alpha,start=nburnin+1))

    ## 3b) Activity coefficient of variation (mu/sd=1+mu/alpha)
    cat("    Activity coefficient of variation (mu/sd)...\n")

    x11()
    par(mfrow=c(2,1))
    matplot(output$trace$cv[,1:4],type="l",
            ylim=c(0,max(output$trace$cv[,1:4])),
            xlab="Iteration",ylab="Coefficient of Variation (Inactive)")
    matplot(output$trace$cv[,5:8],
            type="l",ylim=c(0,max(output$trace$cv[,5:8])),
            xlab="Iteration",ylab="Coefficient of Variation (Active)")

    dev.copy2pdf(file=file.path(plotdir,paste0("trace_cv_",inter,".pdf")))

    summlist$cv <- summary(window(output$trace$cv,start=nburnin+1))

    ## 3) Standard deviation of start time errors (sigma.s)
    cat("    Start time error standard deviations...\n")
    
    x11()
    plot(seq(1,niter+1,nthin),
         window(output$trace$sigma.s,thin=nthin),
         type="l",ylim=c(0,max(output$trace$sigma.s)),
         xlab="Iteration",ylab="sigma_s")

    dev.copy2pdf(file=file.path(plotdir,paste0("trace_sigma_s_",inter,".pdf")))

    summlist$sigma.s <- summary(window(output$trace$sigma.s,start=nburnin+1))

    ## 4) Standard deviation of end time errors (sigma.e)
    cat("    End time error standard deviations...\n")

    x11()
    plot(seq(1,niter+1,nthin),
         window(output$trace$sigma.e,thin=nthin),
         type="l",ylim=c(0,max(output$trace$sigma.e)),
         xlab="Iteration",ylab="sigma_e")
    
    dev.copy2pdf(file=file.path(plotdir,paste0("trace_sigma_e_",inter,".pdf")))
    
    summlist$sigma.e <- summary(window(output$trace$sigma.e,start=nburnin+1))

    ## 5) Activity start times (s)
    cat("    Activity start times...\n")

    x11()
    matplot(seq(1,niter+1,nthin),
            window(t(t(output$trace$s)-apply(output$trace$s,2,mean)),thin=nthin),
            type="l",
            xlab="Iteration",ylab="Start Times",)
    
    summlist$s <- summary(window(output$trace$s,start=nburnin+1))
    dev.copy2pdf(file=file.path(plotdir,paste0("trace_s_",inter,".pdf")))
    
    ## pdf(file=file.path(plotdir,paste0("trace_sstar_",inter,".pdf")))
    
    ## matplot(seq(1,niter+1,nthin),window(t(t(output$trace$sstar)-apply(output$trace$sstar,2,mean)),thin=nthin),
    ##         type="l")
    ## dev.off()

    ## summlist$sstar <- summary(coda::as.mcmc(output$trace$sstar[-(1:nburnin),]))

    ## 5b) Residuals
    ## output$trace$phi1 <- do.call("cbind",output$trace$phi)
    ## output$trace$pred <- cbind(6,output$trace$s[,1:(data$J-1)]+output$trace$phi1+24)
    ## output$trace$resid <- as.matrix((output$trace$s- output$trace$pred)/output$trace$sigma.s)

    ## pdf(file=file.path(plotdir,paste0("boxplot_s_residuals_",inter,".pdf")))

    ## boxplot(output$trace$resid[-(1:nburnin),])
    ## abline(v=cumsum(data$segments)+.5,lty=2,col="blue")
    ## abline(h=c(0,1,-1),col="red",lty=c(1,2,2))

    ## dev.off()

    ## 6) e (Activity end times)
    cat("    Activity end times...\n")

    x11()

    matplot(seq(1,niter+1,nthin),
            window(t(t(output$trace$estar)-apply(output$trace$estar,2,mean)),thin=nthin),
            type="l",
            xlab="Iteration",ylab="End Times")

    dev.copy2pdf(file=file.path(plotdir,paste0("trace_e_",inter,".pdf")))

    summlist$e <- summary(window(output$trace$e,start=nburnin+1))

    ## matplot(seq(1,niter+1,nthin),window(t(t(output$trace$estar)-apply(output$trace$estar,2,mean)),thin=nthin),
    ##         type="l")

    ## dev.copy2pdf(file=file.path(plotdir,paste0("trace_estar_",inter,".pdf")))
    
    ## summlist$estar <- summary(coda::as.mcmc(output$trace$estar[-(1:nburnin),]))

    ##### Plot activity with Summary Statistics #####
    ## Plot actogram
    x11()

    par(mfrow=c(1,1))
    actogram(data,double=TRUE)
    abline(h=max(data$Y)*(data$D-c(8,13,43)),col="blue",lwd=1.5)

    ## Add posterior summary statistics for startpoints
    for(d in 1:data$D){
        tmp <- intersect(which(summlist$s[[1]][,"Mean"]>24*(d-1)),
                         which(summlist$s[[1]][,"Mean"]<24*(d+1)))

        for(j in tmp){
            points(summlist$s[[1]][j,"Mean"]-24*(d-1),
                   (data$D - d + .5)*max(data$Y),pch=16,col="green")

            lines(summlist$s[[2]][j,c("2.5%","97.5%")]-24*(d-1),
                  rep((data$D- d + .5)*max(data$Y),2),col="green")
        }
    }

    tmp <- do.call("c",sapply(1:data$nsegments,function(k) summlist$phi[[k]][[1]][,"Mean"]))
    tmp1 <- summlist$s[[1]][-data$J,"Mean"] + tmp - (0:(data$J-2))*24

    for(j in 1:(data$J-1))
        lines(c(summlist$s[[1]][j,"Mean"]-(j-1)*24,tmp1[j]),
          (data$D-j+c(.5,-.5))*max(data$Y),
              col="blue")

    ## Add posterior summary statistics for endpoints
    for(d in 1:data$D){
        tmp <- intersect(which(summlist$e[[1]][,"Mean"]>24*(d-1)),
                         which(summlist$e[[1]][,"Mean"]<24*(d+1)))

        for(j in tmp){
            points(summlist$e[[1]][j,"Mean"]-24*(d-1),
                   (data$D - d + .5)*max(data$Y),pch=16,col="red")

            lines(summlist$e[[2]][j,c("2.5%","97.5%")]-24*(d-1),
                  rep((data$D- d + .5)*max(data$Y),2),col="red")
        }
    }

    dev.copy2pdf(file=file.path(plotdir,paste0("actogram_c_",inter,".pdf")))

##     ## Histogram of start and end time probabilities
##     pdf(file=file.path(plotdir,paste0("histogram_1_c_",inter,".pdf")))

##         j <- 50

##         par(mfrow=c(2,2))

##         tmp <- seq(min(output$trace$sstar[-(1:nburnin),j]),max(output$trace$sstar[-(1:nburnin),j]),
##                    data$delta)

##         breaks <- seq(min(tmp)-data$delta/2,
##                       max(tmp)+data$delta/2,
##                       data$delta)

##         hist(output$trace$sstar[-(1:nburnin),j],breaks=breaks,
##              ylim=c(0,niter),freq=TRUE,
##              main=paste("Start Point: Day",j),col="grey",
##              xlab="Possible Values",ylab="Percent Weight",axes=FALSE)
##         box(); axis(1)
##         axis(2,at=niter*seq(0,1,.2),label=seq(0,100,20))

##         lines(tmp,niter*data$Y[round(tmp/data$delta)]/max(data$Y),type="h")
##         points(tmp,niter*data$Y[round(tmp/data$delta)]/max(data$Y),pch=16,cex=.5)
##         text(tmp,niter*data$Y[round(tmp/data$delta)]/max(data$Y),
##              paste(data$Y[round(tmp/data$delta)]),pos=3,cex=.75)

##         matplot(as.matrix(cbind(output$trace$s[,j],output$trace$sstar[,j])),type="l")


##         round(table(output$trace$sstar[-(1:nburnin),j])/(niter-nburnin),2)

##         tmp <- seq(min(output$trace$estar[-(1:nburnin),j]),max(output$trace$estar[-(1:nburnin),j]),
##                    data$delta)

##         breaks <- seq(min(tmp)-data$delta/2,
##                       max(tmp)+data$delta/2,
##                       data$delta)

##         hist(output$trace$estar[-(1:nburnin),j],breaks=breaks,ylim=c(0,niter),freq=TRUE,
##              main=paste("End Point: Day",j),col="grey",
##              xlab="Possible Values",ylab="Percent Weight",axes=FALSE)
##         box(); axis(1)
## axis(2,at=niter*seq(0,1,.2),label=seq(0,100,20))

## lines(tmp,niter*data$Y[round(tmp/data$delta)]/max(data$Y),type="h")
## points(tmp,niter*data$Y[round(tmp/data$delta)]/max(data$Y),pch=16,cex=.5)
## text(tmp,niter*data$Y[round(tmp/data$delta)]/max(data$Y),
##      paste(data$Y[round(tmp/data$delta)]),pos=3,cex=.75)

## matplot(as.matrix(cbind(output$trace$e[,j],output$trace$estar[,j])),type="l")

## dev.off()

## round(table(output$trace$sstar[-(1:nburnin),j])/(niter-nburnin),2)

## ## Phi
## pdf(file=file.path(plotdir,paste0("phi_1_c_",inter,".pdf")))

## par(mar=c(5,5,4,2))

## ylim <- c(min(summlist$phi.all[[2]][,"2.5%"]),max(summlist$phi.all[[2]][,"97.5%"]))

## par(mfrow=c(1,1))
## plot(summlist$phi.all[[1]][,"Mean"],pch=16,ylim=ylim,
##      xlab="Index (j)",
##      ylab=expression(paste("Mean Start Time Deviation (",phi[j],")")))

## for(i in 1:nrow(summlist$phi.all[[2]]))
##     lines(rep(i,2),summlist$phi.all[[2]][i,c("2.5%","97.5%")])

## abline(h=0,lty=2,col="grey")

## abline(v=cumsum(data$segments)+.5,col="blue",lwd=1.5)

## text(c(4,10.5,27.5,46),-1.2,labels=c("L:D","D:D","D:D+PB","D:D"))

## dev.off()
    return(summlist)
}
