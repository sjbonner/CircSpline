## Load package
library(CircSpline)

## Analysis settings
niter <- 100                        # Number of MCMC iterations
inter <- 30                          # Observation interval in minutes
model <- "Model1"                    # Which model of phi to consider

## Load raw data directly from csv
## data.raw.all <- read.csv("BirdSongEntrain.TonesFinal.csv",
##                          header=TRUE,stringsAsFactors=FALSE)

## Load data included in package
data(BirdSongEntrain.TonesFinal)

## Remove the first and last days (which are partial) for simplicity
data.raw.all <- data.raw.all[-which(data.raw.all[,1]=="9/4/2012"),]
data.raw.all <- data.raw.all[-which(data.raw.all[,1]=="10/26/2012"),]

## Extract data for individual i
i <- 3
data.raw <- data.raw.all[,i+2]

## Sum data over intervals of length inter minutes
delta <- inter/60                       # Interval in hours
data.raw.1 <- apply(matrix(data.raw,nrow=12*delta),2,sum,na.rm=TRUE)

## Remove NaN values for now
data.raw.1[which(is.nan(data.raw.1))] <- 0

## Construct design list
data <- buildData(Y=data.raw.1,
                  D=51,
                  n=24/delta,
                  J=50,
                  treat=c(8,6,29,8),
                  segments=c(8,6,29,6),
                  skip=9)

## Fit spline model
circfit <- circspline(data,niter,verbose=TRUE)
