#HiatusHistoryTemps.R: Code to do all temperature-only analyses for first ERL paper 2017.

rm(list=ls())
graphics.off()
#setwd("C:/Users/Lewan/Documents/Research Projects/Climate Change/_climate with James/HiatusHistory") #not needed if run via project
source("hiatushistoryfuncs.r") #define all required functions
source("hiatushistorymodelfuncs.r")

library(fArma)
library(lattice)
library(stats)
library(car)
library(diagram)
library(latticeExtra) #for layer
library(sp) #for layer sp.points
library(effects)
library(gplots)
library(rgl)
library(nlme)
library(grid)
library(tidyr)
library(plyr)
library(dplyr)
library(strucchange)

#######################################################################################
#function to read all available versions of a data set, 
#   no matter how many there are.
#   Alphabetical order of data files must correspond to their chronological order
#   Return all versions of data sets in a list.
readallversions <- function(datadir,ds) {
  wd <- paste(datadir,ds,sep="/")
  tbr<- paste (wd,list.files(wd,pattern="*.temp"),sep="/")
  r  <- vector("list",length(tbr))       #build a list of all versions of the data 
  for (i in c(1:length(tbr))) {          #now read each file and return result in a list
    temp <- read.table(tbr[i])
    r[[i]] <- data.frame(t=temp[,1],anom=temp[,2])
  }
  return(r)
}

#function rebaselines anomalies to the desired baseline period 
#  and cleans the data to fall within the desired window of use
cleandata <- function(rawvariants) {
  cln <- rawvariants
  for (i in c(1:length(cln))) {
    cln[[i]] <- dplyr::filter(cln[[i]],t>=yrsused[1] & t<=yrsused[2])
    blmean <- mean(filter(cln[[i]],t>=baselnpd[1] & t<=baselnpd[2])$anom)
    cln[[i]]$anom <- cln[[i]]$anom - blmean
  }
  return(cln)
}

#function annualizes data set for convenience
annualize <- function(ds) {
  ca <- aggregate(anom ~ floor(t), data=ds, FUN=mean)
  names(ca)[1] <- "t"
  return(ca)
}

# function to ensure that year always runs to late December, even if provided as integer
mkyr <- function(y) {return(ifelse(abs(floor(y)-y) > 1e-8, y, y+.99))}
# function to verify that a dataset is annualized to guard against silly errors
verifyannual <- function(ds) {try(if(abs(floor(max(ds$t))-max(ds$t)) > 1e-8) stop("DS not annualized")) }

#http://stackoverflow.com/questions/21011672/automatically-add-variable-names-to-elements-of-a-list
listN <- function(...){
  anonList <- list(...)
  names(anonList) <- as.character(substitute(list(...)))[-1]
  anonList
}

#######################################################################################
#major global driver variables for the entire analysis. Some are packaged into list below.
datadir  <- "./HiatusHistory TemperatureDataV1.0"
outputdir <- "./HiatusHistory Temp Output"
yrsused  <- c(1880,2016.999)  #inclusive
yrslscap <- c(1965,2016.999)  #inclusive for landscaping
baselnpd <- c(1981,2010.999)  #inclusive
boftrnd  <- c(1998,1970)      #[1]=beginning of trend for Cowtan retro graph,
#[2]=beginning of global warming trend for continuous slopes
mintrnd    <- 10              #minimum trend length for histocond figures
maxyrlscap <- 24             #maximum years back in time for landscaping
nMC        <- 1000            #number of Monte Carlo realizations

#obtain data sets by checking directories in the data directory
datasets <- list.dirs(datadir,full.names=FALSE,recursive=FALSE) 
names(datasets) <- datasets
#manually add first-available-date
ds1strel <- rep(1879.999,length(datasets))   #default availability from early on
names(ds1strel) <- datasets
ds1strel["CW"]  <- 2013.80                   #November 2013 (for GT test)
ds1strel["BERKELEY"]  <- 2013.90             #December 2013 (for GT test) [Actually March 2014 but that wouldnt be meaningful]
#now keep track of changepoints
cps <- rep(boftrnd[2],length(datasets))      #just in case some explicitly missing (e.g. BERKELEY)
names(cps) <- datasets
cps["GISTEMP"] <- 1970 #SD=3
cps["CW"]      <- 1974 #SD=2
cps["HADCRUT"] <- 1974 #SD=2
cps["NOAA"]    <- 1970 #SD=3

#plotting variables such as color and so on. Some packaged into list below.
plcol  <- c("red","purple4","darkblue","dark green","deepskyblue")
plthincol <-c("mistyrose","plum","azure3","palegreen3","cyan")
plchar <- c(21,22,23,24,25)
names(plcol) <- names(plthincol) <- names(plchar) <- datasets

#for plotting of datelines
hadreldates        <- c(2010.96,   2012.79, 2013.12,2014.54,2015.38)
names(hadreldates) <- c("HadCRUT3","4.1.1", "4.2",  "4.3",   "4.4")

gissreldates       <- c(2010.21,      2011.96,      2012.96,      2015.46)
names(gissreldates)<- c("GISTEMPv2", "GISTEMPv3", "v3/ERSSTv3", "v3/ERSSTv4")

#######################################################################################
#now read, clean, and rebaseline all available variants of each data set
mcresults <- alldata  <- vector("list",length(datasets))
names(mcresults) <- names(alldata) <- datasets      #permit convenient indexing by name
for (i in c(1:length(datasets))) {
  rawvariants   <- readallversions(datadir,datasets[i])
  alldata[[i]]  <- cleandata(rawvariants)    
}

#######################################################################################
#read the literature corpus, omitting ineligible records
lit <- getlit()
# x11()
# hist(lit$Start)
# x11()
# hist(lit$duration)
summary(lit$Start)
summary(lit$duration)
write.csv(lit,"lit.temp.csv")


#now process each paper
pslopes <- prefslopes <- NULL
for (i in 1:dim(lit)[1]) {
  print(c(i,lit$ds[i,]))
  nds <- unique(lit$ds[i,!is.na(lit$ds[i,])])
  for (ds4histo in nds) {
      apaper <- getstat4paper(alldata,ds4histo,lit$Start[i],lit$End[i])
      pslopes <- c(pslopes,apaper$pslope)
      prefslopes <- c(prefslopes,apaper$refslopes)
  }
}
#plot the results of literature analysis
y1 <- hist(pslopes, plot=FALSE,breaks="FD")
y1$density <- y1$density/sum(y1$density)
x11()
plot( y1, col=rgb(0,0,1,1/4), freq=FALSE,  xlab="K/decade",ylab="Density",main=NULL,
        ylim=c(0,.6),xlim=c(-.3,.5))
os <- hist(prefslopes, plot=FALSE,breaks="FD")
os$density <- os$density/sum(os$density)
plot( os, col=rgb(1,0,1,1/4), freq=FALSE,  add=TRUE)
abline(v=0,lty="dashed")

savePlot(filename = paste(outputdir,"/litHistoSlopes.pdf",sep=""), type = "pdf", device = dev.cur(), restoreConsole = TRUE)
savePlot(filename = paste(outputdir,"/litHistoSlopes.eps",sep=""), type = "eps", device = dev.cur(), restoreConsole = TRUE)
savePlot(filename = paste(outputdir,"/litHistoSlopes.png",sep=""), type = "png", device = dev.cur(), restoreConsole = TRUE)


#######################################################################################
#first plot the latest annualized means 
p2f <- 0
if (p2f) {
  fn <- paste(outputdir,"/AlltempsAnnualized.pdf",sep = "")
  if (file.exists(fn)) {
    file.remove(fn)
  }
  pdf(file = fn,height = 7,width = 7)
} else {
  x11(7, 7)
}
xvals <- seq(yrsused[1],floor(yrsused[2]))
tbpds <- names(datasets) #c("GISTEMP","HADCRUT")
plot(xvals,rep(0,length(xvals)),
     type="n",ylim=c(-0.9,.6),xlab="Year",ylab="GMST Anomaly (K)", yaxt="n")
for (i in tbpds) { 
  curannual <- aggregate(anom ~ floor(t), data=alldata[[i]][[length(alldata[[i]])]], FUN=mean)
  lines(xvals,curannual$anom,col=plcol[i],lwd=1)
  #points(xvals,curannual$anom,pch=plchar[i],bg=plcol[i],cex=0.5)
}
#abline(v=c(1980,1985,1990,1995,2000,2005,2010),lty="dotted")
axis(2,at=seq(from=-.8,to=.6,by=.1), las=1)
legend(1900,.5,datasets[tbpds],col=plcol[tbpds],
       lty="solid",pch=plchar[tbpds],pt.bg=plcol[tbpds])
if (p2f) {
  dev.off()
}
#for EPS use p2f <- 0
savePlot(filename = paste(outputdir,"/AlltempsAnnualized.eps",sep = ""),type = "eps",device = dev.cur(),restoreConsole = TRUE)

#plot historically conditioned major data sets for Nature news and views
p2f <- 0
if (p2f) {
  fn <- paste(outputdir,"/NatureNVannualized.pdf",sep = "")
  if (file.exists(fn)) {
    file.remove(fn)
  }
  pdf(file = fn,height = 7,width = 7)
} else {
  x11(7, 7)
}
xvals <- seq(1970,floor(yrsused[2]))
plot(xvals,rep(0,length(xvals)),
     type="n",ylim=c(-0.7,.6),xlab="Year",ylab="GMST Anomaly (K)", yaxt="n")
  #HadCRUT3
  curannual <- aggregate(anom ~ floor(t), data=subset(alldata[["HADCRUT"]][[2]],t>=1970), FUN=mean)
  #lines(lowess(xvals[1:dim(curannual)[1]],curannual$anom),col="red",lwd=2.5)
  #abline(coef=lm(curannual$anom~xvals[1:dim(curannual)[1]])$coefficients,col="red",lwd=2.5)
  lines(xvals[1:dim(curannual)[1]],curannual$anom,col="red",lwd=2.5)
  #points(xvals[1:dim(curannual)[1]],curannual$anom,pch=21,bg="red",cex=1.)
  #HadCRUT4
  curannual <- aggregate(anom ~ floor(t), data=subset(alldata[["HADCRUT"]][[7]],t>=1970), FUN=mean)
  #lines(lowess(xvals,curannual$anom),col="light green",lwd=2.5)
  #abline(coef=lm(curannual$anom~xvals)$coefficients,col="light green",lwd=2.5)
    lines(xvals,curannual$anom,col="light green",lwd=2.5)
  #points(xvals,curannual$anom,pch=22,bg="light green",cex=1.) 
  #GISTEMP
  curannual <- aggregate(anom ~ floor(t), data=subset(alldata[["GISTEMP"]][[length(alldata[["GISTEMP"]])]],t>=1970), FUN=mean)
  #lines(lowess(xvals,curannual$anom),col="blue",lwd=2.5)
  #abline(coef=lm(curannual$anom~xvals)$coefficients,col="blue",lwd=2.5)
  lines(xvals,curannual$anom,col="blue",lwd=2.5)
  #points(xvals,curannual$anom,pch=23,bg="blue",cex=1.)
axis(2,at=seq(from=-.8,to=.6,by=.1), las=1)
legend(1975,.5,c("HadCRUT3","HadCRUT4","GISTEMP"),col=c("red","light green","blue"),
       lty="solid",lwd=3)#,pch=c(21,22,23),pt.bg=c("red","light green","blue"))
if (p2f) {
  dev.off()
}


#######################################################################################
#now plot historically conditioned temperatures
p2f  <- 0 #1 = to file, 0 = screen only
tofr <- 1 #broken trend = 1, continuous trend = 2

#package relevant variables together into lists
plotvars <- listN(plthincol, plcol, p2f, hadreldates, gissreldates)
timevars <- listN(ds1strel, maxyrlscap, boftrnd, mintrnd, yrsused, yrslscap)

#now create historically-conditioned plots with a suitable subset of data sets
plhistocond(alldata, c("GISTEMP","HADCRUT"), tofr, timevars, plotvars)
plhistocond(alldata, c("GISTEMP","HADCRUT","NOAA"), tofr, timevars, plotvars)
plhistocond(alldata, datasets[-3], tofr, timevars, plotvars)
plhistocond(alldata, datasets, tofr, timevars, plotvars)
plhistocond(alldata, c("CW","HADCRUT"), tofr, timevars, plotvars)
plhistocond(alldata, c("GISTEMP"), tofr, timevars, plotvars)
plhistocond(alldata, c("HADCRUT"), tofr, timevars, plotvars)


#now plot trends for the Nature NV piece
alldata4n  <- vector("list",3)
names(alldata4n) <- c("HadCRUT3","HadCRUT4","GISTEMP")
alldata4n[["GISTEMP"]]  <- list(alldata[["GISTEMP"]][[length(alldata[["GISTEMP"]])]])
alldata4n[["HadCRUT3"]] <- list(alldata[["HADCRUT"]][[2]])
alldata4n[["HadCRUT4"]] <- list(alldata[["HADCRUT"]][[7]])
plhistocond4n(alldata4n, c("HadCRUT3","HadCRUT4","GISTEMP"), tofr, timevars, plotvars)


#######################################################################################
#now plot historically conditioned temperatures using 'santer plots'
p2f <- 0
tofr <- 1
trendlength <- 15
startyr <- 1970
endyr <- 2016
plsanter(alldata,c("GISTEMP","HADCRUT"),tofr,plcol,trendlength,startyr,endyr,p2f)
plsanter(alldata,c("GISTEMP","HADCRUT","NOAA","BERKELEY","CW"),tofr,plcol,trendlength,startyr,endyr,p2f)


#######################################################################################
#landscaping of slopes and p-values
analyses <- c(1,2)  #same as tofr, so 1=OLS and 2=OLS but continuous. 
pgm <- landscaping (alldata, datasets, timevars, plotvars, analyses, 0) #last arg: 0 = hindsight, 1 = historically conditioned
pgm <- landscaping (alldata, datasets, timevars, plotvars, analyses, 1) #last arg: 0 = hindsight, 1 = historically conditioned

#landscaping again but for Chow test or equivalent continuous
analyses <- c(1,2)  #same as tofr, so 1=OLS and 2=OLS but continuous. 
pgm <- chowlandscaping (alldata, datasets, timevars, plotvars, analyses, 0) #last arg: 0 = hindsight, 1 = historically conditioned
pgm <- chowlandscaping (alldata, datasets, timevars, plotvars, analyses, 1) #last arg: 0 = hindsight, 1 = historically conditioned


#######################################################################################
# histogram of slopes using latest available data set, comparing continuous to discontinuous
p2f <- 0
plothistslopes (alldata,c("CW"),cps,15,p2f)


#######################################################################################
#Monte Carlo analysis: first look at targeted snapshots
pauseyears <- c(2001,2016.999)
present <- pauseyears[2]  #normally the present is at end of pause, but we can also do hindsight

for (i in datasets) {
  #call Monte Carlo with a single data set, but pass all versions for historical conditioning
  mcresults[[i]] <- montecarlo(alldata[[i]],cps[i],pauseyears,present,nMC,tofr=1)  #note different start based on cp analysis for each data set
  print(c(i,round(unlist(mcresults[[i]]),3)))
}

#Monte Carlo: now do landscaping analysis
nMC4ls <- 1000
pausegrid <- matrix(NA,110,3)  #manually computed max, way over
k <- 0
earliestvantage <- 2007
present <- 2016
minyr <- 10
yrw<-c(minyr:19) 
for (p2 in c(earliestvantage:present)) {
  for (yback in yrw) {
    if (p2-yback+1 >= 1998) {
      k <- k+1
      pausegrid[k,] <- c(p2-yback+1,p2+.999,yback)
    }
  }
}
pausegrid <- na.omit(pausegrid) 
ncells <-dim(pausegrid)[1]
tofr <- 1   #1==broken, 2==continuous trends
histo <- 1  #1==historically conditioned, 0==hindsight
autocorrs <- 1  #0 = no autocorrelations, 1 = autocorrelations

gridds <-  c("GISTEMP","HADCRUT","NOAA")
gridmcresults <- vector("list",length(gridds))
names(gridmcresults) <- gridds
for (i in gridds) {                     # gridds do all the historically-available data sets (but not CW and BERKELEY)
  gridmcresults[[i]] <- vector("list",ncells)
  for (k in 1:ncells) {                 #run through grid
    print(k)                            #sign of life....
    pauseyears <- pausegrid[k,1:2]
    tpres      <- pauseyears[2]            #this ensures historical conditioning
    gridmcresults[[i]][[k]] <- montecarlo(alldata[[i]],cps[i],pauseyears,tpres,nMC4ls,tofr,histo,autocorrs) 
  }
}
save(gridmcresults,file=paste("gridmcresults_histo_",histo,"_tofr_",tofr,"_autocorrs_",autocorrs,".RData",sep=""))
#load(file=paste("gridmcresults_histo_",histo,"_tofr_",tofr,"_autocorrs_",autocorrs,".RData",sep=""))  #direct entry here without performing MC

#unpack results for landscaping plot
for (i in gridds) {  #do all the historically-available data sets gridds
  pdev <- magdev <- matrix(NA,10,10)
  for (k in 1:ncells) { #run through grid and fill matrix in landscape-conducive order
    row    <-  (floor(pausegrid[k,2]) - earliestvantage + 1)
    column <-  (floor(pausegrid[k,2]) - pausegrid[k,1] - 8)
    #for trend deviation: magdev[row,column] <- gridmcresults[[i]][[k]]$ptrend - gridmcresults[[i]][[k]]$ltrend
    #for p-values (50% white)
    magdev[row,column] <- gridmcresults[[i]][[k]]$pcpause
    pdev[row,column]   <- gridmcresults[[i]][[k]]$pcpause
  }

  #now do landscaping plot--ultimately this should be a function together with the others
  #for trend deviation: minmaxscale <- c(-.3, .3) #symmetrical to force white being at zero
  #for p-values (50%=white)
  minmaxscale <- c(0, 100)
  mags4p <- magdev  #create copy that can be altered for printing
  mags4p[mags4p < minmaxscale[1]] <- minmaxscale[1]
  mags4p[mags4p > minmaxscale[2]] <- minmaxscale[2]
  
  pvals4p <- pdev
  #pvals4p[pvals4p < 5] <- 5
  xysig <- which(pvals4p <= 5, arr.ind = TRUE)    #row & column of significance
  if (length(xysig) == 0) {xysig<-matrix(NA,1,2)} #avoid failure in printing points
  xy4num <- which(!is.na(pvals4p), arr.ind = TRUE) #row & column of printing numbers
  
  x11()  #if embedded in function, levelplot needs a print call
  print(
    levelplot(
      main=i,
      mags4p,
      row.values = c(earliestvantage:present),
      column.values = yrw,
      xlab = "Vantage year",
      ylab = "Years included",
      at = do.breaks(minmaxscale, 100),
      col.regions = colorRampPalette(c("red", "yellow", "yellow3", "white", "green", "green2", "green4"), 
                                     space = "rgb"),
      
      #add monte carlo numbers and 'significance' indicators
      panel = function(...) {
        panel.levelplot(...)
        panel.text(
          xy4num[, 1] + earliestvantage - 1,
          xy4num[, 2] + minyr - 1,
          round((pvals4p[xy4num]),0)
        )
        grid.points(
          xysig[, 1] + earliestvantage - 1,
          xysig[, 2] + minyr - 1,
          pch = 1,
          gp = gpar(col = "yellow",lwd=1.5),
          size = unit(1.9, "char"))
      }
    )
  )
  savePlot(filename = paste(outputdir,"/mcpvals_",i,"_histo_",histo,"_tofr_",tofr,"_autocorrs_",autocorrs,".pdf",sep=""),type = "pdf",
            device = dev.cur(),restoreConsole = TRUE)
  savePlot(filename = paste(outputdir,"/EPS/mcpvals_",i,"_histo_",histo,"_tofr_",tofr,"_autocorrs_",autocorrs,".eps",sep=""),type = "eps",
           device = dev.cur(),restoreConsole = TRUE)
} #end of data sets
