#######################################################################################
#   function illustrates the multiple testing problem
illustratemt <- function(mcresult,trendlength,thisdataset,tofr) {
  #embedded function to do all regressions for reference slopes
  doallregs <- function(ds,duration) {
    ds <- subset(ds,t>=cps[thisdataset]) #only consider modern global warming for comparison
    windows <- min(ds$t):(max(ds$t)-duration+1)
    refslopes <- rep(NA,length(windows))
    for (i in windows) {
      if (tofr==1) {      #broken trends
        refslopes[i-min(windows)+1] <- lm(anom~t, data=subset(ds,t>=i & t<=(i+duration-1)))$coefficients[2] * 10
      } else {            #continuous trends
        spliced <- splice(subset(ds,t<=(i+duration-1)),i)
        refslopes[i-min(windows)+1] <- lm(anom ~ presplice+postsplice, data=spliced)$coefficients[3]*10
      }
    }
    return(listN(refslopes,windows))
  } 
  #embedded function to compare slopes to ensemble
  slopecomp <- function(ensemble,rs,loslopeyr,duration) { 
    do1reg <- function(x) {
      if (tofr==1) {      #broken trends
        return(lm(anom~t, data=data.frame(subset(x,t>=loslopeyr & t<=(loslopeyr+duration-1) )))$coefficients[2] * 10 < critslope)
      } else {            #continuous trends
        spliced <- splice(data.frame(subset(x,t<=(loslopeyr+duration-1) )) ,loslopeyr)
        return(lm(anom ~ presplice+postsplice, data=spliced)$coefficients[3]*10 < critslope)
      }
    }
    print("slopecomp")
    critslope <- rs$refslopes[rs$windows==loslopeyr]
    comps <- lapply(ensemble, do1reg) 
    return (sum(unlist(comps)) / length(ensemble))
  }
  #embedded function to compare slopes to ensemble using correct (i.e., multiple-testing) comparison
  correctcomp <- function(ensemble,rs,loslopeyr,duration) { 
    print("correctcomp")
    critslope <- rs$refslopes[rs$windows==loslopeyr]
    comps <- rep(NA,length(ensemble))
    for (realn in 1:length(ensemble)) {
        forcomp <- doallregs(ensemble[[realn]],duration)
        comps[realn] <- ifelse(min(forcomp$refslopes)<critslope,1,0)
    }
    return (sum(comps)/length(ensemble))
  }
  #begin function illustratemt
  nrealiztns <- length(mcresult)   #number of realizations
  n <- nrealiztns/10               #use some number smaller than that to avoid sampling with replacement
  ltloslope <- ltranslope <- ltloslopecor <- rep(NA,n)
  ltlowishslope <- ltlowishslopecor <- matrix(NA,n,5)  #5 'lowish' slopes that cut off distribution in various places
  
  #pick random realizations for comparison to remainder of ensemble across replications
  sofreals <- sample(1:nrealiztns, n)
  for (i in 1:length(sofreals)) {
    #pick lowest slope and random slope for that realization, and also low-ish slopes
    rs<-doallregs(mcresult[[sofreals[i]]],trendlength)
    loslopeyr   <- rs$windows[which.min(rs$refslopes)]
    ranslopeyr  <- sample(rs$windows,1)
    lowishslopeyrs <- rs$windows[order(rs$refslopes)][c(1,2,as.numeric(floor(quantile(1:length(rs$windows),c(.1,.25,.5)))))]
    
    #comparisons to remainder of ensemble
    print(c(i,sofreals[i],length(rs$windows),date()))
    
    ltloslope[i]    <- slopecomp(mcresult[-sofreals[i]],rs,loslopeyr, trendlength)
    ltranslope[i]   <- slopecomp(mcresult[-sofreals[i]],rs,ranslopeyr,trendlength)
    ltloslopecor[i] <- correctcomp(mcresult[-sofreals[i]],rs,loslopeyr, trendlength)
    k<-0
    for (j in lowishslopeyrs) {
        k<-k+1
        ltlowishslope[i,k]    <- slopecomp(mcresult[-sofreals[i]],  rs,j, trendlength)
        ltlowishslopecor[i,k] <- correctcomp(mcresult[-sofreals[i]],rs,j, trendlength)
    }
  }
  return(listN(ltloslope,ltranslope,ltloslopecor,ltlowishslope,ltlowishslopecor))
}


#######################################################################################
#function to compute trend and reference distribution for a given paper in the corpus,
#  using the CMIP multi model mean or entire ensemble as reference
getmodstat4paper <- function(mmm,ds4histo,startyr,endyr) {
   doallregs <- function(ds,duration) { #embedded function to do all regressions for reference slopes
     verifyannual(ds) #assumes annualized data set. Check if true
     windows <- min(ds$t):(max(ds$t)-duration+1)
     refslopes <- rep(NA,length(windows))
     for (i in windows) {
       refslopes[i-min(windows)+1] <- lm(anom~t, data=subset(ds,t>=i & t<=(i+duration-1)))$coefficients[2] * 10
     }
     return(refslopes)
   } #end embed to do all slopes
  #begin by computing reference slope using the CMIP multi-model mean or ensemble
  if (length(mmm)==1) {
     refslopes <- lm(anom~t, data=subset(mmm[[1]],t>=startyr & t<= endyr))$coefficients[2] * 10 
  } else { #now deal with ensemble members
     refslopes <- as.numeric(unlist(lapply(mmm,FUN=function(x) lm(anom~t, data=subset(annualize(x),t>=startyr & t<= endyr))$coefficients[2] * 10)))
  }
  
  #now find appropriate historical data set and get single pause slope for this paper
  k<-1  
  while (max(alldata[[ds4histo]][[k]]$t) < (endyr+.9)) #last December will exceed and hence trigger exit from while
      {k<-k+1}
  ds4pause <- annualize(subset(alldata[[ds4histo]][[k]],t>=startyr & t<= mkyr(endyr))) #modern AGW only
  pslope <- lm(anom~t,data=ds4pause)$coefficients[2] * 10
  #compute position of pslope (used in paper) in distribution of trends of equal duration in the same data set since beginning of global warming
  ds4comp <- annualize(subset(alldata[[ds4histo]][[k]],t>=cps[ds4histo] & t<= mkyr(endyr)))
  allslopes <- doallregs(ds4comp,(endyr-startyr+1))                       
  q <- sum(allslopes<pslope)/length(allslopes)
  print(sum(pslope==allslopes))
  return(listN(pslope,refslopes,q))
}



#######################################################################################
#function calls all others needed for monte carlo simulation using 
#   the model-generated estimate of internal variability
mobsmontecarlowrapper <- function(allmodels,modelIDs,ds,endyr,nMC,tofm,modwt) {
  #obtain residuals between model and observations 
  mmmplusresid <- moresids(allmodels,
                           tofm,      #type of model (e.g., blended vs not)
                           modwt,     #equal weight (1) or average across all models and runs (0)
                           endyr,modelIDs,ds,
                           1)     #last arg determines historical conditioning of data (not models) up to 'endyr'
  resids <- mmmplusresid$resid
  globalmmm <<- data.frame(t=mmmplusresid$t,anom=mmmplusresid$mmm)
  #now get 'best' ARMA model of residuals and define that model
  bfparms <- getbestARMA(resids)
  #call Monte Carlo to generate realizations with the correct noise structure 
  return(mobsmontecarlo(bfparms,mmmplusresid,nMC))  #nMC = number of Monte Carlo realizations
}

#function performs actual monte carlo simulation
mobsmontecarlo <- function(bfparms,mmmplusresid,nMC) {
  #convert observed variance of residuals into underlying uncorrelated-error variance
  #  by converting model to infinite MA process.
  #  https://ocw.mit.edu/courses/mathematics/18-s096-topics-in-mathematics-with-applications-in-finance-fall-2013/lecture-notes/MIT18_S096F13_lecnote8.pdf (slide 21)
  P <- ARMAtoMA(ar=bfparms[grep("ar",names(bfparms))],ma=bfparms[grep("ma",names(bfparms))],lag.max=100)
  sdepsilons <- sd(mmmplusresid$resid) / sqrt(1+sum(P^2))
  
  #now create Monte Carlo samples ('realizations') for the best ARMA model of residuals
  modelspec  <- list(ar = bfparms[grep("ar",names(bfparms))], ma = bfparms[grep("ma",names(bfparms))])
  ninseries <- length(mmmplusresid$resid)

  # use parallel computing to speed up if machine has multiple cores
  myCluster <- makeCluster(max(1,detectCores()-1), type = "PSOCK") # use all but 1 cores in parallel
  registerDoParallel(myCluster)
  clusterExport(myCluster, c("modelspec", "ninseries", "sdepsilons", "bfparms"), envir=environment())
  realizatns <- foreach(1:nMC, .combine = rbind, .packages="fArma") %dopar% {
    as.numeric(armaSim(model = modelspec, n = ninseries, sd=sdepsilons)) + bfparms["intercept"]
  }
  stopCluster(myCluster)
  # this is returning the residual realization, without mmm being added
  return(list(apply(realizatns,1,FUN=function (x) data.frame(t=mmmplusresid$t,anom=x))))
}

#function determines which type of model to use for given data set and year
gettofm <- function(dsn,year,adjname) {
  hadyrs <- list(y=c(1879,2013,2014,2016),
                 tofms=c("masked_had3_tas_hemi","masked_had4_tas_hemi",paste("masked_had4_tas_hemi_",adjname,sep=""),
                         paste("masked_had4_blend_hemi_",adjname,sep="")))
  othyrs <- list(y=c(1879,2014,2016), 
                 tofms=c("global_tas",paste("global_tas_",adjname,sep=""),paste("global_blend_",adjname,sep=""))) 
  if (dsn=="HADCRUT") {
    return(hadyrs$tofms[max(which(hadyrs$y<=year))])
  } else {  
    return(othyrs$tofms[max(which(othyrs$y<=year))])
  } 
}

#function takes all realizations of model monte carlo and provides pause statistics
examineMC <- function(fullds,agwonset,mcrealizations,pauseyears,present,nMC,tofr) {
  #embedded function to do all regs on a monte carlo realization of models
  doallregs <- function(onerealiztn) {
    idxs <- c(1:(length(onerealiztn)-pdur+1))
    ts<- rep(NA,length(idxs))
    for (i in idxs) {
      if (tofr==1) {      #broken trends
        ts[i] <- lm(onerealiztn[i:(i+pdur-1)] ~ c(1:pdur))$coefficients[2]*10
        
      } else {            #continuous trends
        spliced <- splice(data.frame(anom=onerealiztn[1:(i+pdur-1)], 
                                     t=c(1:(i+pdur-1))), i)
        ts[i] <- lm(anom ~ presplice+postsplice, data=spliced)$coefficients[3]*10
      }
    }
    return(ifelse((min(ts) <= ptrend),1,0))
  }
  #embedded function to get pause trend (with mmm subtracted)
  getptrend <- function(ds,pauseyears) {
    if (tofr==1) {      #broken trends
      pds      <- annualize(dplyr::filter(ds, t>=pauseyears[1] & t<= mkyr(pauseyears[2]))) 
      pds$anom <- pds$anom - dplyr::filter(globalmmm, t>=pauseyears[1] & t<= mkyr(pauseyears[2]))$anom
      ptrend <- lm(anom~t,data=pds)$coefficients[2]*10 #'pause' trend
    } else {           #continuous trends
      s4pt   <- splice(annualize(dplyr::filter(ds, t<= mkyr(pauseyears[2]))),pauseyears[1])
      m4pt   <- dplyr::filter(globalmmm, t>=agwonset & t<= mkyr(pauseyears[2]))
      verifyannual(m4pt)
      # x11()
      # plot(m4pt$t,m4pt$anom)
      # lines(m4pt$t,m4pt$anom)
      # points(s4pt$t,s4pt$anom,pch=21,bg="red")
      # abline(v=pauseyears[1])
      # print (pauseyears)
      s4pt$anom <- s4pt$anom - m4pt$anom   #having spliced obs, subtract models from obs to work on deviations
      ptrend <- lm(anom ~ presplice+postsplice, data=s4pt)$coefficients[3]*10
    }
    return(ptrend)
  }
  # find historically appropriate data set
  ptr2cords <- 1
  while (max(fullds[[ptr2cords]]$t) < (present+.9)) { #time just shy of december guaranteed to trigger last one
    ptr2cords <- ptr2cords + 1  
  }  
  #having obtained correct historical data set, now filter out pause 
  pdur   <- floor(pauseyears[2]-pauseyears[1]+1)   #pause years are inclusive
  ptrend <- getptrend(subset(fullds[[ptr2cords]],t>=agwonset), pauseyears) 
  print(c(pauseyears,ptrend))

  #now compare observed pause trend to realizations, and compute percentage in realizations below pause trend
  #  first all regs from beginning of record onward ... 
  ltptrend <- lapply(mcrealizations[[1]],FUN=function(x) doallregs(x$anom[x$t<=mkyr(pauseyears[2])]))
  pcpause  <- sum(unlist(ltptrend))/nMC*100
  #  ... then all regs from beginning of AGW onward ...
  ltptrendagwonly <- lapply(mcrealizations[[1]],FUN=function(x) doallregs(x$anom[x$t>=agwonset & x$t<=mkyr(pauseyears[2])]))
  pcpauseagwonly  <- sum(unlist(ltptrendagwonly))/nMC*100
  #  ... then against 'pause' window only
  ltptrendponly <- lapply(mcrealizations[[1]],FUN=function(x) doallregs(x$anom[x$t>=pauseyears[1] & x$t<=mkyr(pauseyears[2])]))
  pcpauseponly  <- sum(unlist(ltptrendponly))/nMC*100
  return(listN(pcpause,pcpauseponly,pcpauseagwonly,ptrend))
}

#function that plots a monte carlo ensemble with CMIP multimodel mean and the data
plotmobsmc <- function(alldata,mmm,mcresult,tbpds,endyr,tofm) {
  if (p2f) {
    fn <- paste(outputdir,"/mobsmc_",substr(tbpds,1,3),"_",endyr,"_",tofm,".pdf",sep = "")
    if (file.exists(fn)) {
      file.remove(fn)
    }
    pdf(file = fn,height = 7,width = 7)
  } else {
    x11(7, 7)
  }
  xvals <- seq(yrsused[1],endyr)
  ylims2 <- c(-1.2,1.2)
  xlims2 <- c(min(xvals),max(xvals)+2)
  plot(xvals,rep(0,length(xvals)),
       type="n",xlim=xlims2,ylim=ylims2,xlab="Year",ylab="GMST Anomaly (K)", yaxt="n")
  
  #first obtain the envelope
  realizatns <- matrix(unlist(lapply(mcresult,FUN=function(x) return(x$anom))),nrow=nMC,byrow=TRUE)
  ensq25  <- apply(realizatns,2,quantile,0.025)
  ensq975 <- apply(realizatns,2,quantile,0.975)
  ensmin  <- apply(realizatns,2,min)
  ensmax  <- apply(realizatns,2,max)
  
  poly.y <- c(ensq975, ensq25[c(length(ensq25):1)])
  poly.x <- c(xvals, xvals[c(length(xvals):1)])
  polygon(poly.x,poly.y,col="thistle",border=NA)  
  lines(xvals,ensmax,lwd=1,col="red",lty="dashed")
  lines(xvals,mmm,lwd=2,col="red")
  lines(xvals,ensmin,lwd=1,col="red",lty="dashed")
  
  #now the observations
  for (i in tbpds) { 
    curannual <- annualize(getobs(alldata[[i]],endyr,1)) #last argument indicates historical conditioning to end year
    lines(xvals,curannual$anom,col=plcol[i],lwd=1.)
    points(xvals,curannual$anom,pch=plchar[i],bg=plcol[i],cex=0.5)
  }
  axis(2,at=round(seq(from=ylims2[1],to=ylims2[2],by=.2),1), las=1)
  legend(1900,ylims2[2]-.2,c(datasets[tbpds],"CMIP5 MMM"),col=c(plcol[tbpds],"red"),
         lty="solid",lwd=c(rep(1,length(datasets[tbpds])),2),pch=c(plchar[tbpds],NA),pt.bg=c(plcol[tbpds],NA))
  lines(c(max(xvals),max(xvals)),c(ylims2[1],.9),lty="dashed")
  text(max(xvals)-1,1,as.character(max(xvals)))
  if (p2f) {
    dev.off()
  }
}

#function to obtain data for comparison to models
#   argument histo determines if historically conditioned (1) or not (0) up to 'present'
getobs <- function(ds,present,histo) {
  if (histo) {  #obtain data in historically conditioned manner
    ptr2cords <- 1
    while (max(ds[[ptr2cords]]$t) < (present+.9)) { #time just shy of december guaranteed to trigger last one
      ptr2cords <- ptr2cords + 1  
    }
  } else {      #no historical conditioning
    ptr2cords <- length(ds) 
  }
  return(subset(ds[[ptr2cords]],t <= mkyr(present)))
}

#function to return multi-model mean for a given version of the model output
#  if equalweight == 1 then average computed within each model first,
#  if equalweight == 0 then average computed across all runs
getmmm <- function(allmodels,modelversion,present,modelIDs,equalweight) {
  getmodav <- function(thisone) {  #embedded function to average runs for one model
    partavg <- mv[modelIDs==thisone]
    partnum4avg <- matrix(unlist(lapply(partavg,FUN=function(x) return(x$anom))),nrow=sum(modelIDs==thisone),byrow=TRUE)
    return(apply(partnum4avg,2,mean))
  }
  #getmmm function commences here
  mv   <- allmodels[[modelversion]]
  nmod <- length(mv)
  if (equalweight) {  #average within models, then across models
    num4avg <- t(sapply(unique(modelIDs), getmodav))
  } else {            #average all runs irrespective of model 
    num4avg <- matrix(unlist(lapply(mv,FUN=function(x) return(x$anom))),nrow=nmod,byrow=TRUE)
  }
  return(subset(data.frame(t=mv[[1]]$t, anom=apply(num4avg,2,mean)), t<=mkyr(present)))
}

#function to compute discrepancy between multi-model mean and observations and return residuals
moresids <- function(allmodels,modelversion,equalwt,endyr,modelIDs,dataset,histo) {
  mmm <- annualize(getmmm (allmodels,modelversion,endyr,modelIDs,equalwt)) #last argument determines type of model average
  ooo <- annualize(getobs(dataset,endyr,histo))
  return(data.frame(t=ooo$t,resid=mmm$anom-ooo$anom,mmm=mmm$anom))
}

#explore a number of candidate ARMA models for a given set of residuals
getbestARMA <- function(resids) {
  maxnmod <- 15  
  aics <- rep(NA,maxnmod)
  allcoefs <- vector("list",maxnmod)
  k<-0
  for (p in 0:3) {
    for (q in 0:3) {
      if( p == 0 && q == 0 ) {next}    #skip impossible parameters
      k<-k+1
      formula <- as.formula( paste( sep="", "~ arma(", p, ",", q, ")" ) )
      mf      <- armaFit(formula,data=resids,include.mean=TRUE,method="ML") #mean is the intercept being returned 
      aics[k]      <-mf@fit$aic
      allcoefs[[k]]<-mf@fit$coef
    }
  }
  #return the 'best' ARMA model
  return(allcoefs[[which.min(aics)]])
}


####################################################################################
#provide annualized CMIP envelope for a model variant 
cmipenvelope <- function(modelvariant) {
  yrs <- annualize(modelvariant[[1]])$t
  
  modmeans <- matrix(0,length(yrs),length(modelvariant))
  for (i in c(1:length(modelvariant))) {
    annmod <- annualize(modelvariant[[i]])
    if (identical(annmod$t,yrs)) {
      modmeans[,i] <- annmod$anom
    } else {
      modmeans[,i] <- rep(NA,length(yrs))
    }
  }
  ensmean <- apply(modmeans,1,mean)
  enssd   <- apply(modmeans,1,sd)
  ensq25  <- apply(modmeans,1,quantile,0.025)
  ensq975 <- apply(modmeans,1,quantile,0.975)
  ensmin  <- apply(modmeans,1,min)
  ensmax  <- apply(modmeans,1,max)
  return(listN(ensmean,enssd,ensq25,ensq975,ensmin,ensmax))
}


#function that plots a particular type of CMIP runs together with the data
plotcmip <- function(whichversion,tbpds) {
  if (p2f) {
    fn <- paste(outputdir,"/cmip_",whichversion,"_bl",baselnpd[1],".pdf",sep = "")
    if (file.exists(fn)) {
      file.remove(fn)
    }
    pdf(file = fn,height = 7,width = 7)
  } else {
    x11(7, 7)
  }
  xvals <- seq(yrsused[1],floor(yrsused[2]))
  ylims2 <- c(-1.2,1.2)
  plot(xvals,rep(0,length(xvals)),
       type="n",ylim=ylims2,xlab="Year",ylab="GMST Anomaly (K)", yaxt="n")
  
  #first obtain the envelope
  cmip4plotting <- cmipenvelope(allmodels[[whichversion]])
  ub<- cmip4plotting$ensq975
  lb<- cmip4plotting$ensq25
  poly.y <- c(ub, lb[c(length(lb):1)])
  poly.x <- c(xvals, xvals[c(length(xvals):1)])
  polygon(poly.x,poly.y,col="thistle",border=NA)  
  lines(xvals,cmip4plotting$ensmax,lwd=1,col="black",lty="dashed")
  lines(xvals,cmip4plotting$ensmean,lwd=2.5,col="gray22")
  lines(xvals,cmip4plotting$ensmin,lwd=1,col="black",lty="dashed")
  
  #now the observations (latest available, no conditioning)
  for (i in tbpds) { 
    curannual <- annualize(alldata[[i]][[length(alldata[[i]])]])
    lines(xvals,curannual$anom,col=plcol[i],lwd=1.)
    points(xvals,curannual$anom,pch=plchar[i],bg=plcol[i],cex=0.5)
  }
  axis(2,at=round(seq(from=ylims2[1],to=ylims2[2],by=.2),1), las=1)
  legend(1900,ylims2[2]-.2,c(datasets[tbpds],"CMIP5 multi-model mean","Extreme projections","95% envelope"),
         col=c(plcol[tbpds],"gray22","black",NA),
         lty=c(rep("solid",length(tbpds)),"solid","dashed",NA),
         lwd=c(rep(1,length(tbpds)),2.5,1,NA),
         pch=c(plchar[tbpds],NA,NA,22),pt.bg=c(plcol[tbpds],NA,NA,"thistle"),
         pt.cex=c(rep(1,length(tbpds)),NA,NA,3))
  if (p2f) {
    dev.off()
  }
}


####################################################################################
#function to plot "santer plots" of models and data
plmodsanter <- function (alldata,datanames,allmodels,modelversion,adjname,
                         tofr,plcol,trendlength,startyr,endyr,p2f,mc=0) {
  #set up frame for plot and determine file name if appropriate
  if (p2f) {
    if (mc) {
      prefix <- "/MCsanter_"
    } else {
      prefix <- "/santer_"
    }
    fn <- paste(outputdir,prefix,trendlength,"_",endyr,"_",tofr,"_",
                paste(substr(datanames,1,3),collapse="_"),"_",
                modelversion,"_",adjname,".pdf",sep = "")
    if (file.exists(fn)) {
      file.remove(fn)
    }
    pdf(file = fn,height = 7,width = 7)
  } else {
    x11(7, 7)
  }
  modxvals <- xvals <- seq(startyr + trendlength - 1, endyr) #for X axis
  ylimits <- c(-0.2, .7)
  plot(xvals,rep(0, length(xvals)),
       type = "n",
       ylim = ylimits,
       xlab = "Vantage year", ylab = paste(ifelse(tofr==1,"Broken ","Continuous "),
                                           "GMST trend (K/decade)",sep=""),
       las = 1,yaxt = "n")
  
  #obtain envelope for CMIP (j indexes model run)
  if (length(grep("histo",modelversion))) {  #if historical, splice different versions
    if (length(grep("concat",modelversion))) {  #we want truly historical with different masks
      t2bused <- c("concatmasked_had34_tas_hemi", paste("masked_had4_tas_hemi_",adjname,sep=""),  
                   paste("masked_had4_blend_hemi_",adjname,sep=""))      
    } else {
      if (length(grep("had4",modelversion))) {  #we want masked historical for had4
        t2bused <- c("masked_had4_tas_hemi", paste("masked_had4_tas_hemi_",adjname,sep=""),  
                     paste("masked_had4_blend_hemi_",adjname,sep=""))
      } else {
        if (length(grep("had3",modelversion))) {  #we want masked historical for had3
          t2bused <- c("masked_had3_tas_hemi", paste("masked_had3_tas_hemi_",adjname,sep=""),  
                       paste("masked_had3_tas_hemi_",adjname,sep="")) #2nd & 3rd options will not be plotted ...
          modxvals <- seq(startyr + trendlength - 1, mkyr(had3endyr)) #...because had3 ends in 2012
        } else {  #historical unmasked
          t2bused <- c("global_tas",paste("global_tas_",adjname,sep=""),paste("global_blend_",adjname,sep="")) 
        }
      }
    }
    nmod <- length(allmodels[["global_blend"]])
    mptrend <- matrix(0,length(xvals),nmod)
    
    xv1 <- seq(startyr + trendlength - 1,floor(fadjdate)-1)
    for (j in 1:nmod) {
      modelvariant <- allmodels[[t2bused[1]]][[j]]
      mptrend [1:length(xv1), j] <- sapply(xv1,FUN=function(x) gettrends(x,startyr,trendlength,modelvariant)) 
    } 
    xv2 <- seq(floor(fadjdate),floor(blenddate))
    for (j in 1:nmod) {
      modelvariant <- allmodels[[t2bused[2]]][[j]]
      mptrend [(length(xv1)+1):(length(xv1)+length(xv2)), j] <- 
        sapply(xv2,FUN=function(x) gettrends(x,startyr,trendlength,modelvariant)) 
    } 
    xv3 <- seq(floor(blenddate)+1,endyr)
    for (j in 1:nmod) {
      modelvariant <- allmodels[[t2bused[3]]][[j]]
      mptrend [(length(xv1)+length(xv2)+1):(length(xv1)+length(xv2)+length(xv3)), j] <- 
        sapply(xv3,FUN=function(x) gettrends(x,startyr,trendlength,modelvariant)) 
    } 
  } else {    #models are not historical
    nmod <- length(allmodels[[modelversion]])
    mptrend <- matrix(0,length(xvals),nmod)
    for (j in 1:nmod) {
      print(j)
      modelvariant <- allmodels[[modelversion]][[j]]
      mptrend [, j] <- sapply(xvals,FUN=function(x) gettrends(x,startyr,trendlength,modelvariant)) 
    }    
  }
  ensmean <- apply(mptrend[1:length(modxvals),],1,mean)  #use modxvals not xvals in case had3
  enssd   <- apply(mptrend[1:length(modxvals),],1,sd)
  ensq25  <- apply(mptrend[1:length(modxvals),],1,quantile,0.025)
  ensq975 <- apply(mptrend[1:length(modxvals),],1,quantile,0.975)
  ensmin  <- apply(mptrend[1:length(modxvals),],1,min)
  ensmax  <- apply(mptrend[1:length(modxvals),],1,max)
  
  ub<- ensq975
  lb<- ensq25
  poly.y <- c(ub, lb[c(length(lb):1)])
  poly.x <- c(modxvals, modxvals[c(length(modxvals):1)])
  polygon(poly.x,poly.y,col="thistle",border=NA)  
  if (mc) {  #if it's monte carle realizations, print a few as spaghetti graph
    for (j in 1:min(30,nmod)) {
      lines(modxvals,mptrend[1:length(modxvals),j],lwd=0.5,col="pink")
      abline(v=floor(fadjdate),col="white",lwd=2)
      abline(v=floor(blenddate),col="white",lwd=2)
    } 
  }  else {
    lines(modxvals,ensmean,lwd=3,col="gray22") #only print mmm if the envelope is CMIP 
  }
  print(ensmean)
  lines(modxvals,ensmax,lwd=1,col="black",lty="dashed")
  lines(modxvals,ensmin,lwd=1,col="black",lty="dashed")
  abline(h = 0, lty = "dashed", col = "dark gray")
  
  #if historical, add separating lines to envelope to indicate updates of model output
  if (length(grep("histo",modelversion))) { 
      abline(v=floor(fadjdate)-.5,col="gray",lwd=2)
      abline(v=floor(blenddate)-.5,col="gray",lwd=2)
    }
  
  #for each data set, plot trends of specified lengths at each vantage point
  for (i in datanames) { #datanames is character array of names
    #grab a data set and work through all versions
    k<-0
    ptrend <- lastptrend <- rep(NA, length(xvals))
    for (vp in xvals) { #go through all to be plotted years
      ptr2cords <- 1    #find correct data set available at the time in question
      while (max(alldata[[i]][[ptr2cords]]$t) < (vp+.9)) { #time just shy of december guaranteed to trigger last one
        ptr2cords <- ptr2cords + 1  
      } 
      k <- k+1
      ptrend[k]     <- gettrends(vp,startyr,trendlength,alldata[[i]][[ptr2cords]])
      lastptrend[k] <- gettrends(vp,startyr,trendlength,alldata[[i]][[length(alldata[[i]])]])
    } #finish all to-be-plotted years
    
    #if historical, 
    if (length(grep("histo",modelversion))) { 
       lines(c((startyr+trendlength-1):endyr),ptrend,col = plcol[i],lwd = 1,lty="solid")
       points(c((startyr+trendlength-1):endyr),ptrend,pch=plchar[i],bg=plcol[i],cex=1.)
    } else {
      lines(c((startyr+trendlength-1):endyr),lastptrend,col = plcol[i],lwd = 1,lty="solid")
      points(c((startyr+trendlength-1):endyr),lastptrend,pch=plchar[i],bg=plcol[i],cex=1.)
    }
    print(ptrend)
    print(lastptrend)
  }
  axis(2,at = round(seq(from = ylimits[1], to = ylimits[2], by = .1),1),las = 1)
  if (mc) {
    legend(startyr+trendlength+2,ylimits[2],
           c(datanames,"Extreme realizations","95% envelope"),
           col=c(plcol[datanames],"black",NA),
           lty=c(rep("solid",length(datanames)),"dashed",NA),
           lwd=c(rep(1,length(datanames)),1,NA),
           pch=c(plchar[datanames],NA,22),pt.bg=c(plcol[datanames],NA,"thistle"),
           pt.cex=c(rep(1,length(datanames)),NA,3))
  } else {
    legend(startyr+trendlength+2,ylimits[2],
           c(datanames,"CMIP5 multi-model mean","Extreme projections","95% envelope"),
           col=c(plcol[datanames],"gray22","black",NA),
           lty=c(rep("solid",length(datanames)),"solid","dashed",NA),
           lwd=c(rep(1,length(datanames)),2.5,1,NA),
           pch=c(plchar[datanames],NA,NA,22),pt.bg=c(plcol[datanames],NA,NA,"thistle"),
           pt.cex=c(rep(1,length(datanames)),NA,NA,3))
  }
  if (p2f) {
    dev.off()
  }
}


#create data set for all models that concatenates the had3 and had4 masks into one
concatmasks <- function(allmodels,had3endyr) {
  nmod <- length(allmodels[["masked_had3_tas_hemi"]])
  conmask <- vector("list",nmod)
  for (i in 1:nmod) {
    conmask[[i]] <- rbind(subset(allmodels[["masked_had3_tas_hemi"]][[i]],t<=mkyr(had3endyr)),
                          subset(allmodels[["masked_had4_tas_hemi"]][[i]],t>=mkyr(had3endyr)))
  }
  return(conmask)
}


#adjust all models with updated forcings (Schmidt et al)
makeadj <- function(allmodels,forceadjf) {
  mk1adj <- function(amodv) {   #embedded function to adjust one type of models
    mk1mod <- function(amodel) { #embedded function for single model
      for (i in 1:length(adj$t)) {
        ptrs <- which(floor(amodel$t)==adj$t[i])
        amodel$anom[ptrs] <- amodel$anom[ptrs]+adj$forcedelta[i]
      }
      return(amodel)
    }
    return(lapply(amodv,FUN=function(x) mk1mod(x)))
  }
  adj <- getforceadj(forceadjf) #get the forces for adjustment
  return(lapply(allmodels,FUN=function(x) mk1adj(x)))
}

#function reads adjustment from file (Schmidt or Huber) and returns annualized data frame containing
#   vector of adjustments
getforceadj <- function(fn) {
  if (grepl("schmidt",fn)) {
    schmidt <- read.table(paste(modeldir,fn,sep="/"),skip=4,
                          row.names=NULL,header=FALSE,
                          colClasses = rep("numeric",12), 
                          col.names=c("t","c5lo","c5hi","giss","cw13","c5raw",
                                      "c5ensocor","c5ensoforcecor",
                                      "ensodelta","forcedelta","c5corlo","c5corhi"))
  }
  if (grepl("huber",fn)) {  #this could be either GISSvolc or GISSback, depending on file name
    schmidt <- read.table(paste(modeldir,fn,sep="/"),
                          row.names=NULL,header=TRUE,
                          colClasses = rep("numeric",2), 
                          col.names=c("t","forcedelta"))
    
  }
  tandf <- data.frame(t=schmidt$t,forcedelta=schmidt$forcedelta)
  for (ay in (max(schmidt$t)+1):yrsused[2]) {
    tandf <- rbind(tandf,data.frame(t=ay,forcedelta=schmidt$forcedelta[dim(schmidt)[1]]))
  }
  return(tandf)
}


