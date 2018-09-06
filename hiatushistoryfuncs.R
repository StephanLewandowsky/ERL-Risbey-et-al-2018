#functions for HiatusHistory -- commenced late February 2017
#  called from hiatushistoryTemps.R

#######################################################################################
#function to plot the trend since 1998 or whatever (in variable `boftrnd') up to various vantage points,
#  using data available at each particular time
#
plhistocond <-  function(alldata,ds4histo,tofr,tv, pv) {
  #set aside work space
  pt4p <- d4p <- vector("list", length(ds4histo))
  names(pt4p) <- names(d4p) <- ds4histo
  
  #set up frame for plot and determine file name if appropriate
  if (pv$p2f) {
    fn <- paste(outputdir,"/HHhistocond_",tofr,"_",tv$boftrnd[1],"_",
                paste(substr(ds4histo, 1, 1),sep = "",collapse = "_"),".pdf",sep = "")
    if (file.exists(fn)) {
      file.remove(fn)
    }
    pdf(file = fn,height = 7,width = 7)
  } else {
    x11(7, 7)
  }
  xvals <- seq(tv$boftrnd[1] + tv$mintrnd - 1, (tv$yrsused[2]), by = 1 / 12) #for X axis
  ylimits <- c(-0.05, .35)
  plot(xvals,rep(0, length(xvals)),
       type = "n",
       ylim = ylimits,
       xlab = "Vantage year", ylab = paste(ifelse(tofr==1,"Broken ","Continuous "),
                                           "GMST trend (K/decade) since ",tv$boftrnd[1],sep=""),
       las = 1,yaxt = "n",xaxt = "n")
  abline(h = 0, lty = "dashed", col = "dark gray")
  
  #for each data set, construct thick and thin lines in turn
  for (i in ds4histo) { #ds4histo is character array of names
    #grab a data set and work through all versions
    d4p[[i]] <-arrange4histocond(alldata[[i]], tv$boftrnd, tv$yrsused[2], tofr)
    #thin lines
    lapply(d4p[[i]], FUN = function(x) lines(xvals[c(1:length(x$slopes))],x$slopes,  
                                             col = pv$plthincol[i],lwd = 0.1,lty = "dotted") )
    #save distinct points to overlay later
    pt4p[[i]] <-lapply(d4p[[i]], FUN = function(x) c(xvals[length(x$slopes)], x$slopes[length(x$slopes)]))
  }
  for (i in ds4histo) {
    #single thick line plus version points for each data set
    thickline <- rep(NA, length(xvals))
    lapply(d4p[[i]][seq(length(d4p[[i]]), 2)], FUN = function(x) thickline[1:length(x$slopes)] <<- x$slopes)
    
    idx21st <- which(d4p[[i]][[1]]$vpts > tv$ds1strel[i])[1] #point to first release date
    thickline[1:(idx21st-1)] <- NA    #data predating the earliest release has to be masked
    #add labels for the two major datasets
    if (i == "HADCRUT") {
       tofl <- c(-.035,-.035,-.049,-.035,-.035)
       x4lines <- sapply(hadreldates,FUN=function(x) which(abs(x-xvals)<.05)[1])
       sapply(x4lines,FUN=function(x) lines(c(xvals[x],xvals[x]),c(-0.028,thickline[x]),col="gray"))
       for (k in 1:length(hadreldates)) {
         text(xvals[x4lines[k]],tofl[k],names(hadreldates)[k],col="dark gray")
       }
    }
    if (i == "GISTEMP") {
      tofl <- c(0.21, 0.21, 0.198, 0.21)+.015
      x4lines <- sapply(gissreldates,FUN=function(x) which(abs(x-xvals)<.05)[1])
      sapply(x4lines,FUN=function(x) lines(c(xvals[x],xvals[x]),c(.205,thickline[x]),col="gray"))
      for (k in 1:length(gissreldates)) {
        text(xvals[x4lines[k]],tofl[k],names(gissreldates)[k],col="dark gray")
      }
    }
    lines(xvals, thickline, col = plcol[i], lwd = 4)
    lapply(pt4p[[i]], FUN = function(x) points(x[1],x[2],
                                               pch = 21,col = "black",bg = pv$plcol[i],cex = 1.2))
  }
  axis(1,at = seq(from = xvals[1],to = floor(tv$yrsused[2]),by = 1) + .5,
       labels = seq(from = xvals[1],to = floor(tv$yrsused[2]),by = 1))
  axis(2,at = seq(from = ylimits[1], to = ylimits[2], by = .05),las = 1)
  legend(xvals[1] + 5,ylimits[2],
         ds4histo,
         col = pv$plcol[ds4histo],lty = "solid",lwd = 4,
         pch = 21,pt.bg = pv$plcol[ds4histo],pt.cex = 1.5)
  #sapply(hadreldates,FUN=function(x) lines(c(x,x),c(-0.02,.1),col="gray"))
  if (p2f) {
    dev.off()
  } else {savePlot(filename = paste(outputdir,"/HHhistocond_",tofr,"_",tv$boftrnd[1],"_",
                                    paste(substr(ds4histo, 1, 1),sep = "",collapse = "_"),".eps",sep = ""),
                   type = "eps",device = dev.cur(),restoreConsole = TRUE) 
  }

}

#function takes a list of variants of one data set and
#  returns decadal trend values for all variants
arrange4histocond <- function(thisd, boftrnd, lastt, tofr) {
  #embedded function to work on single variant of data set
  runtrends <- function(onevariant) {
    #embedded function to perform desired type of regression (continuous or broken)
    doreg <- function(vpts,x2) {
      slopes <- rep(NA,length(vpts))
      for (i in 1:length(vpts)) {
        subs4reg <- subset(x2, t <= vpts[i])  #only consider data preceding this vantage point
        if (tofr==1) { #broken trend
          slopes[i]<-lm(anom ~ t, data = subs4reg)$coefficients[2] * 10
        } else {       #continuous trend
          ctrd <- splice(subs4reg,boftrnd[1])
          slopes[i]<-lm(anom ~ presplice+postsplice, data=ctrd)$coefficients[3] * 10
        }
      }
      return(listN(slopes,vpts)) #associate date with each slope
    } #end embedded function doreg
    #subset this variant to beginning of AGW or start year for slopes, depending on tofr
    x2 <- dplyr::filter(onevariant, t >= boftrnd[tofr] & t <= mkyr(lastt))
    vpts <- x2$t[x2$t >= boftrnd[1] + mintrnd - 1 & x2$t <= mkyr(lastt)]
    if (length(vpts)>0) {
      return(doreg(vpts,x2))
    } else { 
      return(NULL)
    }
  } #end embedded function runtrends
  a4h <- lapply(thisd,FUN = function(x) runtrends(x))
  return(delete.NULLs(a4h))
}

#auxiliary function to center data set on 0 for continuous splicing
splice <- function (dat, yr) {
  dat$contint <- dat$t - yr
  dat$postsplice <- ifelse(dat$contint < 0, 0, dat$contint)
  dat$presplice <- ifelse(dat$contint >= 0, 0, dat$contint)
  return(dat)
}
#auxiliary function to strip empty elements of list
# from: https://stat.ethz.ch/pipermail/r-help/2006-August/111896.html
delete.NULLs  <-  function(x.list) {
  x.list[unlist(lapply(x.list, length) != 0)]
}



#######################################################################################
#landscaping function for any data set received as argument
#  pauselit not implemented in this version as yet
#
landscaping <- function (alldata, ds4histo, tv, pv, analyses, histo) {
  #Control the analysis in one place
  for (dsname in ds4histo) {
    for (ac in analyses) {
      ffname <- paste(outputdir,"/Lscape_",ac,"_",histo,"_",
                      paste(substr(dsname, 1, 1),sep = "",collapse = "_"),".eps",sep = "")
      if (file.exists(ffname)) {
        file.remove(ffname)
      }
      pffname <- paste(outputdir,"/pvLscape_",ac,"_",histo,"_",
                       paste(substr(dsname, 1, 1),sep = "",collapse = "_"),".eps",sep = "")
      if (file.exists(pffname)) {
        file.remove(pffname)
      }
      if (ac == 1) { #OLS -- same as tofr
        f2c <- getpv
        d2u <- alldata[[dsname]]
      }
      if (ac == 2) { #OLS but continuous -- same as tofr
        f2c <- getpvcont
        d2u <- alldata[[dsname]]
      } 

      #define computational and graphical parameters
      startyr <- floor(tv$yrslscap[1])
      present <- floor(tv$yrslscap[2])
      minyr <- 2
      earliestvantage <- startyr + tv$maxyrlscap
      yrw <- c(minyr:tv$maxyrlscap)
      mags <- bfvals <- pvals  <- matrix(0, length(c(present:earliestvantage)), tv$maxyrlscap - minyr + 1)
      
      
      #compute regressions from all vantage points
      k <- 0
      for (tpp in c(present:earliestvantage)) {
        k <- k + 1
        pp <- f2c (d2u, startyr, tpp, yrw, "anom", histo)  #f2c is getpv or getpvcont, depending on analysis
        pvals[k, ] <- pp$pv
        pvals[is.na(pvals)]<-0  #this is needed because there is no 'presplice' segment for longest trend and earliest vantage point
        mags[k, ] <- pp$mag
      }
      
      #2D color plots first
      minmaxscale <- c(-1., 1.) #symmetrical to force white being at zero
      mags4p <- mags  #create copy that can be altered for printing
      mags4p[mags4p < minmaxscale[1]] <- minmaxscale[1]
      mags4p[mags4p > minmaxscale[2]] <- minmaxscale[2]
      allrsig <- (dim(pvals)[2] - min(apply(pvals,1,FUN = function(x) which(x[c(length(x):1)] > .05)[1]), 
                                      na.rm = TRUE)) + minyr + 1.5
      print(allrsig)
      pvals4p <- pvals
      pvals4p[pvals4p < .05] <- .05
      xysig <- which(pvals4p[c(nrow(pvals4p):1), ] == .05, arr.ind = TRUE)
      
      x11()  #if embedded in function, levelplot needs a print call
      print(
        levelplot(
          main=dsname, #paste(dsname,ifelse(histo," historical"," hindsight"),
                       #   ifelse(ac==1," discontinuous"," continuous"),sep=""),
          mags4p[c(nrow(mags4p):1), ],
          row.values = c(earliestvantage:present),
          column.values = yrw + 1,
          xlab = "Vantage year",
          ylab = "Years included",
          at = do.breaks(minmaxscale, 100),
          col.regions =
            colorRampPalette(c(
              "purple", "blue", "white", "red", "brown"
            ), space = "rgb"),
          
          #add points
          panel = function(...) {
            panel.levelplot(...)
            #panel.2dsmoother(...,args=list(span=0.014),n=200)
            panel.abline(h = allrsig, lty = "dashed")
            
            grid.points(
              xysig[, 1] + earliestvantage - 1,
              xysig[, 2] + minyr,
              pch = 19,
              size = unit(.17, "char")
            )
            #                      grid.points(pauselit$yrto,
            #                                 pauselit$yrto-pauselit$yrfrom+1,pch=1,size=unit(1.17,"char"))
          }
        )
      )
      savePlot(filename = ffname,type = "eps",device = dev.cur(),restoreConsole = TRUE)
      
      #same for p-values
      x11()
      jet.colors <-colorRampPalette(c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000"))
      YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
      
      #custom scale for levelplot
      brks <-
        c(.15, .1, .09, .08, .07, .06, .05, .01,  0)  #breakpoints for plot
      pvalcolors <- c ("beige","gray85","gray80","gray75","gray70","gray65","gray60","coral1","coral2","coral3")
      #"darkgoldenrod1",
      #colorRampPalette(c("orange","green"))
      pvals4p <- pvals
      pvals4p[pvals4p > .1] <- .15
      print(
        levelplot(
          main=dsname, #paste(dsname,ifelse(histo," historical"," hindsight"),
                       #ifelse(ac==1," discontinuous"," continuous"),sep=""),
          pvals4p[c(nrow(pvals4p):1), ],
          row.values = c(earliestvantage:present),
          column.values = yrw + 1,
          xlab = "Vantage year",
          ylab = "Years included",
          at = brks,
          col.regions = pvalcolors[c(length(pvalcolors):1)],
          panel = function(...) {
            panel.levelplot(...)
            panel.abline(h = allrsig, lty = "dashed")
            
            if ((substr(dsname,1,1) == "G") && (ac == 1)) {
              panel.lines(
                c(earliestvantage - .5, (earliestvantage + tv$maxyrlscap - 15 + 1) - .5),
                c(15.5, tv$maxyrlscap + 1.5),
                lty = "solid",
                col.line = "black",
                lwd = 1
              )
              panel.text(earliestvantage + 4, tv$maxyrlscap - 3., "1970", srt =
                           45)
              panel.lines(
                c(earliestvantage - .5, (earliestvantage + tv$maxyrlscap - 10 + 1) - .5),
                c(10.5, tv$maxyrlscap + 1.5),
                lty = "solid",
                col.line = "black",
                lwd = 1
              )
              panel.text(earliestvantage + 6.56, tv$maxyrlscap - 5.26, "1975", srt =
                           45)
              panel.lines(
                c(earliestvantage - .5, (earliestvantage + tv$maxyrlscap - 20 + 1) - .5),
                c(20.5, tv$maxyrlscap + 1.5),
                lty = "solid",
                col.line = "black",
                lwd = 1
              )
              panel.text(earliestvantage + 1.4, tv$maxyrlscap - .5, "1965", srt =
                           45)
            }
          }
        )
      )
      savePlot(filename = pffname,type = "eps",device = dev.cur(),restoreConsole = TRUE)
    } #end of analysis loop
  } #end of data set loop
} #end of landscaping function

#compute p-values from arbitrary 'present' into the past (discontinuous trends)
getpv <- function (cw, startyr, present, yrw, dv, histo) {
  pv <- NULL
  sy <- NULL
  mag <- NULL
  #data set is determined by 'present' (unless hindsight)
  if (histo) {  
    ptr2cords <- 1
    while (present+.9 > max(cw[[ptr2cords]]$t)) { #time just shy of december guaranteed to trigger last one
      ptr2cords <- ptr2cords + 1  
    }
  } else {
    ptr2cords <- length(cw) 
  }
  #annualize the whole lot for efficiency
  cw4pres <- annualize(subset(cw[[ptr2cords]],t>=startyr))
  for (cutoff in yrw) {
    cwsub <- subset(cw4pres, t >= (present - cutoff) & t <= present)
    r4c <- lm(formula = cwsub[[dv]] ~ t, data = cwsub)
    pv  <- c(pv, summary(r4c)$coefficients[8])
    sy  <- c(sy, present - cutoff)
    mag <- c(mag, r4c$coefficients[2] * 10)
  }
  return(list(pv = pv, sy = sy, mag = mag))
}

#compute p-values from arbitrary 'present' into the past (*continuous* trends)
getpvcont <- function (cw, startyr, present, yrw, dv, histo) {
  pv <- NULL
  sy <- NULL
  mag <- NULL
  if (histo) {
    ptr2cords <- 1
    while (present+.9 > max(cw[[ptr2cords]]$t)) { #time just shy of december guaranteed to trigger last one
      ptr2cords <- ptr2cords + 1   
    }
  } else {  # take the latest data set for all but splice as needed
    ptr2cords <- length(cw) 
  }
  #annualize the whole lot for efficiency
  cw4pres <- annualize(subset(cw[[ptr2cords]],t>=startyr & t < mkyr(present)))  
  
  for (cutoff in yrw) {
    cwsub <- splice(cw4pres, present-cutoff) 
    r4c <- lm(formula = cwsub[[dv]] ~ presplice+postsplice, data = cwsub)
    pv  <- c(pv, summary(r4c)$coefficients[12]) 
    sy  <- c(sy, present - cutoff)
    mag <- c(mag, r4c$coefficients[3] * 10)
  }
  return(list(pv = pv, sy = sy, mag = mag))
}


#######################################################################################
#Chow landscaping function for any data set received as argument
#  pauselit not implemented in this version as yet
#
chowlandscaping <- function (alldata, ds4histo, tv, pv, analyses, histo) {
  #Control the analysis in one place
  for (dsname in ds4histo) {
    print(dsname)
    for (ac in analyses) {
      pffname <- paste(outputdir,"/chowpvLscape_",ac,"_",histo,"_",
                       paste(substr(dsname, 1, 1),sep = "",collapse = "_"),".pdf",sep = "")
      if (file.exists(pffname)) {
        file.remove(pffname)
      }
      if (ac == 1) { #do Chow conventional with broken trends
        f2c <- getpvChow
        d2u <- alldata[[dsname]]
      }
      if (ac == 2) { #do Chow with continuous trends
        f2c <- getpvChowcont
        d2u <- alldata[[dsname]]
      }
      
      #define computational and graphical parameters
      startyr <- cps[dsname]
      present <- floor(tv$yrslscap[2])
      minyr <- 2
      earliestvantage <- startyr + tv$maxyrlscap
      signofp <- pvals  <- matrix(NA, length(c(present:earliestvantage)), tv$maxyrlscap - minyr + 1)
      
      #compute regressions from all vantage points
      k <- 0
      for (tpp in c(present:earliestvantage)) {
        k <- k + 1
        yrw <- c(minyr:min(tv$maxyrlscap,tpp-startyr-15)) #15=min length of pre-break trend
        pp <- f2c (d2u, startyr, tpp, yrw, "anom", histo)  #f2c is one of the Chowpv functions
        pvals[k, 1:length(pp$pv)] <- pp$pv
        signofp[k, 1:length(pp$sofp)] <- pp$sofp
      }
      xysig <- which(signofp[c(nrow(signofp):1), ]>0, arr.ind = TRUE)
      if (length(xysig)==0) {xysig<-matrix(0,1,2)}
      
      #only p-values for Chow
      x11()
      jet.colors <-colorRampPalette(c("#00007F","blue","#007FFF","cyan","#7FFF7F","yellow","#FF7F00","red","#7F0000"))
      YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
      
      #custom scale for levelplot
      brks <-
        c(.15, .1, .09, .08, .07, .06, .05, .01,  0)  #breakpoints for plot
      pvalcolors <- c ("beige","gray85","gray80","gray75","gray70","gray65","gray60","coral1","coral2","coral3")
      #"darkgoldenrod1",
      #colorRampPalette(c("orange","green"))
      pvals4p <- pvals
      pvals4p[pvals4p > .1] <- .15
      print(
        levelplot(
          main=paste(dsname,ifelse(histo," historical"," hindsight"),
                     ifelse(ac==1," discontinuous"," continuous"),sep=""),
          pvals4p[c(nrow(pvals4p):1), ],
          row.values = c(earliestvantage:present),
          column.values = c(minyr:tv$maxyrlscap) +1,
          xlab = "Vantage year",
          ylab = "Years included",
          at = brks,
          col.regions = pvalcolors[c(length(pvalcolors):1)],
          
          #add points
          panel = function(...) {
            panel.levelplot(...)
             grid.points(
              xysig[, 1] + earliestvantage - 1,
              xysig[, 2] + minyr,
              pch = 3,
              size = unit(.5, "char")
            )
          }
        )
      )
      savePlot(filename = pffname,type = "pdf",device = dev.cur(),restoreConsole = TRUE)
    } #end of analysis loop
  } #end of data set loop
} #end of landscaping function


#compute p-values from arbitrary 'present' into the past, using Chow test to examine whether there is a breakpoint in the past
# this is discontinuous trends
getpvChow <- function (cw, startyr, present, yrw, dv, histo) {
  pv <- NULL
  sy <- NULL
  sofp <- NULL
  #data set is determined by 'present' (unless hindsight)
  if (histo) {  
    ptr2cords <- 1
    while (present+.9 > max(cw[[ptr2cords]]$t)) { #time just shy of december guaranteed to trigger last one
      ptr2cords <- ptr2cords + 1  
    }
  } else {
    ptr2cords <- length(cw) 
  }
  #annualize the whole lot for efficiency
  cw4pres <- annualize(subset(cw[[ptr2cords]],t>=startyr & t<= present+.999))
  for (cutoff in yrw) {
     cp <- which(cw4pres$t==(present - cutoff -1 ))
     
     #using package
     # chow<-sctest(cw4pres[[dv]] ~ t, data = cw4pres, type = "Chow", 
     #             point = cp) #point _before_ begin of second splice
     # pv  <- c(pv, as.numeric(chow$p.value))

     fac <- cw4pres$t > (present - cutoff -1)
     f0 <- lm(formula = cw4pres[[dv]] ~ t,     data = cw4pres)
     f1 <- lm(formula = cw4pres[[dv]] ~ fac/t, data = cw4pres)
     xx <- anova(f0, f1)
     pval <- xx[length(xx)][[1]][2]
     if (pval < .10) {
       #print(f1)
       print(c(present,cutoff,round(f0$coefficients[2],3),round(f1$coefficients,3),ifelse(f1$coefficients[4]>f1$coefficients[3],as.numeric(1),as.numeric(0)),
               pval,round(cw4pres$anom[cw4pres$t>(present - cutoff -1)],2)),digits=3)       
     }
     pv  <- c(pv, pval)
     sofp <- c(sofp, pval<.1 & as.numeric(f1$coefficients[4]>f1$coefficients[3]))
     sy  <- c(sy, present - cutoff)
    }
  return(listN(pv,sy,sofp))
}

#compute p-values from arbitrary 'present' into the past, using Chow test to examine whether there is a breakpoint in the past
# this is for continuous trends
getpvChowcont <- function (cw, startyr, present, yrw, dv, histo) {
  pv <- NULL
  sy <- NULL
  sofp <- NULL
  #data set is determined by 'present' (unless hindsight)
  if (histo) {  
    ptr2cords <- 1
    while (present+.9 > max(cw[[ptr2cords]]$t)) { #time just shy of december guaranteed to trigger last one
      ptr2cords <- ptr2cords + 1  
    }
  } else {
    ptr2cords <- length(cw) 
  }
  #annualize the whole lot for efficiency
  cw4pres <- annualize(subset(cw[[ptr2cords]],t>=startyr & t<= present+.999))
  for (cutoff in yrw) {
    f0 <- lm(anom ~ t, data = cw4pres)
    ctrd <- splice(cw4pres, (present - cutoff))
    f1 <-lm(anom ~ presplice+postsplice, data=ctrd)
    xx <- anova(f0, f1)
    pval <- xx[length(xx)][[1]][2]
    if (pval < .10) {
      print(f1)
      print(c(present,cutoff,round(f0$coefficients[2],3),round(f1$coefficients,3),ifelse(f1$coefficients[3]>f1$coefficients[2],as.numeric(1),as.numeric(0)),
              pval,round(cw4pres$anom[cw4pres$t>(present - cutoff -1)],2)),digits=3)       
    }
    pv  <- c(pv, pval)
    sofp <- c(sofp, pval<.1 & as.numeric(f1$coefficients[3]>f1$coefficients[2]))
    sy  <- c(sy, present - cutoff)
  }
  return(listN(pv,sy,sofp))
}


#some graphical cludge...
mylevelplot <- function(...) {
  panel.levelplot(...)
  panel.abline(h = 15.5, lty = "dashed")
}


#######################################################################################
#function takes data set and performs monte carlo analysis on the entire period, using
# trend of the prescribed duration.
montecarlo <- function(fullds,agwbegin,pauseyears,present,nMC,tofr,histo=1,autocorrs=0) {  #by default, do historical conditioning and no autocorrelations
  #embedded function to do all regs
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
  #take modern global warming period as data up to begin of pause and annualize
  # using appropriate historical conditioning
  present <- floor(ifelse(present>pauseyears[2],present,pauseyears[2])) #check for silly error and make sure it's annualized
  if (histo) {  
    ptr2cords <- 1
    while (present+.9 > max(fullds[[ptr2cords]]$t)) { #time just shy of december guaranteed to trigger last one
      ptr2cords <- ptr2cords + 1  
    }  
  } else {
    ptr2cords <- length(fullds) 
  }
  #having obtained correct historical data set, now filter into pause and other years
  ds  <- annualize(dplyr::filter(fullds[[ptr2cords]], t>=agwbegin & t < pauseyears[1])) 
  pds <- annualize(dplyr::filter(fullds[[ptr2cords]], t>=pauseyears[1] & t<= pauseyears[2]))
  
  pdur   <- floor(pauseyears[2]-pauseyears[1]+1)   #pause years are inclusive
  if (tofr==1) {      #broken trends
    ptrend <- lm(anom~t,data=pds)$coefficients[2]*10 #'pause' trend
  } else {           #continuous trends
    s4pt <- splice(rbind(ds,pds),pauseyears[1])
    ptrend <- lm(anom ~ presplice+postsplice, data=s4pt)$coefficients[3]*10
  }
  
  ldur   <- floor(pauseyears[1]-agwbegin)
  ltreg  <- lm(anom~t,data=ds)
  print(ncvTest(ltreg))
  ltrend <- ltreg$coefficients[2]*10   #long-term trend
  ltsdev <- sd(ltreg$residuals)        #standard deviation for Monte Carlo
  
  #depending on whether or not auto correlations are to be modeled....
  if (autocorrs==1) {
    bfparms <- getbestARMAreg(ltreg$residuals)
    print((bfparms))
    #convert observed variance of residuals into underlying uncorrelated-error variance
    #  by converting model to infinite MA process.
    #  https://ocw.mit.edu/courses/mathematics/18-s096-topics-in-mathematics-with-applications-in-finance-fall-2013/lecture-notes/MIT18_S096F13_lecnote8.pdf (slide 21)
    P <- ARMAtoMA(ar=bfparms[grep("ar",names(bfparms))],ma=bfparms[grep("ma",names(bfparms))],lag.max=100)
    sdepsilons <- ltsdev / sqrt(1+sum(P^2))
    
    #now create Monte Carlo samples ('realizations') for the best ARMA model of residuals
    modelspec  <- list(ar = bfparms[grep("ar",names(bfparms))], ma = bfparms[grep("ma",names(bfparms))])
    noise <- matrix(NA,nMC,ldur+pdur)
    for (i in 1:nMC) {
      noise[i,] <- (as.numeric(armaSim(model = modelspec, n = ldur+pdur, sd=sdepsilons)) + 
                           bfparms["intercept"]) 
    } 
  } else { #now monte carlo instantiations of the long-term series without autocorrelations
    noise <- matrix(rnorm(nMC*(ldur+pdur),sd=ltsdev),nMC,ldur+pdur)
  }
  synds <- noise + matrix(ltreg$coefficients[1]+ltreg$coefficients[2]*c(agwbegin:floor(pauseyears[2])),nMC,ldur+pdur,byrow=TRUE)
  
  ltptrend <- apply(synds,1,FUN=function(x) doallregs(x))
  pcpause <- sum(ltptrend)/nMC*100
  return(listN(pcpause,ptrend,ltrend,ltsdev))
}


#explore a number of candidate ARMA models for a given set of residuals
getbestARMAreg <- function(resids) {
  maxnmod <- 15  
  aics <- rep(NA,maxnmod)
  allcoefs <- vector("list",maxnmod)
  k<-0
  for (p in 0:1) {
    for (q in 0:1) {
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

#######################################################################################
#function to plot the trend since 1998 or whatever (in variable `boftrnd') up to various vantage points,
#  using data available at each particular time
# This version only plots thick lines for the three data sets highlighted in the Nature news piece
plhistocond4n <-  function(alldata,ds4histo,tofr,tv, pv) {
  #set aside work space
  pt4p <- d4p <- vector("list", length(ds4histo))
  plthincol <- c("red","light green","blue")
  names(plthincol) <- names(pt4p) <- names(d4p) <- ds4histo
  
  #set up frame for plot and determine file name if appropriate
  if (pv$p2f) {
    fn <- paste(outputdir,"/HHhistocond4N.pdf",sep = "")
    if (file.exists(fn)) {
      file.remove(fn)
    }
    pdf(file = fn,height = 7,width = 7)
  } else {
    x11(7, 7)
  }
  xvals <- seq(tv$boftrnd[1] + tv$mintrnd - 1, (tv$yrsused[2]), by = 1 / 12) #for X axis
  ylimits <- c(-0.05, .35)
  plot(xvals,rep(0, length(xvals)),
       type = "n",
       ylim = ylimits,
       xlab = "Vantage year", ylab = paste(ifelse(tofr==1,"Discontinuous ","Continuous "),
                                           "GMST trend (K/decade) since ",tv$boftrnd[1],sep=""),
       las = 1,yaxt = "n",xaxt = "n")
  abline(h = 0, lty = "dashed", col = "dark gray")
  
  #for each data set, plot slopes from 1998
  for (i in ds4histo) { #ds4histo is character array of names
    #grab a data set and work through all versions
    d4p[[i]] <-arrange4histocond(alldata[[i]], tv$boftrnd, tv$yrsused[2], tofr)
    #thin lines
    lapply(d4p[[i]], FUN = function(x) lines(xvals[c(1:length(x$slopes))],x$slopes,  
                                             col = plthincol[i],lwd = 3,lty = "solid") )
  }
  axis(1,at = seq(from = xvals[1],to = floor(tv$yrsused[2]),by = 1) + .5,
       labels = seq(from = xvals[1],to = floor(tv$yrsused[2]),by = 1))
  axis(2,at = seq(from = ylimits[1], to = ylimits[2], by = .05),las = 1)
  legend(xvals[1] + 5,ylimits[2],
         ds4histo,
         col = plthincol[ds4histo],lty = "solid",lwd = 3)
  #pch = 21,pt.bg = plthincol[ds4histo],pt.cex = 1.5)
  if (p2f) {
    dev.off()
  }
}


#######################################################################################
#function to plot histogram of slopes of certain duration 
plothistslopes <- function (alldata,datasets,beginyrs,trendlgth,p2f) {
  #another embedded function to do all regs
  doallregs <- function(thisds,tofr) {
    idxs <- c(1:(length(thisds$anom)-trendlgth+1))
    ts<- rep(NA,length(idxs))
    for (i in idxs) {
      if (tofr==1) {      #broken trends
        ts[i] <- lm(thisds$anom[i:(i+trendlgth-1)] ~ c(1:trendlgth))$coefficients[2]*10
      } else {            #continuous trends
        spliced <- splice(data.frame(anom=thisds$anom[1:(i+trendlgth-1)], 
                                     t=c(1:(i+trendlgth-1))), i)
        ts[i] <- lm(anom ~ presplice+postsplice, data=spliced)$coefficients[3]*10
      }
    }
    return(ts)
  }
  #loop over requested data sets
  for (i in datasets){
    if (p2f) {
      fn <- paste(outputdir,"/SlopeHisto_",i,".pdf",sep = "")
      if (file.exists(fn)) {
        file.remove(fn)
      }
      pdf(file = fn,height = 7,width = 7)
    } else {
      x11(7, 7)
    }
    tbpslopes <- doallregs(annualize(subset(alldata[[i]][[length(alldata[[i]])]],t>=beginyrs[i])),1)
    y1 <- hist(tbpslopes, plot=FALSE, breaks=8)
    #y1$density <- y1$density/sum(y1$density)
    plot( y1, col=rgb(0,0,1,1/4), freq=FALSE,  xlab="K/decade",ylab="Frequency",main=NULL,ylim=c(0,14),
            xlim=c(-.12,.42))
    
    tbpslopes2 <- doallregs(annualize(subset(alldata[[i]][[length(alldata[[i]])]],t>=beginyrs[i])),2)
    os <- hist(tbpslopes2, plot=FALSE)#,breaks=10)
    #os$density <- os$density/sum(os$density)
    plot( os, col=rgb(1,0,1,1/4), freq=FALSE,  add=TRUE)
    
    abline(v=0,lty="dashed",col="gray")
    if (p2f) {
      dev.off()
    }
  } #end of dataset loop
}


#######################################################################################
# read literature from disk and remove ineligible papers
getlit <- function() {
  fn <- "./Hiatus History Ari Literature/Hiatus_R17_mod_obs_added.csv"
  lit <- read.csv(fn,header=TRUE,stringsAsFactors = FALSE) %>% dplyr::filter(St.type %in% c("1":"8")) %>%
            dplyr::mutate(Start=as.numeric(Start)) %>% dplyr::mutate(End=as.numeric(End))
         
  lit$duration <- lit$End-lit$Start+1
  lit$ds <- cbind(ifelse(lit$D_G==1,"GISTEMP"), ifelse(lit$D_N==1,"NOAA"), 
              ifelse(lit$D_H3==1,"HADCRUT"), ifelse(lit$D_H4==1,"HADCRUT"),
              ifelse(grepl("best",lit$D_Other,ignore.case=TRUE ),"BERKELEY",NA),
              ifelse(grepl("cowtan",lit$D_Other,ignore.case=TRUE ),"CW",NA))
   tbk<-apply(lit$ds,1,FUN=function(x) sum(is.na(x))<6)
   return(subset(lit,tbk==1 & !is.na(Start) & !is.na(End) & Start < 2012))
}


#######################################################################################
#function to compute trend and reference distribution for a given paper in the corpus,
#  using data available at each particular time
#
getstat4paper <- function(alldata,ds4histo,startyr,endyr) {
  doallregs <- function(ds,duration) { #embedded function to do all regressions for reference slopes
      windows <- min(ds$t):(max(ds$t)-duration+1)
      refslopes <- rep(NA,length(windows))
      for (i in windows) {
          refslopes[i-min(windows)+1] <- lm(anom~t, data=dplyr::filter(ds,t>=i & t<=(i+duration-1)))$coefficients[2] * 10
      }
      return(refslopes)
  } #end embed to do all slopes
  k<-1
  while (max(alldata[[ds4histo]][[k]]$t) < endyr+.9) 
      {k<-k+1}
  ds <- annualize(dplyr::filter(alldata[[ds4histo]][[k]],t>=cps[ds4histo])) #modern AGW only
  ds4pause <- dplyr::filter(ds, t>=startyr & t<= endyr)
  pslope <- lm(anom~t,data=ds4pause)$coefficients[2] * 10
  refslopes <- doallregs(ds,endyr-startyr+1)
  return(listN(pslope,refslopes))
}



#######################################################################################
#function to plot trend of specified length at various vantage points,
#  using data available at each particular time
# This version only plots thick lines (no counterfactuals)
plsanter <-  function(alldata,ds4histo,tofr,plthincol,trendlength,startyr,endyr,p2f) {
  #set up frame for plot and determine file name if appropriate
  if (p2f) {
    fn <- paste(outputdir,"/santer_",trendlength,".pdf",sep = "")
    if (file.exists(fn)) {
      file.remove(fn)
    }
    pdf(file = fn,height = 7,width = 7)
  } else {
    x11(7, 7)
  }
  xvals <- seq(startyr + trendlength - 1, endyr) #for X axis
  ylimits <- c(-0.05, .4)
  plot(xvals,rep(0, length(xvals)),
       type = "n",
       ylim = ylimits,
       xlab = "Vantage year", ylab = paste(ifelse(tofr==1,"Discontinuous ","Continuous "),
                                           "GMST trend (K/decade)",sep=""),
       las = 1,yaxt = "n")
  abline(h = 0, lty = "dashed", col = "dark gray")
  
  #for each data set, plot trends of specified lengths at each vantage point
  for (i in ds4histo) { #ds4histo is character array of names
    #grab a data set and work through all versions
    k<-0
    ptrend <- lastptrend <- rep(NA, length(xvals))
    for (vp in xvals) { #go through all to be plotted years
      ptr2cords <- 1    #find correct data set available at the time in question
      while (vp+.9 > max(alldata[[i]][[ptr2cords]]$t)) { #time just shy of december guaranteed to trigger last one
        ptr2cords <- ptr2cords + 1  
      }  
      k <- k+1
      ptrend[k]     <- gettrends(vp,startyr,trendlength,alldata[[i]][[ptr2cords]])
      lastptrend[k] <- gettrends(vp,startyr,trendlength,alldata[[i]][[length(alldata[[i]])]])
    } #finish all to-be-plotted years
    lines(c((startyr+trendlength-1):endyr),lastptrend,col = plthincol[i],lwd = 0.5,lty="dashed")
    lines(c((startyr+trendlength-1):endyr),ptrend,col = plthincol[i],lwd = 3,lty="solid")
  }
  axis(2,at = seq(from = ylimits[1], to = ylimits[2], by = .05),las = 1)
  legend(startyr+trendlength+2,ylimits[2],
         ds4histo,
         col = plthincol[ds4histo],lty = "solid",lwd = 3)
  #pch = 21,pt.bg = plthincol[ds4histo],pt.cex = 1.5)
  if (p2f) {
    dev.off()
  }
}


#function produces 'pause' trends for Santer plots (for temp analysis as well as models)
gettrends <- function(vp,startyr,trendlength,dat) {
  #having obtained correct historical data set, now filter into pause vs other years
  if (tofr==1) {      #broken trends
    pds    <- annualize(dplyr::filter(dat, t>=(vp-trendlength+1) & t <=mkyr(vp)))
    ptrend <- lm(anom~t,data=pds)$coefficients[2]*10 #'pause' trend
  } else {           #continuous trends
    ds     <- annualize(dplyr::filter(dat, t>=startyr & t <=mkyr(vp)))
    s4pt   <- splice(ds,vp-trendlength+1)
    ptrend <- lm(anom ~ presplice+postsplice, data=s4pt)$coefficients[3]*10
  }
  return(ptrend)
}


