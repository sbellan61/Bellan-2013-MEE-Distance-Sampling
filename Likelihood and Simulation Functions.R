## Calculate likelihood for given parameters of detection function
## given carcass & road data.
library(distr)                          # for function truncate
library(Matrix)
allcues <- c("avian","mamm","cstate25" ,"cstate67","cstate8")
## file.folder <- "~/Documents/R files/dist samp/frequentist from cefe/"
file.folder <- ""

## ## prepare real data
## source("carc prepPC.r")
## start <- as.POSIXct("2010-02-01")
## end <- as.POSIXct("2010-05-31")
## carc.file <- "~/Documents/R files/dist samp/UCB Mortality Data 110702-cln111220.xls"
## dat <- carc.prep(carc.file = carc.file, period = TRUE, start = start,end = end, max.dist = 1.2, browse=F)
## carc <- dat[[1]]; head(carc,2)
## seff <- dat[[2]]
## seff <- seff[seff$Driver != "RZ",]
## carc <- carc[carc$driver != "RZ",]
## ## choose Roads
## load("roads.df.Rdata")
## roads.df.all <- roads.df                #for backup
## roads.df <- roads.df.all[roads.df.all$rdID %in% rrs,]
## roads.df$p <- roads.df$LENGTH / sum(roads.df$LENGTH)
## ## Choose time period
## d.min <- min(carc$nday)
## d.max <- max(carc$nday)
## dseq <- (d.min-btrack):d.max
## d.pdf <- data.frame(date = dseq, p = 1/length(dseq))
## ## select subsets
## seff <- seff[seff$nday >= (d.min-btrack) & seff$nday <= d.max,]
## rrs <- unique(seff$Road)
## seff <- seff[seff$Road %in% rrs,]
## ## correct time in of 2:15 instead of 1415
## seff[seff$tin < seff$tout,"tin"] <- as.POSIXct("2009-03-19 14:15:00")
## seff <- seff[order(seff$nday),]
## out.real <- list(carc = carc, seff = seff, d.min = d.min, d.max = d.max,
##             roads.df = roads.df, d.pdf = d.pdf, rrs = rrs)
## dr.out.master <- ddrr.fun(out.real, btrack = 5, browse = F)
## save(dr.out.master, file ="~/Documents/R files/dist samp/frequentist from cefe/all dr.out.Rdata")
## save(out.real, file = "realcarc120105.Rdata")
## save(seff, file = "seff.Rdata")
load("all dr.out.Rdata")                # load surveillance
load("roads.df.Rdata")                  # load road network GIS data
btrack <- 5                             # carcasses only visible for 5 days
load("realcarc120105.Rdata")            # load real carcass data
load("seff.Rdata")                      # load prepped surveillance effort data
roads.df.all <- roads.df                # for backup

## ## Fix time zone
## seff$tout <- format(seff$tout - 3600*9, usetz=F)
## seff$tin <- format(seff$tin - 3600*9, usetz=F)
## save(seff, file = "seff.Rdata")

######################################################################
##  Detection probability given distance and sighting cue
######################################################################
h.fun <- function(xx,                   # distance from road
                  cc,                   # sighting cue
                  browse = FALSE,
                  sig.v, sig.m, sig25, d.shape,
                  ...)
{
    if(browse) browser()
    allcues <- c("avian","mamm","cstate25" ,"cstate67","cstate8")
    allcues <- allcues[1:3]
    whicher <- allcues==cc
    denom <- c(sig.v,sig.m,sig25)[whicher]
    out <- exp( -(xx/denom)^d.shape )
    out
}
######################################################################
## Vectorize h.fun
######################################################################
vec.h <- function(xx,cc,
                  sig.v, sig.m, sig25, d.shape)
{
    Vectorize(h.fun, vectorize.args = c("xx","cc", "sig.v", "sig.m", "sig25", "d.shape"))
}
h.fun.v <- vec.h()

## Extract tripmat that shows # of times each trip type
## (tt=0,-1,-2,...) was driven in the last 5 days.
carcddrr <- function(out, btrack, verbose = F, browse = F) # d.min/max??
{
    carc <- out$carc
    seff <- out$seff
    if(browse) browser()
    tripmat <- matrix(0, nrow = nrow(carc), ncol = btrack+1)
    colnames(tripmat) <- 0:-btrack
    if(nrow(carc)>0)
    {
        for(ii in 1:nrow(carc))
        {
            on.road <- seff$Road == carc$roadwhich[ii]
            temp.seff <- seff[on.road,]
            temp.seff <- temp.seff[rev(order(temp.seff$nday)),]
            temp.ddiff <-  carc[ii,"nday"] - temp.seff$nday
            in.window <- temp.ddiff >= 0 & temp.ddiff <= btrack # changed second exp to <=, not sure if it is right
            temp.seff <- temp.seff[in.window,]
            temp.seff$dsf <- temp.seff$nday - carc$nday[ii]
            temp.tab <- xtabs(~temp.seff$dsf)
            ## feed the # of times each trip type (tt=0,...,-5) was driven into tripmat
            tripmat[ii, match(names(temp.tab), colnames(tripmat))] <- temp.tab
            ## calculate how many trips on day it was found were before the trip
            ## note this removes trip of detection too because time in > time out
            samedayaft <- sum(temp.seff$dsf == 0 & temp.seff$todin > carc$tod[ii])
            tripmat[ii,1] <- tripmat[ii,1] -  samedayaft
        }
    }
    return( tripmat)
}

## Construct matrix of ddrr.output that corresponds to each carcass
## (if carc dd=7, rr = 192, we want tripmat for dd=2:7 rr=192, and
## repeat down vector for each carc. gets fed into g.fun.  This
## handles the summation over d_i = (l_i-5, l_i) in the numerator

## NOTE: in the pr.functions I basically redid this again just copying
## the info from ddrr in a more straightforward way)
carcdrforward <- function(out, ddrr.output, btrack, verbose = F, browse = F)
{
    if(browse) browser()
    carc <- out$carc
    ddrr <- ddrr.output$ddrr
    tripmat <- ddrr.output$tripmat
    ind <- NULL
    for(ii in 1:nrow(carc))
    {
        for(tt in 0:btrack)
        {
            if(length(which(ddrr$dd == (carc$ll[ii]-tt) & ddrr$rr == carc$roadwhich[ii]))==0) browser()
            ind <- c(ind, which(ddrr$dd == (carc$ll[ii]-tt) & ddrr$rr == carc$roadwhich[ii]))

        }
    }
    tr.for <- tripmat[ind,]
    ddrr.for <- ddrr[ind,]
    tr.for <- tr.for[,1:max(ddrr.for$numtrips)]
    ddrr.for$carc <- rep(1:nrow(carc), each = btrack + 1)
    ddrr.for$yy <- rep(carc$roaddist, each = btrack + 1)
    ddrr.for$cc <- rep(carc$cue, each = btrack + 1)
    ddrr.for$ll <- rep(carc$ll, each = btrack + 1)
    return(list(ddrr.for = ddrr.for, tr.for = tr.for))
}

## Calculate numerator of likelihood. Two main steps, first must
## calculat [t_i][d_i], then must calculate g(y,c,d,r) for each
## potential d,t combination
g.fun <- function(out, cdr.output, foreff.output,
                  sig.v, sig.m, sig25, d.shape,
                  gamma.y = F, st.like = T,
                  temp.gamm,
                  max.dist,
                  cue.trips,
                  browse = F)
{
    if(browse) browser()
    carc <- out$carc
    seff <- out$seffs
    ## Calculate p det | each cue type at that distance
    pdet.allcues <- h.fun.v(xx = rep(carc$roaddist, each = length(allcues[1:3])), cc = rep(allcues[1:3], nrow(carc)),
                            sig.v = sig.v, sig.m = sig.m, sig25 = sig25, d.shape = d.shape)
    pdet.allcues <- matrix(pdet.allcues, nr = length(allcues[1:3]), nc = nrow(carc))
    ## calculate p det | given trip day (using camera trap p(cues | trip day)
    p.margbytrip <- cue.trips %*% pdet.allcues
    p.margbytrip.vec <- as.vector(p.margbytrip)
    nr <- btrack + 1
    nc <- nrow(carc)
    ctripmat <- cdr.output
    ######################################################################
    ## Part 1
    ######################################################################
    ## Calculate [t] given road effort: probability of detection on
    ## trip * p of und on all previous days carcass existed
    p.undmat <- matrix(1, nrow = nr, ncol = nc)
    ## for each possible age of the carcass (0,btrack)
    for(ii in 1:nr)
    {
        ## take the probability of undetection on all previous days
        for(jj in 1:ii)
        {
            p.undmat[ii,] <- p.undmat[ii,] * (1-p.margbytrip[jj,])^ctripmat[,ii-jj+1]
        }
    }
    ## probability of detection given the cue (independent of time bc cue known)
    pdet.givecue.vec <- h.fun.v(xx = carc$roaddist, cc = carc$cue,
                                sig.v = sig.v, sig.m = sig.m, sig25 = sig25, d.shape = d.shape)
    pdet.givecue <- matrix(rep(pdet.givecue.vec, each = nrow(p.undmat)), nrow = nr, ncol = nc)
    pdet.givecue.vec <- as.vector(pdet.givecue)
    ## prob of cue given t
    p.cue <- cue.trips[,match(carc$cue, allcues)]
    ## [d] for each day
    d.pdf <- out$d.pdf
    ind.mat <- matrix(rep(carc$nday, each = nr) - rep(0:btrack, nc), nrow = nr, ncol = nc)
    dps <- matrix(d.pdf[match(ind.mat, d.pdf$date), "p"], nr, nc)
    ## [d][det | cue seen, t, & undetection <t]
    pdets.give.td <- p.undmat * p.cue * dps # don't need pdet.givecue here because its the same for all trips
    ## now we take colsums to give the probability of detecting a
    ## carcass with this cue on day dd over all carcasses (where we
    ## ignore [rr] since we can just multiply by it at the end
    norm <- matrix(rep(colSums(pdets.give.td), each = nr), nr, nc)
    p.t <- pdets.give.td / norm         #note [d] is already in [t] now, so don't do it again later
    ######################################################################
    ## Part 2
    ######################################################################
    ## Now calculate the probability of detecting a carcass at that
    ## dist yy and with that cue cc for all the possible t's
    ## similar to the g.dk function below
    ddrr.for <- foreff.output$ddrr.for
    tr.for <- foreff.output$tr.for
    ## initialize gdrs (row is for each dd-rr combo, col is for subdivs)
    undet <- matrix(1, nr = nr, nc = nc)
    gdrs <- matrix(0, nr = nr, nc = nc)
    left.trunc  <- matrix(0, nr = nr, nc = nc)
    for(ii in 1:(ncol(tr.for)-1))         # -1 because last trip cannot be counted as an undetection (must be a detection)
    {
        is.trip <- ddrr.for$numtrips >= ii
        which.trips <- matrix(tr.for[,ii] + 1, nr, nc)
        is.trip.mat <- !is.na(which.trips)
                                        # + 1 is because trip 0's correspond to 1st row of
                                        # pmarg.bytrip, need to generalize if add tod
        ## need vector giving probability of sighting cue on that trip
        sel.mat <- cbind(as.vector(which.trips), match(ddrr.for$cc, colnames(cue.trips)))
        p.cue.mat <- matrix(cue.trips[sel.mat], nr, nc)
        if(ii > 1)
        {                           # need to add undetection prob for last trip
            last.which.trips <- tr.for[is.trip,ii-1] + 1
            ## because p.margybytip is a nc*nr vector now, need to
            ## make sure each carcass is matched with its own
            ## p.margbytrip and not the first one
            sel.mat.marg <- which(is.trip.mat, T)
            sel.mat.marg[,1] <- last.which.trips
            undet[is.trip.mat] <- undet[is.trip.mat]*(1 - p.margbytrip[sel.mat.marg])
        }
        gdrs[is.trip.mat] <- gdrs[is.trip] + undet[is.trip.mat] *  p.cue.mat[is.trip.mat]
        ## if(sum(gdrs>1)>1) browser()
        ## Need to left truncate effort outside time window. To do
        ## this have a separate matrix that will save a snapshot of
        ## gdrs at the last trip before d.min (the cumulative
        ## probability of having detected a carcass between d.min)
        ## which will then be subtracted from gdrs at the end.
        ## befores <- as.numeric(ddrr.for$dd) + tr.for[,ii] < out$d.min #trips being considered are before
        befores <- ind.mat + tr.for[,ii] < out$d.min
        if(ii<ncol(tr.for)) #if not at the max numtrips yet see if the next trip is after d.min
        {
            next.after <- ind.mat + tr.for[,ii+1] >= out$d.min
            ## next.after <- as.numeric(ddrr.for$dd) + tr.for[,ii + 1] >= out$d.min #& it is the last trip before d.min
        }else{ #if at max numtrips, then nextafter is TRUE because we want to subtract all previous prob
            next.after <- matrix(TRUE, nr ,nc) # this way, if all trips were before d.min, left.trunc becomes the full marg prob
        }
        next.after[is.na(next.after)] <- TRUE
        left.trunc[is.trip.mat & befores & next.after] <- gdrs[is.trip.mat & befores & next.after]
    }
    gdrs <- gdrs - left.trunc
    ## now multiply gdrs by the probability of detecting with that cue * [c|t]
    ## browser()
    gdrs <- gdrs * pdet.givecue.vec
    gdrs <- matrix(gdrs, nr = nr, nc = nc)
    ## multiply by [t][d] but only [d] if doing st.like
    pdets.and.td <- gdrs * p.t
    if(st.like) pdets.and.td <- pdets.and.td* dps
    ## sum across t
    pdets <- colSums(pdets.and.td)
    ## add [rr] and [yy]
    if(gamma.y)
    {                                   # use truncated gamma to get density
        y.pdf <- d(temp.gamm)(carc$roaddist)
    }else{
        y.pdf <- dunif(carc$roaddist,0,max.dist)
    }
    ## multiply by [r]'s (if st.like only)
    pdets.and.ry <- pdets * y.pdf
    if(st.like) pdets.and.ry <- pdets.and.ry * out$roads.df[match(carc$roadwhich, out$roads.df$rdID), "p"]
    ## multiply pdets.and.ry across carcasses in the likelihood function (take logsum really)
    return(pdets.and.ry)
    ## return(list(p.t = p.t, pdets.and.ry = pdets.and.ry)    )
}
## note gdrs>1 bc they are prob densities at yy

## Create effort matrices, ddrr giving road-day combos & number of
## subsequent trips w/in btrack days, and tripmat giving which trips
## were done for all d in [d.min-5, d.max]
ddrr.fun <- function(out, btrack, browse = F)
{
    seffset <- out$seff
    d.min <- out$d.min
    d.max <- out$d.max
    if(browse) browser()
    seffset$nday.fac <- factor(seffset$nday, levels = (d.min-btrack):d.max)
    stab <- xtabs(~Road + nday.fac, seffset)
    nc <- ncol(stab)
    nr <- nrow(stab)
    tripmat <-matrix(NA, nrow = nr*nc, ncol = 100) #100 is arbitrary assumed uperbound on # of trips in next btrack days
    ddrr <- data.frame(rr = rep(NA, nr*nc), dd = rep(NA, nr*nc), numtrips = rep(NA, nr*nc),
                       dd.p = rep(NA, nr*nc), rr.p = rep(NA, nr*nc))
    stepper <- 0
    for(r.ind in 1:nr)
    {
        for(d.ind in 1:nc)
        {
            stepper <- stepper + 1
            ddrr[stepper,"rr"] <- rownames(stab)[r.ind]
            ddrr[stepper,"dd"] <- colnames(stab)[d.ind]
            ddrr[stepper,"dd.p"] <- out$d.pdf[out$d.pdf$date == ddrr[stepper, "dd"], "p"]
            ddrr[stepper,"rr.p"] <- out$roads.df[out$roads.df$rdID  ==ddrr[stepper,"rr"],"p"]
            ## This ifelse deals with days at the end of time window so rep() statement works
                                        # don't need to worry about truncating effort at end if
                                        # data fed in is from d.min - btrack : d.max. Still need
                                        # to truncate effort at beginning (do it in gd.fun)
            if(d.ind + btrack > nc)
            {
                temp.dmax <- nc
                temp.btrack <- btrack - (d.ind + btrack - nc)
            }else{
                temp.dmax <- d.ind + btrack
                temp.btrack <- btrack
            }
            temp.trip.days <- rep(0:temp.btrack, stab[r.ind, d.ind:temp.dmax])
            ddrr[stepper,"numtrips"] <- length(temp.trip.days)
            if(ddrr[stepper,"numtrips"] > 0)
            {
                tripmat[stepper, 1:length(temp.trip.days)] <- temp.trip.days
            }
        }
    }
    maxnumtrips <- max(ddrr$numtrips)
    tripmat <- tripmat[,1:maxnumtrips]
    output <- list(ddrr = ddrr, tripmat = tripmat)
    return(output)
}


## ## Create effort matrices, ddrr giving road-day combos & number of
## ## subsequent trips w/in btrack days, and tripmat giving which trips
## ## were done for all d in [d.min-5, d.max]
## ddrr.fun <- function(out, btrack, browse = F)
## {
##     seffset <- out$seff
##     d.min <- out$d.min
##     d.max <- out$d.max
##     if(browse) browser()
##     seffset$nday.fac <- factor(seffset$nday, levels = (d.min-btrack):d.max)
##     stab <- xtabs(~Road + nday.fac, seffset)
##     nc <- ncol(stab)
##     nr <- nrow(stab)
##     tripmat <-matrix(NA, nrow = nr*nc, ncol = 100) #100 is arbitrary assumed uperbound on # of trips in next btrack days
##     ddrr <- data.frame(rr = rep(NA, nr*nc), dd = rep(NA, nr*nc), numtrips = rep(NA, nr*nc),
##                        dd.p = rep(NA, nr*nc), rr.p = rep(NA, nr*nc))
##     ddrr[,"rr"] <- rep(rownames(stab), each = nc)
##     ddrr[,"dd"] <- as.numeric(rep(colnames(stab), nr))
##     matcher <- match(ddrr[,"dd"], out$d.pdf$date)
##     ddrr[,"dd.p"] <- out$d.pdf[matcher, "p"]
##     matcher <- match(ddrr[,"rr"], out$roads.df$rdID)
##     ddrr[,"rr.p"] <- out$roads.df[matcher,"p"]
## ## This ifelse deals with days at the end of time window so rep() statement works
##                                         # don't need to worry about truncating effort at end if
##                                         # data fed in is from d.min - btrack : d.max. Still need
##                                         # to truncate effort at beginning (do it in gd.fun)
##     temp.dmax <- ddrr$dd + btrack
##     temp.btrack <- rep(btrack, nc*nr)
##     grtr <- (as.numeric(ddrr$dd) + btrack) > out$d.max
##     temp.dmax[grtr] <- out$d.max
##     temp.btrack[grtr] <- out$d.max - ddrr$dd[grtr]
##     temp.trip.days <- sapply
##             temp.trip.days <- rep(0:temp.btrack, stab[r.ind, d.ind:temp.dmax])
##             ddrr[stepper,"numtrips"] <- length(temp.trip.days)
##             if(ddrr[stepper,"numtrips"] > 0)
##             {
##                 tripmat[stepper, 1:length(temp.trip.days)] <- temp.trip.days
##             }
##         }
##     }
##     maxnumtrips <- max(ddrr$numtrips)
##     tripmat <- tripmat[,1:maxnumtrips]
##     output <- list(ddrr = ddrr, tripmat = tripmat)
##     return(output)
## }

## Select ddrr from master file
ddrr.select <- function(out, btrack,
                        master.file = "all dr.out.Rdata",
                        truncate.end = T,
                        browse = F)
{
    if(browse) browser()
    load(master.file)
    ddrr <- dr.out.master$ddrr
    tripmat<- dr.out.master$tripmat
    ## select for road-days in out
    sel <- ddrr[,"rr"] %in% out$roads.df$rdID
    sel <- sel & ddrr[,"dd"] %in% (out$d.min - btrack):out$d.max
    ddrr <- ddrr[sel,]
    tripmat <- tripmat[sel,]
    if(truncate.end)
    {
        for(ii in 0:btrack)
        {
            temp.ind <- ddrr[,"dd"] == out$d.max - ii
            grtr <- tripmat[temp.ind, ] > ii
            if(sum(temp.ind)>1)
            {
                less.trips <- rowSums(grtr, na.rm = T)
            }else{
                less.trips <- sum(grtr, na.rm = T)
            }
            ddrr[temp.ind, "numtrips"] <- ddrr[temp.ind, "numtrips"] - less.trips
            tripmat[temp.ind,][grtr] <- NA # remove all trips after d.max
        }
    }
    output <- list(ddrr = ddrr, tripmat = tripmat)
    return(output)
}


## ## Create ddrr matrix for pr.functions. Just calculate ddrr for all dd-rr within btrack days of a carc
## ddrr.carc <- function(out, btrack, browse = F)
## {
##     carc <- out$carc
##     seffset <- out$seff
##     d.min <- out$d.min
##     d.max <- out$d.max
##     if(browse) browser()
##     seffset$nday.fac <- factor(seffset$nday, levels = (d.min-btrack):d.max)
##     stab <- xtabs(~Road + nday.fac, seffset)
##     tripmat <-matrix(NA, nrow = (btrack+1)*nrow(carc), ncol = 100) #100 is arbitrary assumed uperbound on # of trips in next btrack days
##     ddrr <- data.frame(rr = rep(NA, (btrack+1)*nrow(carc)), dd = rep(NA, (btrack+1)*nrow(carc)), numtrips = rep(NA, (btrack+1)*nrow(carc)),
##                        dd.p = rep(NA, (btrack+1)*nrow(carc)), rr.p = rep(NA, (btrack+1)*nrow(carc)))
##     for(ii in 1:nrow(carc))
##     {
##         start.ind <- (ii-1)*(btrack+1)+1
##         ddrr[start.ind:(start.ind+btrack), "rr"] <- carc$roadwhich[ii]
##         ddrr[start.ind:(start.ind+btrack), "dd"] <- carc$ll[ii] - 0:btrack
##     }
##     ddrr.str <- paste(ddrr$rr, ddrr$dd) #for matching to n.in.eff below
##     unq.rows <- !duplicated(ddrr.str)
##     unq.ddrr <- ddrr[unq.rows,]
##     unq.tripmat <- tripmat[unq.rows,]
##     nc <- ncol(stab)
##     nr <- nrow(stab)
##     for(ii in 1:sum(unq.rows))
##     {
##         unq.ddrr[ii,"dd.p"] <- out$d.pdf[out$d.pdf$date == unq.ddrr[ii, "dd"], "p"]
##         unq.ddrr[ii,"rr.p"] <- out$roads.df[out$roads.df$rdID  ==unq.ddrr[ii,"rr"],"p"]
##         if(unq.ddrr[ii,"dd"] + btrack > out$d.max)
##         {
##             temp.dmax <- out$d.max
##             temp.btrack <- btrack - (unq.ddrr[ii,"dd"] + btrack - out$d.max)
##         }else{
##             temp.dmax <- unq.ddrr[ii,"dd"] + btrack
##             temp.btrack <- btrack
##         }
##         temp.trip.days <- rep(0:temp.btrack, stab[rownames(stab)==unq.ddrr[ii,"rr"], colnames(stab) %in% unq.ddrr[ii,"dd"]:temp.dmax])
##         unq.ddrr[ii,"numtrips"] <- length(temp.trip.days)
##         if(unq.ddrr[ii,"numtrips"] > 0)
##         {
##             unq.tripmat[ii, 1:length(temp.trip.days)] <- temp.trip.days
##         }
##     }
##     maxnumtrips <- max(unq.ddrr$numtrips)
##     unq.tripmat <- unq.tripmat[,1:maxnumtrips]
##     matcher <- match(ddrr.str, ddrr.str[unq.rows])
##     ddrr <- unq.ddrr[matcher,]
##     tripmat <- unq.tripmat[matcher,]
##     output <- list(ddrr = ddrr, tripmat = tripmat)
##     return(output)
## }

## Calculate the probability of detecting each carcass that was
## detected given only r and d, for use in the likelihood that does
## not account for spatiotemporal patterns (i.e. no
## [r,d]). Denominator for non-sptemp likelihood
pr.fun <- function(out, cdr.output, ddrr.output,
                   sig.v, sig.m, sig25, d.shape,
                   yys = yy.seq,       # quadrature subdivisions
                   gamma.y = F,
                   temp.gamm, max.dist, cue.trips,
                   browse = F)
{
    if(browse) browser()
    carc <- out$carc
    ## Calculate p det | each cue type at quadtrature distances
    pdet.allcues <- h.fun.v(xx = rep(yys, each = length(allcues[1:3])), cc = rep(allcues[1:3], length(yys)),
                            sig.v = sig.v, sig.m = sig.m, sig25 = sig25, d.shape = d.shape)
    pdet.allcues <- matrix(pdet.allcues, nr = length(allcues[1:3]), nc = length(yys))
    ## calculate p det | given trip day (using camera trap p(cues | trip day)
    p.margbytrip <- cue.trips %*% pdet.allcues
    ## p.margbytrip.vec <- as.vector(p.margbytrip)
    nr <- btrack + 1
    nc <- nrow(carc)
    ctripmat <- cdr.output
######################################################################
    ## Part 1
######################################################################
    ## Calculate [t] given road effort: probability of detection on
    ## trip * p of und on all previous days carcass existed
    p.undmat <- array(1, dim =c(nr, nc, length(yys)))
    ## for each possible age of the carcass (0,btrack)
    for(ii in 1:nr)
    {
        ## take the probability of undetection on all previous days
        for(jj in 1:ii)
        {
            is.trip <- which(ctripmat[,ii-jj+1]>0) #which carcasses had trip on that day (for speed)
            temp.trips <- ctripmat[rep(is.trip, length(yys)),ii-jj+1]
            p.undmat[ii,is.trip,] <- p.undmat[ii,is.trip,] * (1-p.margbytrip[rep(jj,length(is.trip)),]) ^temp.trips
        }
    }
    ## Calculate probability of detection on that day
    ## probability of detection given the cue (independent of time bc cue known)
    ## pdet.givecue <- pdet.allcues[match(carc$cue, allcues[1:3]),]
    ## prob of cue given t
    p.cue <- cue.trips[,match(carc$cue, allcues[1:3])]
    ## [d] for each day
    d.pdf <- out$d.pdf
    ind.mat <- matrix(rep(carc$nday, each = nr) - rep(0:btrack, nc), nrow = nr, ncol = nc)
    dps <- matrix(d.pdf[match(ind.mat, d.pdf$date), "p"], nr, nc)
    ## [d][det | cue seen, t, & undetection <t]
                                        # don't need pdet.givecue here because its the same for all trips
                                        # and so it divides out when we normalize below (basicaly [t] is
                                        # just the prob of not observing previously for each day
    pdets.give.td <- abind(lapply(1:length(yys), function(kkk) p.undmat[,,kkk] * p.cue * dps), along = 3)
    ## now we take colsums to give the probability of detecting a
    ## carcass with this cue on day dd over all carcasses (where we
    ## ignore [rr] since we can just multiply by it at the end
    norm <- array(rep(colSums(pdets.give.td), each = nr), dim = c(nr, nc, length(yys)))
    p.t <- pdets.give.td / norm         #note [d] is already in [t] now, so don't do it again later
######################################################################
    ## Part 2
######################################################################
    ## Now calculate the probability of detecting a carcass at each
    ## yys dist for each trip and for each possible days old of the
    ## carcass
    ## Select the future btrack effort from the ddrr matrix that
    ## matches to the btrack days behind each carcass
    ddrr <- ddrr.output$ddrr
    tripmat <- ddrr.output$tripmat
    ddrr.str <- paste(ddrr$rr, ddrr$dd) #for matching to n.in.eff below
    ## match ddrr.for and tr.for rows to corresponding rr-dd from carc
    carc.str <- paste(rep(carc$roadwhich, each = nr), as.vector(ind.mat))
    matcher <- match(carc.str, ddrr.str)
    ddrr <- ddrr[matcher,]
    tripmat <- tripmat[matcher,]
    ## initialize gdrs (row is for each dd-rr combo, col is for subdivs)
    gdrs <- matrix(0, nr = nrow(ddrr), ncol = length(yys))
    left.trunc <- matrix(0, nr = nrow(ddrr), ncol = length(yys))
    for(ii in 1:ncol(tripmat))         # -1 because last trip cannot be counted as an undetection (must be a detection)
    {
        ## if there are still trips for that dd-rr combo
        is.trip <- ddrr$numtrips >= ii
        ## which trips are they (classified by t_i) but we add + 1 to
        ## allow us to index the t=0 trips as the first row of the
        ## p.margybytrip matrix
        which.trips <- tripmat[is.trip,ii] + 1
        ## the probability of detection is the accumulated probability
        ## of detection + (1-the accumulated probability of
        ## detection)*(the marginal probability of detection on the
        ## next trip classified by t)
        gdrs[is.trip,] <- gdrs[is.trip,] + (1-gdrs[is.trip,]) * p.margbytrip[which.trips,]
        ## Need to left truncate effort outside time window. To do
        ## this have a separate matrix that will save a snapshot of
        ## gdrs at the last trip before d.min (the cumulative
        ## probability of having detected a carcass between d.min)
        ## which will then be subtracted from gdrs at the end.
        befores <- as.numeric(ddrr$dd[is.trip]) + tripmat[is.trip,ii] < out$d.min #trips being considered are before
        if(ii<ncol(tripmat)) #if not at the max numtrips yet see if the next trip is after d.min
        {
            next.after <- as.numeric(ddrr$dd[is.trip]) + tripmat[is.trip,ii + 1] >= out$d.min #& it is the last trip before d.min
        }else{ #if at max numtrips, then nextafter is TRUE because we want to subtract all previous prob
            next.after <- TRUE # this way, if all trips were before d.min, left.trunc becomes the full marg prob
        }
        next.after[is.na(next.after)] <- TRUE
        if(sum(is.trip)==1)             #separate if statement if it is a vector not a matrix
        {
            left.trunc[is.trip,][befores & next.after] <- gdrs[is.trip,][befores & next.after]
        }else{
            left.trunc[is.trip,][befores & next.after,] <- gdrs[is.trip,][befores & next.after,]
        }
    }
    gdrs <- gdrs - left.trunc
    ## re-order so it matches the p.t dimensions
    gdrs.arr <- array(gdrs, dim = c(nr, nc, length(yys)))
    ## multiply by p.t to collapse across the nr (potential days since death) dimension
    gdrs.t <- gdrs.arr * p.t
    ## multiply by y.pdf to collapse across the length(yys) dimension
    if(gamma.y)
    {
        ## y.pdf <- d(temp.gamm)(yy.seq)
        y.pdf <- p(temp.gamm)(yy.seq + int/2) - p(temp.gamm)(yy.seq - int/2)
    }else{
        y.pdf <- punif(yy.seq + int/2, 0, max.dist) - punif(yy.seq - int/2, 0, max.dist)
        ## y.pdf <- dunif(yy.seq, 0, max.dist)
    }
    ## now have marginal probability of detecting the carcass that was detected
    carc.probs <- colSums(gdrs.t) %*% y.pdf
    return(carc.probs)
}


## Calculate the probability of detecting each carcass that was
## detected given rr, dd, AND yy. For use in Horvitz-Thompson
## estimator without spatiotemporal patterns (i.e. no [r,d]).
pr.fun.yspec <- function(out, cdr.output, ddrr.output,
                         sig.v, sig.m, sig25, d.shape,
                         yys = yy.seq,       # quadrature subdivisions
                         gamma.y = F,
                         temp.gamm = NULL, cue.trips,
                         browse = F, late.br = F)
{
    if(browse) browser()
    carc <- out$carc
    ## Calculate p det | each cue type at quadtrature distances
    pdet.allcues <- h.fun.v(xx = rep(carc$roaddist, each = length(allcues[1:3])), cc = rep(allcues[1:3], nrow(carc)),
                            sig.v = sig.v, sig.m = sig.m, sig25 = sig25, d.shape = d.shape)
    pdet.allcues <- matrix(pdet.allcues, nr = length(allcues[1:3]), nc = nrow(carc))
    ## calculate p det | given trip day (using camera trap p(cues | trip day)
    p.margbytrip <- cue.trips %*% pdet.allcues
    ## p.margbytrip.vec <- as.vector(p.margbytrip)
    nr <- btrack + 1
    nc <- nrow(carc)
    ctripmat <- cdr.output
######################################################################
    ## Part 1
######################################################################
    ## Calculate [t] given road effort: probability of detection on
    ## trip * p of und on all previous days carcass existed
    p.undmat <- matrix(1, nrow = nr, ncol = nc)
    ## for each possible age of the carcass (0,btrack)
    for(ii in 1:nr)
    {
        ## take the probability of undetection on all previous days
        for(jj in 1:ii)
        {
            p.undmat[ii,] <- p.undmat[ii,] * (1-p.margbytrip[jj,])^ctripmat[,ii-jj+1]
        }
    }
    ## prob of cue given t
    p.cue <- cue.trips[,match(carc$cue, allcues[1:3])]
    ## [d] for each day
    d.pdf <- out$d.pdf
    ind.mat <- matrix(rep(carc$nday, each = nr) - rep(0:btrack, nc), nrow = nr, ncol = nc)
    dps <- matrix(d.pdf[match(ind.mat, d.pdf$date), "p"], nr, nc)
    ## [d][det | cue seen, t, & undetection <t]
    pdets.give.td <- p.undmat * p.cue * dps # don't need pdet.givecue here because its the same for all trips
    ## now we take colsums to give the probability of detecting a
    ## carcass with this cue on day dd over all carcasses (where we
    ## ignore [rr] since we can just multiply by it at the end
    norm <- matrix(rep(colSums(pdets.give.td), each = nr), nr, nc)
    p.t <- pdets.give.td / norm         #note [d] is already in [t] now, so don't do it again later
######################################################################
    ## Part 2
######################################################################
    ## Now calculate the probability of detecting a carcass given its
    ##  rr, dd, yy for each trip and for each possible days old of the
    ##  carcass. Select the future btrack effort from the ddrr matrix
    ##  that matches to the btrack days behind each carcass
    ddrr <- ddrr.output$ddrr
    tripmat <- ddrr.output$tripmat
    ddrr.str <- paste(ddrr$rr, ddrr$dd) #for matching to n.in.eff below
    ## match ddrr.for and tr.for rows to corresponding rr-dd from carc
    carc.str <- paste(rep(carc$roadwhich, each = nr), as.vector(ind.mat))
    matcher <- match(carc.str, ddrr.str)
    ddrr <- ddrr[matcher,]
    tripmat <- tripmat[matcher,]
    ## initialize gdrs (row is for each dd-rr combo, col is for subdivs)
    gdrs <- rep(0, nrow(ddrr))
    left.trunc <- rep(0, nrow(ddrr))
    for(ii in 1:ncol(tripmat))         # -1 because last trip cannot be counted as an undetection (must be a detection)
    {
        ## if there are still trips for that dd-rr combo
        is.trip <- ddrr$numtrips >= ii
        ## which trips are they (classified by t_i) but we add + 1 to
        ## allow us to index the t=0 trips as the first row of the
        ## p.margybytrip matrix
        which.trips <- tripmat[is.trip,ii] + 1
        ## the probability of detection is the accumulated probability
        ## of detection + (1-the accumulated probability of
        ## detection)*(the marginal probability of detection on the
        ## next trip classified by t)
        trip.carc.ind <- cbind(which.trips, rep(1:nc, each = nr)[is.trip])
        gdrs[is.trip] <- gdrs[is.trip] + (1-gdrs[is.trip]) * p.margbytrip[trip.carc.ind]
        ## Need to left truncate effort outside time window. To do
        ## this have a separate matrix that will save a snapshot of
        ## gdrs at the last trip before d.min (the cumulative
        ## probability of having detected a carcass between d.min)
        ## which will then be subtracted from gdrs at the end.
        befores <- as.numeric(ddrr$dd[is.trip]) + tripmat[is.trip,ii] < out$d.min #trips being considered are before
        if(ii<ncol(tripmat)) #if not at the max numtrips yet see if the next trip is after d.min
        {
            next.after <- as.numeric(ddrr$dd[is.trip]) + tripmat[is.trip,ii + 1] >= out$d.min #& it is the last trip before d.min
        }else{ #if at max numtrips, then nextafter is TRUE because we want to subtract all previous prob
            next.after <- TRUE # this way, if all trips were before d.min, left.trunc becomes the full marg prob
        }
        next.after[is.na(next.after)] <- TRUE
        left.trunc[is.trip][befores & next.after] <- gdrs[is.trip][befores & next.after]
    }
    ## subtract off probability accumulated before d.min
    gdrs <- gdrs - left.trunc
    ## re-order so it matches the p.t dimensions
    gdrs.arr <- matrix(gdrs, nr, nc)
    ## multiply by p.t to collapse across the nr (potential days since death) dimension
    gdrs.t <- 1 /gdrs.arr * p.t
    gdrs.t[is.na(gdrs.t)] <- 0
    ## now have probability of detecting observed carcasses given dd, rr, yy
    carc.prob <- 1 / colSums(gdrs.t)         # effective carcasses
    if(late.br) browser()
    return(carc.prob)
                                        #    sum(p.t / grds.arr )
}

## Calculate denominator
g.dk <- function(ddrr.output,                 # output of ddrr
                 out,                         #
                 only.obs = FALSE, # modify to calculate p of det given it could have been detected (i.e. within range of seff)
                 yys = yy.seq,       # quadrature subdivisions
                 sig.v, sig.m, sig25, d.shape,
                 gamma.y = F, temp.gamm, max.dist, cue.trips,
                 verbose = F, browse = F)
{
    if(browse) browser()
    if(verbose) print(list(sig.v,sig.m,sig25,d.shape))
    pdet.allcues <- h.fun.v(xx = rep(yys, each = length(allcues[1:3])), cc = rep(allcues[1:3], length(yys)),
                            sig.v = sig.v, sig.m = sig.m, sig25 = sig25, d.shape = d.shape)
    pdet.allcues <- matrix(pdet.allcues, nr = length(allcues[1:3]), nc = length(yys))
    ## Calculate marginal probability of finding carcass for each trip type (t = 0, 1 , etc, or by time later on)
    pmarg.bytrip <- cue.trips %*% pdet.allcues
    ddrr <- ddrr.output$ddrr
    tripmat <- ddrr.output$tripmat
    ## initialize gdrs (row is for each dd-rr combo, col is for subdivs)
    gdrs <- matrix(0, nr = nrow(ddrr), ncol = length(yys))
    ## initialize left.trunc matrix which will let us subtract off
    ## probability accumulated from effort <d.min at the end.
    left.trunc <- matrix(0, nr = nrow(ddrr), ncol = length(yys))
    for(ii in 1:ncol(tripmat))
    {
        if(verbose) print(ii)
        ## if there are still trips for that dd-rr combo
        is.trip <- ddrr$numtrips >= ii
        ## which trips are they (classified by t_i) but we add + 1 to
        ## allow us to index the t=0 trips as the first row of the
        ## p.margybytrip matrix
        which.trips <- tripmat[is.trip,ii] + 1
        ## the probability of detection is the accumulated probability
        ## of detection + (1-the accumulated probability of
        ## detection)*(the marginal probability of detection on the
        ## next trip classified by t)
        gdrs[is.trip,] <- gdrs[is.trip,] + (1-gdrs[is.trip,]) * pmarg.bytrip[which.trips,]
        ## Need to left truncate effort outside time window. To do
        ## this have a separate matrix that will save a snapshot of
        ## gdrs at the last trip before d.min (the cumulative
        ## probability of having detected a carcass between d.min)
        ## which will then be subtracted from gdrs at the end.
        befores <- as.numeric(ddrr$dd[is.trip]) + tripmat[is.trip,ii] < out$d.min #trips being considered are before
        if(ii<ncol(tripmat)) #if not at the max numtrips yet see if the next trip is after d.min
        {
            next.after <- as.numeric(ddrr$dd[is.trip]) + tripmat[is.trip,ii + 1] >= out$d.min #& it is the last trip before d.min
        }else{ #if at max numtrips, then nextafter is TRUE because we want to subtract all previous prob
            next.after <- TRUE # this way, if all trips were before d.min, left.trunc becomes the full marg prob
        }
        next.after[is.na(next.after)] <- TRUE
        if(sum(is.trip)==1)             #separate if statement if it is a vector not a matrix
        {
            left.trunc[is.trip,][befores & next.after] <- gdrs[is.trip,][befores & next.after]
        }else{
            left.trunc[is.trip,][befores & next.after,] <- gdrs[is.trip,][befores & next.after,]
        }
    }
    ## subtract off probability accumulated before d.min
    gdrs <- gdrs - left.trunc
    ## multiply by [y] and sum rows so we now have sum across dy for each ddrr combo
    ## first calculate y.pdf using gamma or unifrm
    if(gamma.y)
    {
        y.pdf <- d(temp.gamm)(yy.seq)
        ## y.pdf <- p(temp.gamm)(yy.seq + int/2) - p(temp.gamm)(yy.seq - int/2)
    }else{
        ## y.pdf <- punif(yy.seq + int/2, 0, max.dist) - punif(yy.seq - int/2, 0, max.dist)
        y.pdf <- dunif(yy.seq, 0, max.dist)
    }
    gdrs <- gdrs %*% y.pdf
    ## for each dd-rr combo, multiply by [dd][rr]
    dd.p <- ddrr$dd.p
    rr.p <- ddrr$rr.p
    ddrr.p <- dd.p*rr.p
    if(only.obs)                        # normalize by rr-dd in window
    {
        gdrs <- gdrs[ddrr$numtrips>0,]
        ddrr.p <- ddrr.p[ddrr$numtrips>0]
        ddrr.p <- ddrr.p / sum(ddrr.p)
    }
                                        # i can only do this because gdrs dd-rr organization matches that
                                        # of ddrr. and where gdrs = 0, it doesn't matter that we're
                                        # normalizing. so sum(ddrr.p) > 0 when there are trips included
                                        # that didn't have any effort in ddrr, but many are already
                                        # deleted since empty roads don't show up from non-factor usage of
                                        # xtabs above
    gdrs <- gdrs * ddrr.p
    quadsum <- sum(gdrs)*int
}

############################################################
## 9. Nll function
nll.fun <- function(pars = c(           # optimizing over pars
                    lsig.v = log(.8),
                    lsig.m = log(.3),
                    lsig25 = log(.2),
                    ld.shape = log(2)),
                    st.like = TRUE,    # use spatiotemporally explicit form of likelihood
                    gamma.y = F,
                    max.dist,
                    only.obs = F,
                    temp.gamm, cue.trips,
                    out,                # carc, seff, seffs, d.min, d.max, d.pdf, roads.df
                    ddrr.output,        # output from ddrr.fun()
                    cdr.output,         # output from carc.ddrr()
                    foreff.output,
                    show.pars = F,
                    lower.shape.bound = 2,
                    pen.scale = 8*10^4,    # scale for penalty of boundary
                    verbose = FALSE,
                    timer = TRUE,
                    browse = FALSE)
{
    if(timer) start.time <- proc.time()[3]
    if(browse) browser()
    if(show.pars) print(exp(pars))
                                        # initialize negative log likelihood
    nll <- 0
    ## add smooth penalty to ridiculous parameters
    if(pars["lsig.v"] > log(1.5))
    {
        nll <- nll + pen.scale * (pars["lsig.v"] - log(1.5))^2
        pars["lsig.v"] <- log(1.5)
    }
    if(pars["lsig.m"] > log(1.5))
    {
        nll <- nll + pen.scale * (pars["lsig.m"] - log(1.5))^2
        pars["lsig.m"] <- log(1.5)
    }
    if(pars["lsig25"] > log(1.5))
    {
        nll <- nll + pen.scale * (pars["lsig25"] - log(1.5))^2
        pars["lsig25"] <- log(1.5)
    }
    if(pars["ld.shape"] > log(5))
    {
        nll <- nll + pen.scale * (pars["ld.shape"] - log(5))^2
        pars["ld.shape"] <- log(5)
    }
    if(pars["ld.shape"] < log(lower.shape.bound))
    {
        nll <- nll + pen.scale * (pars["ld.shape"] - log(lower.shape.bound))^2
        pars["ld.shape"] <- log(lower.shape.bound)
    }
    indiv.carc.L <- g.fun(out = out, cdr.output = cdr.output, foreff.output = foreff.output,
                          sig.v = exp(pars["lsig.v"]), sig.m = exp(pars["lsig.m"]), sig25 = exp(pars["lsig25"]),
                          d.shape = exp(pars["ld.shape"]),
                          gamma.y = gamma.y, temp.gamm = temp.gamm, cue.trips = cue.trips,
                          max.dist = max.dist, st.like = st.like,
                          browse = F)
    if(st.like)
    {
        ## The denominator of the likelihood (proportion of carcasses detected)
        quadsum <- g.dk(ddrr.output, out, yys = yy.seq,
                        sig.v = exp(pars["lsig.v"]), sig.m = exp(pars["lsig.m"]), sig25 = exp(pars["lsig25"]),
                        d.shape = exp(pars["ld.shape"]),
                        gamma.y = gamma.y, temp.gamm = temp.gamm, max.dist = max.dist, only.obs = only.obs,
                        cue.trips = cue.trips,
                        browse = F)
        ## The numerator of the likelihood: comes from the carcass data
        nll <- nll + nrow(out$carc) * log(quadsum) - sum(log(indiv.carc.L))
    }else{
        denoms <-  pr.fun(out = out, cdr.output = cdr.output,
                          ddrr.output = ddrr.output,
                          sig.v = exp(pars["lsig.v"]), sig.m = exp(pars["lsig.m"]), sig25 = exp(pars["lsig25"]),
                          d.shape = exp(pars["ld.shape"]), yys = yy.seq,
                          gamma.y = gamma.y, temp.gamm = temp.gamm, max.dist = max.dist,
                          cue.trips = cue.trips,
                          browse = F)
        nll <- nll + sum(log(denoms)) - sum(log(indiv.carc.L))
    }
    if(timer) end.time <- proc.time()[3]
    if(timer) print(paste("nll calc took", end.time - start.time, "s"))
    return(nll)
}

## ## Vectorize nll over parameters if interested in looking at contours
## nll.w <- function(x,y, pars,...)
## {
##     ## pars["lsig.v"] <- log(x)
##     pars["ld.shape"] <- log(y)
##     nll.fun(pars, out, ddrr.output = dr.out, cdr.output = cdr.out, show.pars = F, verbose = F, timer = F, browse = F)
## }
## nll.v <- Vectorize(nll.w, vectorize.args=c("x","y"))

## xs <- seq(.2,1,l=10)
## ys <- seq(.5,9,l=10)
## contour(outer(xs,ys,nll.v, pars = init.true))

## Simulate carcass data given real effort data in a given time window, rrs, etc
data.sim <- function(rrs = 192,
                     only.in.eff = TRUE,
                     ddrr.str = NULL,
                     d.min = min(seff$nday) + 5,
                     d.max = min(seff$nday) + 50,
                     set.dist = F,     # fix dist for troubleshooting
                     to.set = .1,
                     set.dd = F,
                     rand.dd.p = F,     # add temporal variability in carcass distribution
                     rand.rr.p = F,     # add spatial variability in carcass distribution (note no s-t variability at this point)
                     wrong.rr.p = F,    # use uniform rr.p to fit heterogenous rr.p
                     dd.set = NULL,
                     get.eff=F,
                     sim.vals = c(sig.v = NULL, sig.m = NULL, sig25 = NULL, d.shape = NULL),
                     max.dist, btrack = 5,
                     gamma.y = F, temp.gamm, cue.trips,
                     N = 300, browse = F)
{
    if(browse) browser()
    ## Load road length and calculate [r] proportional to road length (for now)
    roads.df <- roads.df.all[roads.df.all$rdID %in% rrs,]
    roads.df$p <- roads.df$LENGTH / sum(roads.df$LENGTH)
    seff.int <- seff
    seff.int <- seff.int[seff.int$Road %in% rrs,]
    seff.int <- seff.int[seff.int$nday <= d.max & seff.int$nday >= (d.min- btrack),]
    if(rand.rr.p)
    {
        ## multiply by gamm distr scaling factor (retain original 1/length to preserve that longer roads still cover more area)
        roads.df$p <- roads.df$p * rgamma(nrow(roads.df), shape = 1, scale = .5) # aggregated carcass distr in time
        roads.df$p <- roads.df$p / sum(roads.df$p)
        ## roads.df$p[1] <- .99
        ## roads.df$p[-1] <- rep(.01, nrow(roads.df)-1)
        ## roads.df$p <- roads.df$p / sum(roads.df$p)
    }
    dseq <- (d.min-btrack):d.max
    d.pdf <- data.frame(date = dseq, p = 1/length(dseq))
    d.pdf$p[dseq<d.min] <- 0
#    d.pdf$p[(d.max - dseq)<6 ] <- 0
    d.pdf$p <- d.pdf$p / sum(d.pdf$p)
    if(rand.dd.p)
    {
        d.pdf$p <- rgamma(length(dseq), shape = 1, scale = .5) # aggregated carcass distr in time
        d.pdf$p[dseq<d.min] <- 0
        d.pdf$p <- d.pdf$p / sum(d.pdf$p)
    }
    ## d.pdf <- rbind(data.frame(date = (min(dseq)-btrack):(min(dseq)-1), p = 0),d.pdf)
    ## Set up matrix to use for runif(1) easily so that if x < first val
    ## => avian, etc
    ## temp.frame <- day.frame[1:(btrack+1),allcues[1:3]]
    temp.frame <- cue.trips
    crand.frame <- cue.trips
    crand.frame[,2] <- rowSums(temp.frame[,1:2])
    crand.frame[,3] <- rowSums(temp.frame)
    crand.frame <- data.frame(day = 0:btrack, crand.frame)
    ## Initialize the covariates
    if(gamma.y)
    {             # simulate from truncated gamma using library(distr)
        yy <- r(temp.gamm)(N)
    }else{                              # uniform dist
        yy <- runif(N, 0, max.dist)
    }
    if(!only.in.eff)
    {
        dd <- rep(dseq, as.vector(rmultinom(1, N, prob = d.pdf$p))) # floor(runif(N, min(dseq), max(dseq)+1))
        rr <- rep(roads.df$rdID, as.vector(rmultinom(1,N, prob = roads.df$p)))
    }else{
        ## if only generating carcasses in the period with effort
        dds.rep <- d.pdf[rep(1:nrow(d.pdf), each = nrow(roads.df)),]
        rrs.rep <- roads.df[rep(1:nrow(roads.df), nrow(d.pdf)),]
        rrs.rep <- rrs.rep[,5:6]
        ddrr.dist <- cbind(dds.rep, rrs.rep)
        ddrr.dist$drp <- ddrr.dist[,2]*ddrr.dist[,3]
        ddrr.all <- paste(ddrr.dist[,"rdID"], ddrr.dist[,"date"])
        in.eff <- ddrr.all %in% ddrr.str
        ddrr.dist <- ddrr.dist[in.eff,]
        ddrr.dist$drp <- ddrr.dist$drp / sum(ddrr.dist$drp)
        indices <- rep(1:nrow(ddrr.dist),  as.vector(rmultinom(1, N, prob = ddrr.dist$drp)))
        dd <- ddrr.dist$date[indices]
        rr <- ddrr.dist$rdID[indices]
    }
                                        #if(set.dist) yy <- rep(to.set, N)
                                        #if(set.dd) dd <- rep(dd.set, N)
    ## if analyzing with the wrong rr.p, reset to based on length
    if(wrong.rr.p) roads.df$p <- roads.df$LENGTH / sum(roads.df$LENGTH)
    det <- rep(FALSE,N)                         #if detected
    cc <- rep(NA,N)
    ll <- rep(NA,N)
    tt <- rep(NA,N)
    driver <- factor(rep(NA,N), levels = levels(seff.int$Driver))
    tin <- as.POSIXct(rep(NA,N))
    tout <- as.POSIXct(rep(NA,N))
    tms <- rep(NA,N)
    obs <- rep(NA,N)
    dateTime <- as.POSIXct(rep(NA,N) )
    tod <- rep(NA,N)
    ## loop through carcasses
    for(ii in 1:N)
    {                                       # find seff.int in window of dd
        day.diff <-  seff.int$nday - dd[ii]
        in.window <- day.diff >= 0 & day.diff <= btrack
        in.window <- in.window & seff.int$Road == rr[ii]
        if(sum(in.window)>0)
        {                                   # if there is surv eff in window, order it
            rest.seff.int <- seff.int[in.window,]
            rest.seff.int$dsd <- rest.seff.int$nday - dd[ii]
            rest.seff.int <- rest.seff.int[order(rest.seff.int$tout),]
            for(jj in 1:nrow(rest.seff.int))
            { # if there is seff and carcass hasn't been detected yet
                if(!det[ii])
                { # determine which cue is present with random uniform variable
                    rand <- runif(1)
                    if(rand < crand.frame[crand.frame$day == rest.seff.int$dsd[jj],"avian"])
                    {temp.cue <- "avian"
                 }else{
                     if(rand < crand.frame[crand.frame$day == rest.seff.int$dsd[jj],"mamm"])
                     {temp.cue <- "mamm"
                  }else{
                      if(rand < crand.frame[crand.frame$day == rest.seff.int$dsd[jj],"cstate25"])
                      {temp.cue <- "cstate25"
                   }else{
                       temp.cue <- NA
                   }}}
                    ## if there was a cue, attempt detection
                    if(!is.na(temp.cue))
                    { # calls h.fun to get det probs
                        temp.p <- h.fun(yy[ii], temp.cue, sig.v = sim.vals["sig.v"], sig.m = sim.vals["sig.m"],
                                        sig25 = sim.vals["sig25"], d.shape= sim.vals["d.shape"])
                        det[ii] <- rbinom(1,1, temp.p)==1
                        if(det[ii]) # if detected assign covariates relating to detection
                        {
                            cc[ii] <- temp.cue
                            ll[ii] <- rest.seff.int$nday[jj]
                            tt[ii] <- rest.seff.int$dsd[jj]
                            driver[ii] <- rest.seff.int$Driver[jj]
                            tin[ii] <- rest.seff.int$tin[jj]
                            tout[ii] <- rest.seff.int$tout[jj]
                            tms[ii] <- rest.seff.int$tms[jj]
                            obs[ii] <- rest.seff.int$obs[jj]
                            dateTime[ii] <- format(as.POSIXct(rest.seff.int$tin[jj]), "%Y-%m-%d")
                            ## pick random time during trip for carcass to have been discovered
                            tin.t <- as.numeric(format(tin[ii], "%H")) + as.numeric(format(tin[ii], "%M"))/60
                            tout.t <- as.numeric(format(tout[ii], "%H")) + as.numeric(format(tout[ii], "%M"))/60
                            tod[ii] <- runif(1,tout.t, tin.t)
                                        #                            if(is.na(tod[ii])) print("tod na"); browser()
                          }#end det if
                      }    #end cue if
                  }        #end cue assignment
              }            #end seff.int loop
          }                #end if seff.int exists statement
      }                    #end loop through carcasses
    carc.all <- data.frame(det = det, nday = ll, cue = cc, roaddist = yy,
                           roadwhich = rr, dd = dd,ll = ll, tt = tt, driver = driver, tms = tms, obs = obs,
                           tout = tout,tin = tin, tod = tod, dateTime = dateTime)
    carc <- carc.all[det & ll >= d.min,] # detected after d.min is the subst within the seff window
    dat <- list(carc, seff.int)
    seffs <- NA
    if(get.eff)    seffs <- eff.list(dat, btrack = btrack, browse = F)
    ## Figure out non empty effort matrix
    return(list(carc = carc, seff =  seff.int, seffs = seffs, d.pdf = d.pdf, roads.df = roads.df, carc.all = carc.all,
                d.min = d.min, d.max = d.max))
  }                                       # end

## Optimize detectaility (and other) parameters on simulated data sets
optnll <- function(true.vals = NULL, #seed = 1, # for multicore
                   real.dat = F,
                   only.obs = F,        # only estimate N for carcasses in surv eff
                   only.in.eff = T,     # only generate simulated carcasses in d-r with effort
                   ddrr.str = NULL,     # day-roads with effort in next btrack days, fed in if only.in.eff=T
                   out.real =  NULL,     # feed in real dat if using it
                   which.pars = 1:4,    # indices of parameters to fit
                   st.like = T,         # maximize st likelihood or not if F
                   st.est = T,          # use spatiotemporally explicity estimator or not if F
                   spec.fix.pars = NULL, # values of parameters to fix
                   rand.dd.p = F,
                   rand.rr.p = F, wrong.rr.p = F,
                   hess.conf = T, boot = F, nboot = 1, # bootstrap confidence intervals
                   boot.y = F, # resample from fitted [y] during bootstrap
                   gamm.boot.pars = paste(file.folder, "bframe.wet.Rdata", sep = ""), # file with [y] gamm pars
                   cue.trips.file = paste(file.folder, "cue.trips.Rdata",sep = ""), # [c|t]
                   boot.c = F, # resample from fitted [c|t] during bootstrap                  
                   cam.boot.frame = paste(file.folder, "day.frame.b.Rdata",sep = ""),
                   lower = .025, upper = .975,
                   make.pdf = F,
                   file.name = paste(file.folder, "output", sep = ""),
                   sd.pert = .1,        # std dev of perturbation on log scale from true par values for initial values during fit
                   d.min, d.max, rrs, N = 600, # time window, roads, and # carcasses
                   nn = 4,                     # number of simulations
                   max.dist,
                   quantsims = 300,            # number of quadratures to calculate to get confidence intervals
                   vectorized.quad = T,        # calculate quadratures in vectorized format (can't do for big ones)
                   load.dr = NA,               # give path to dr.out if already done
                   do.sann = T,                            #optimize with sann first
                   hess = T,                               #calculate hessian (necessary for conf int)
                   gamma.y = F,                            #use gamma for [y]
                   g.shape, g.scale,                       #gamma pars
                   nm.its = 350,                           #nelder mead iterations
                   tr = 6,                                 #optim verbose par
                   breaks = seq(0,1.2, by = .05),          #hist breaks for plots
                   like.pts = T,                           #plot pdf of y,c on histogram
                   sann.its = 60,                          #sann iterations
                   reltol = 10^-6,                         #tolerance for optimization
                   show.pars = F, verbose = F, timer = F, browse = F) # verbose parameters
{
  graphics.off()
  if(browse) browser()
  load(cue.trips.file)                 # loads cue.trips [c|t]
  cue.trips <- as.matrix(cue.trips)
  temp.gamm <- Truncate(Gammad(shape = g.shape, scale = g.scale), lower = 0, upper = max.dist) # create truncated gamma distribution for [y]
  ## create a new file name for output (to make sure we never copy over old sims)
  stepper <- 0
  file.original <- file.name
  while(file.exists(paste(file.name, ".Rdata", sep = "")))
    {
      stepper <- stepper + 1
      file.name <- paste(file.original, "-", stepper, sep = "")
    }
  pdf.name <- paste(file.name, ".pdf", sep = "")
  save.file <- paste(file.name, ".Rdata", sep = "")
  ## create nll wrapper to feed fixed & fitted pars appropriately
  nll.chosen <- function(pars, fix.pars,...)
    {
      pars <- c(pars, fix.pars)
      nll.fun(pars, ...)
    }
  ## initialize plot
  if(make.pdf)     pdf(pdf.name)
  par(mfrow=c(2,1), mar = c(3,3,3,1.5), oma = rep(0,4))
  ymax <- 1.5*N
  if(real.dat) ymax <- 8*nrow(out.real$carc)
  plot(0,0, xlim = c(1, nn), ylim = c(0,2.5*N), type = "n", bty = "n", xlab = "simulation", xaxt="n",
       ylab = "# carcasses")
  if(!real.dat) abline(h=N)
  ## start sims
  sims <- rbind(data.frame(N = rep(N, nn), n.in.eff = NA, n.obs = NA,
                           ml.est = NA, ml.l = NA, ml.u = NA,
                           ml.est.b = NA, ml.l.b = NA, ml.u.b = NA,
                           min = NA, conv = NA, its = NA,
                           nll = NA, true.est = NA))
  pars.frame <- data.frame(rep(NA,nn),rep(NA,nn),rep(NA,nn),rep(NA,nn))
  inits.frame <- data.frame(rep(NA,nn),rep(NA,nn),rep(NA,nn),rep(NA,nn))
  init.true <- log(true.vals)
  names(init.true) <- c("lsig.v","lsig.m","lsig25","ld.shape")
  names(pars.frame) <- names(init.true)
  names(inits.frame) <- names(init.true)
  outs <- list()
  for(jj in 1:nn)
    {
      temp.start <- proc.time()[3]
      print(paste("working on sim", jj,"of",nn))
      if(!real.dat)
        {
          out <- data.sim(sim.vals = true.vals, ddrr.str = ddrr.str, only.in.eff = only.in.eff,
                          d.max = d.max, d.min = d.min, get.eff = F, rrs = rrs, N = N,
                          rand.dd.p = rand.dd.p, rand.rr.p = rand.rr.p, wrong.rr.p = wrong.rr.p, max.dist = max.dist,
                          cue.trips = cue.trips,
                          gamma.y = gamma.y, temp.gamm = temp.gamm)
          outs[[jj]] <- out
        }else{
          out.real$carc <- out.real$carc[out.real$carc$roaddist < max.dist,]
          outs[[jj]] <- out.real
          out <- out.real
        }
      if(jj==1)
        {#only need to do this once for each rr-dd combo (even if doing multiple sims)
                                        #            dr.out <- ddrr.fun(out, btrack = 5)
          dr.out <- ddrr.select(out, btrack = 5, truncate.end = T)
        }
      ddrr.str <- paste(dr.out[["ddrr"]]$rr, dr.out[["ddrr"]]$dd) #for matching to n.in.eff below
      ddrr.str <- ddrr.str[dr.out[["ddrr"]]$numtrips>0]
      p.marg <- g.dk(ddrr.output = dr.out, out = out, yys = yy.seq, sig.v = true.vals["sig.v"],
                     sig.m = true.vals["sig.m"], sig25 = true.vals["sig25"], d.shape = true.vals["d.shape"],
                     gamma.y = gamma.y, temp.gamm = temp.gamm, max.dist = max.dist, cue.trips = cue.trips,
                     verbose = F, browse = F)
                                        #        if(!real.dat) abline(h=N*p.marg, col = "blue")
      cdr.out <- carcddrr(out, btrack = 5, verbose = F, browse = F) # d.min/max??
      foreff.out <- carcdrforward(out = out, ddrr.output = dr.out, btrack = btrack, verbose = F, browse = F)
      sims[jj,"n.obs"] <- nrow(out$carc)
      ##how many carcasses occurred in rr-dd areas (in the simulation) that were in {d,r} st r was driven in [d,d+5]?
      sims[jj,"n.in.eff"] <- sum(paste(out$carc.all$roadwhich,out$carc.all$dd) %in% ddrr.str)
      ## random perturbed initial values
      inits <- init.true[which.pars]
      inits <- inits + rnorm(length(which.pars), 0, sd.pert)
      inits.frame[jj,which.pars] <- inits
      fix.pars <- init.true[!1:4 %in% which.pars]
      for(pp in 1:length(spec.fix.pars)) #if fixing any of the pars at the wrong value
        {
          fix.pars[names(spec.fix.pars)[pp]] <- spec.fix.pars[names(spec.fix.pars)[pp]]
          inits.frame[jj,names(spec.fix.pars)] <- spec.fix.pars[names(spec.fix.pars)[pp]]
        }
######################################################################
######################################################################
      ## Begin optimzation
######################################################################
      if(do.sann)
        {
          ## Optimize with SANN first
          opt.sann <- optim(par = inits,
                            nll.chosen,
                            fix.pars = fix.pars,
                            timer=F, verbose=F,
                            out = out, ddrr.output = dr.out, cdr.output = cdr.out, foreff.output = foreff.out,
                            gamma.y = gamma.y, temp.gamm = temp.gamm, max.dist = max.dist, cue.trips = cue.trips,
                            st.like = st.like,
                            method = "SANN", show.pars = show.pars,
                            hessian = hess,
                            control = list(trace = tr, maxit = sann.its, reltol = reltol))
          next.pars <- opt.sann$par
        }else{
          next.pars <- inits
        }
      ## Optimize with Nelder-Mead
      opt <- optim(par = next.pars,
                   nll.chosen,
                   fix.pars = fix.pars,
                   timer=F, verbose=F,
                   out = out, ddrr.output = dr.out, cdr.output = cdr.out, foreff.output = foreff.out,
                   gamma.y = gamma.y, temp.gamm = temp.gamm, max.dist = max.dist, cue.trips = cue.trips,
                   st.like = st.like,
                   method = "Nelder-Mead", show.pars = show.pars,
                   hessian = hess,
                   control = list(trace = tr, maxit = nm.its, reltol = reltol))
######################################################################
      realval <- nll.fun(init.true, out = out, ddrr.output = dr.out, cdr.output = cdr.out, foreff.output = foreff.out,
                         gamma.y = gamma.y, temp.gamm = temp.gamm, st.like = st.like, max.dist = max.dist,
                         cue.trips = cue.trips,
                         show.pars = F, verbose = F, timer = F, browse = F)
      print(paste("nll at true pars is", realval))
      sims[jj,"nll"] <- opt$value
      sims[jj,"its"] <- opt$count[1]
      sims[jj,"conv"] <- opt$conv
      pars.frame[jj,names(opt$par)] <-  opt$par
      pars.frame[jj, names(fix.pars)] <- fix.pars
      print(exp(pars.frame[jj,]))
      temp.covmat <- as.matrix(nearPD(solve(opt$hessian))$mat)
######################################################################
######################################################################
      ## Hessian confidence intervals (calls confinter())
      if(hess.conf)
        {
######################################################################
          sims[jj,c("ml.est","ml.l","ml.u")] <- confinter(outs[[jj]], cdr.output = cdr.out, ddrr.output = dr.out,
                                                          pars = opt$par, fix.pars = fix.pars,
                                                          only.obs = only.obs, st.est = st.est,
                                                          gamma.y = gamma.y, temp.gamm = temp.gamm, vectorized.quad = vectorized.quad,
                                                          max.dist = max.dist, cue.trips = cue.trips,
                                                          covmat = temp.covmat, numsim = quantsims, browse=F)
          points(jj-.05, sims$ml.est[jj], col="gray", pch = 19, cex=.7)
          arrows(jj-.05, sims$ml.l[jj], jj-.05, sims$ml.u[jj], angle = 90, length = .05, code = 3, col = "gray")
        }
######################################################################
######################################################################
                                        #bootstrap confidence intervals
######################################################################
      if(boot)
        {
          mcpr.boot <- rep(NA,nboot)
          if(real.dat)
            {
              par.boot <- data.frame(x1 = rep(NA,nboot), x2 = NA, x3 = NA, x4 = NA)
              names(par.boot) <- names(init.true)
            }else{
              par.boot <- NULL
            }
          for(bb in 1:nboot)
            {
              if(bb%%10 == 1) print(paste("Bootstrap", bb, "of", nboot))
                                        #            browser()
              out.temp <- out
              boot.sel <- sample(1:nrow(out$carc), nrow(out$carc), replace = TRUE)
              out.temp$carc <- out.temp$carc[boot.sel,]
              cdr.temp <- cdr.out[boot.sel,]
              ## selecting from foreff is slightly tricky given that each carc has (btrack + 1) rows
              inflater <- (boot.sel-1)*(btrack+1) + 1
              boot.sel.inflater <- as.vector(sapply(inflater, function(x) {x:(x+btrack)}))
              foreff.temp <- foreff.out
              foreff.temp[[1]] <- foreff.temp[[1]][boot.sel.inflater,]
              foreff.temp[[2]] <- foreff.temp[[2]][boot.sel.inflater,]
              start.pars <- pars.frame[jj,which.pars]      # start at parameters at mle to save computing time
              ## if resampling over parameters from the gamma [y] during bootstrap
              if(boot.y)
                {
                  load(gamm.boot.pars) # should load bframe.wet
                  rand.choice <- sample(1:nrow(bframe.wet), 1)
                  g.shape.boot <- bframe.wet[rand.choice, "shape"]
                  g.scale.boot <- bframe.wet[rand.choice, "scale"]
                  temp.gamm.boot <- Truncate(Gammad(shape = g.shape.boot, scale = g.scale.boot), lower = 0, upper = max.dist)
                }else{                  #otherwise use the original gamma distribution
                  temp.gamm.boot <- temp.gamm
                }
              if(boot.c)
                {
                  load(cam.boot.frame) # loads day.frame.b, an array of cue.trips ([c|t]) fitted while bootstrapping camera trap data
                  rand.choice <- sample(1:dim(day.frame.b)[3], 1)
                  cue.frame.temp <- day.frame.b[,,rand.choice] #select the matrix within the array of [c|t] as loaded
                }else{
                  cue.frame.temp <- cue.trips
                }
              if(do.sann)
                {
                  ## Optimize with SANN first
                  opt.sann <- optim(par = start.pars,
                                    nll.chosen,
                                    fix.pars = fix.pars,
                                    timer=F, verbose=F,
                                    out = out.temp, ddrr.output = dr.out, cdr.output = cdr.temp, foreff.output = foreff.temp,
                                    gamma.y = gamma.y, temp.gamm = temp.gamm.boot, max.dist = max.dist,
                                    cue.trips = cue.frame.temp,
                                    st.like = st.like,
                                    method = "SANN", show.pars = show.pars,
                                    hessian = FALSE, # don't need Hessian for bootstrap
                                    control = list(trace = tr, maxit = sann.its, reltol = reltol))
                  next.pars <- opt.sann$par
                }else{
                  next.pars <- inits
                }
              ## Optimize with Nelder-Mead
              opt <- optim(par = next.pars,
                           nll.chosen,
                           fix.pars = fix.pars,
                           timer=F, verbose=F,
                           out = out.temp, ddrr.output = dr.out, cdr.output = cdr.temp, foreff.output = foreff.temp,
                           gamma.y = gamma.y, temp.gamm = temp.gamm.boot, max.dist = max.dist,
                           cue.trips = cue.frame.temp,                             
                           st.like = st.like,
                           method = "Nelder-Mead", show.pars = show.pars,
                           hessian = FALSE,   # don't need for bootstrap
                           control = list(trace = tr, maxit = nm.its, reltol = reltol))
              if(real.dat)
                {
                  par.boot[bb,which.pars] <- opt$par
                  par.boot[bb, !1:4 %in% which.pars] <- init.true[!1:4 %in% which.pars]
                }
              pars.temp <- opt$par
              pars.temp[!1:4 %in% which.pars] <- init.true[!1:4 %in% which.pars]
              names(pars.temp)[!1:4 %in% which.pars] <- names(init.true)[!1:4 %in% which.pars]
              mcpr.boot[bb] <- sum(1/ pr.fun.yspec(out.temp, cdr.temp, dr.out,
                                                   sig.v = exp(pars.temp["lsig.v"]), sig.m = exp(pars.temp["lsig.m"]),
                                                   sig25 = exp(pars.temp["lsig25"]), d.shape= exp(pars.temp["ld.shape"]),
                                                   yys = yy.seq,       # quadrature subdivisions
                                                   gamma.y = gamma.y,
                                                   cue.trips = cue.frame.temp,
                                                   temp.gamm = temp.gamm, browse = F))
            }
          if(real.dat) par.boot <- data.frame(par.boot, N = mcpr.boot)
          sims[jj,c("ml.l.b","ml.u.b")] <-  quantile(mcpr.boot, c(lower, upper))
          sims[jj,"ml.est.b"] <- sum(1/ pr.fun.yspec(out, cdr.out, dr.out,
                                                     sig.v = exp(pars.frame[jj,"lsig.v"]), sig.m = exp(pars.frame[jj,"lsig.m"]),
                                                     sig25 = exp(pars.frame[jj,"lsig25"]), d.shape= exp(pars.frame[jj,"ld.shape"]),
                                                     yys = yy.seq,       # quadrature subdivisions
                                                     gamma.y = gamma.y,
                                                     cue.trips = cue.trips,
                                                     temp.gamm = temp.gamm, browse = F))
          points(jj+.05, sims$ml.est.b[jj], pch = 19, cex=.7)
          arrows(jj+.05, sims$ml.l.b[jj], jj+.05, sims$ml.u.b[jj], angle = 90, length = .05, code = 3)
######################################################################
        }                                   # end bootstrap
######################################################################
######################################################################
      if(!st.est) sims[jj, "true.est"] <- sum(1/pr.fun.yspec(out, cdr.output = cdr.out, ddrr.output = dr.out,
                                                             sig.v = exp(init.true["lsig.v"]), sig.m = exp(init.true["lsig.m"]),
                                                             sig25 = exp(init.true["lsig25"]), d.shape= exp(init.true["ld.shape"]),
                                                             yys = yy.seq,       # quadrature subdivisions
                                                             gamma.y = gamma.y,
                                                             cue.trips = cue.trips,
                                                             temp.gamm = temp.gamm, browse = F))
      sims[jj, "min"] <- (proc.time()[3] - temp.start)/60
      points(jj, sims$n.obs[jj], col="blue", pch = 19)
      if(!st.est & !real.dat) points(jj, sims$true.est[jj], col="purple", pch = 21)
      temp.hess <<- opt$hess
      if(!real.dat)  points(jj, sims$n.in.eff[jj], col="brown", pch = 21)
      if(!real.dat)
        {
          legend("topleft", c("obs carc","carc in eff", "estimator at true vals", "estimator at mle"), pch = c(19,21,21,19),
                 col  = c("blue","brown","purple","black"), cex = .7, bty = "n")
        }else{
          legend("topleft", c("obs carc","estimator at mle"), pch = c(19,19),
                 col  = c("blue","black"), cex = .7, bty = "n")
        }
      legend("topright", c("hessian CI", "bootstrap CI"), col = c("gray","black"), lty = 1, cex = .7, bty = "n")
    }
  if(real.dat) print(pars.frame)
  test.n <- N
  if(only.obs) test.n <- sims$n.in.eff
  p.val <- sum(sims$ml.l > test.n | sims$ml.u < test.n) / nrow(sims)
  mtext(paste("P =", signif(p.val,2)), 3, 1)
  plot(0,0, xlim = c(1,4), ylim = c(0,1.5), type = "n", xaxt="n", bty = "n", xlab = "parameter")
  axis(1, 1:4, names(true.vals))
  points(1:3, exp(init.true)[1:3], pch = 21, col = "black")
  points(4, exp(init.true)[4]/4, pch = 21, col = "black")
  for(ii in 1:nn)
    {
      xs <- 1:3 + (ii-1)/(nn*2)
      points(xs,exp(pars.frame[ii,1:3]) , pch = 19, cex = .5, col="red")
      xs <- 4 + (ii-1)/(nn*2)
      points(xs,exp(pars.frame[ii,4])/4 , pch = 19, cex = .5, col = "blue")
    }
  axis(4, seq(0,1.5, by = .5), seq(0,1.5, by = .5)*4, col = "blue")
######################################################################
  pars.frame <- exp(pars.frame)
  inits.frame[names(fix.pars)] <- fix.pars
  inits.frame <- exp(inits.frame)
  colnames(pars.frame) <- names(true.vals)
  if(make.pdf)
    {
      print("plotting conditional detection pdfs vs data")
      for(jj in 1:nrow(pars.frame))
        {
          out <- outs[[jj]]
          cdr.out <- carcddrr(out, btrack = 5, verbose = F, browse = F) # d.min/max??
          foreff <- carcdrforward(out = out, ddrr.output = dr.out, btrack = btrack, verbose = F, browse = F)
          if(st.like)
            {
              denom <- g.dk(ddrr.output = dr.out, out = out, yys = yy.seq,
                            sig.v = true.vals["sig.v"], sig.m = true.vals["sig.m"],
                            sig25 = true.vals["sig25"], d.shape = true.vals["d.shape"],
                            gamma.y = gamma.y, temp.gamm = temp.gamm, max.dist = max.dist,
                            cue.trips = cue.trips,
                            only.obs = only.obs,
                            verbose = F, browse = F)
              denom.opt <- g.dk(ddrr.output = dr.out, out = out, yys = yy.seq,
                                sig.v = pars.frame[jj,"sig.v"], sig.m = pars.frame[jj,"sig.m"],
                                sig25 = pars.frame[jj,"sig25"], d.shape = pars.frame[jj,"d.shape"],
                                gamma.y = gamma.y, temp.gamm = temp.gamm, max.dist = max.dist,
                                cue.trips = cue.trips,
                                only.obs = only.obs,
                                verbose = F, browse = F)
            }else{
              denom <- pr.fun(out, cdr.output = cdr.out, ddrr.output = dr.out,
                              sig.v = true.vals["sig.v"], sig.m = true.vals["sig.m"],
                              sig25 = true.vals["sig25"], d.shape = true.vals["d.shape"],
                              yys = yy.seq,       # quadrature subdivisions
                              gamma.y = gamma.y, max.dist = max.dist,
                              cue.trips = cue.trips,                                
                              temp.gamm = temp.gamm, browse = F)
              denom.opt <- pr.fun(out, cdr.output = cdr.out, ddrr.output = dr.out,
                                  sig.v = pars.frame[jj,"sig.v"], sig.m = pars.frame[jj,"sig.m"],
                                  sig25 = pars.frame[jj,"sig25"], d.shape = pars.frame[jj,"d.shape"],
                                  yys = yy.seq,       # quadrature subdivisions
                                  gamma.y = gamma.y, max.dist = max.dist,
                                  cue.trips = cue.trips,                                    
                                  temp.gamm = temp.gamm, browse = F)
            }
          gs <- g.fun(out, cdr.out, foreff, sig.v = true.vals["sig.v"], sig.m = true.vals["sig.m"], sig25 = true.vals["sig25"],
                      d.shape = true.vals["d.shape"], max.dist = max.dist, cue.trips = cue.trips, st.like = st.like, browse = F) / as.vector(denom)
          gs.opt <- g.fun(out, cdr.out, foreff, sig.v = pars.frame[jj,"sig.v"], sig.m = pars.frame[jj,"sig.m"],
                          sig25 = pars.frame[jj,"sig25"], d.shape = pars.frame[jj,"d.shape"],
                          max.dist = max.dist, st.like = st.like,  cue.trips = cue.trips, browse = F) / as.vector(denom.opt)
          ## and prob of observing each carc
          mcprs <- pr.fun.yspec(out, cdr.output = cdr.out, ddrr.output = dr.out,
                                sig.v = true.vals["sig.v"], sig.m = true.vals["sig.m"],
                                sig25 = true.vals["sig25"], d.shape = true.vals["d.shapee"],
                                yys = yy.seq,       # quadrature subdivisions
                                gamma.y = gamma.y,
                                cue.trips = cue.trips,
                                temp.gamm = temp.gamm, browse = F)
          mcprs.opt <- pr.fun.yspec(out, cdr.output = cdr.out, ddrr.output = dr.out,
                                    sig.v = pars.frame[jj,"sig.v"], sig.m = pars.frame[jj,"sig.m"],
                                    sig25 = pars.frame[jj,"sig25"], d.shape = pars.frame[jj,"d.shape"],
                                    yys = yy.seq,       # quadrature subdivisions
                                    gamma.y = gamma.y,
                                    cue.trips = cue.trips,
                                    temp.gamm = temp.gamm, browse = F)
          par(mfrow=c(4,2), oma = c(1.5,0,3,.5), mar = c(2.5,3,2,1.7))
          ymax <- 0
          for(cc in allcues[1:3])
            {
              hist1 <- hist(out$carc$roaddist[out$carc$cue==cc], breaks = breaks, plot = F)
              hist1$counts <- hist1$counts/nrow(out$carc)
              ymax <- max(ymax,hist1$counts)*1.05
            }
          for(cc in allcues[1:3])
            {
              xs <- seq(.01,1.2, l = 200)
              hist1 <- hist(out$carc$roaddist[out$carc$cue==cc], breaks = breaks, plot = F)
              hist1$counts <- hist1$counts/nrow(out$carc)
              plot(hist1, bty = "n", col = "yellow", ylim = c(0,ymax),
                   xlab = "", main = cc)
              if(cc == "avian") scalar <- max(hist1$counts)/max(gs.opt[names(gs)==cc])
              if(like.pts)
                {
                  if(!real.dat) points(out$carc$roaddist[out$carc$cue == cc], gs[names(gs) == cc]*scalar, pch = 21, cex = .8,
                                       col = "purple")
                  points(out$carc$roaddist[out$carc$cue == cc], gs.opt[names(gs.opt) == cc]*scalar, pch = 21, cex = .8,
                         col = "red")
                }
              ## add det fxn at true pars
              xs <- seq(0,max.dist, l=100)
              ys <- h.fun(xs, cc = cc, sig.v = true.vals["sig.v"], sig.m = true.vals["sig.m"], sig25 = true.vals["sig25"],
                          d.shape = true.vals["d.shape"])
              ys <- ys/max(ys)*ymax/1.05
              if(!real.dat) lines(xs, ys, col = "purple")
              ## add det fxn at init pars
              ys.inits <- h.fun(xs, cc = cc, sig.v = inits.frame[jj,"lsig.v"], sig.m = inits.frame[jj,"lsig.m"],
                                sig25 = inits.frame[jj,"lsig25"], d.shape = inits.frame[jj,"ld.shape"])
              ys.inits <- ys.inits/max(ys.inits)*ymax/1.05
              lines(xs, ys.inits, col = "green")
              ## add det fxn at opt pars
              ys.opt <- h.fun(xs, cc = cc, sig.v = pars.frame[jj,"sig.v"], sig.m = pars.frame[jj,"sig.m"],
                              sig25 = pars.frame[jj,"sig25"], d.shape = pars.frame[jj,"d.shape"])
              ys.opt <- ys.opt/max(ys.opt)*ymax/1.05
              lines(xs, ys.opt, col = "red")
              axis(4, at = seq(0,max(ys), l = 3), seq(0,1, l=3), col = "red", las = 2)
              if(cc == "cstate25") legend("topright", c("true det funx","fitted det funx", "initial det funx for optim"),
                   lwd = 2,bty="n", col = c("purple","red","green"))
              ## axis(4, tics*scalar, tics, las = 2)
              plot(hist1, bty = "n", col = "yellow", ylim = c(0,ymax),
                   xlab = "", main = cc)
              if(cc == "avian") scalar <- max(hist1$counts)/max(mcprs.opt[names(gs)==cc])
              if(like.pts)            # add prob of observing each carc
                {
                  if(!real.dat) points(out$carc$roaddist[out$carc$cue == cc], mcprs[names(gs) == cc]*scalar, pch = 21, cex = .8,
                                       col = "purple")
                  points(out$carc$roaddist[out$carc$cue == cc], mcprs.opt[names(gs.opt) == cc]*scalar, pch = 21, cex = .8,
                         col = "red")
                }
              abline(h = max(ys) / 100, lty = 3, col = "black")
              axis(4, at = seq(0,max(ys), l = 3), seq(0,1, l=3), col = "red", las = 2)
              ## tics <- pretty(c(0,ymax*.8/scalar),4)
              ## if(cc == "cstate25") legend("topright", c("fitted det funx", "initial det funx for optim"), lwd = 2,bty="n",
              ##    col = c("red","green"))
              ## axis(4, tics*scalar, tics, las = 2)
            }
          if(!real.dat)
            {
              hist(out$carc.all$roaddist, breaks = breaks, freq = T, bty = "n", col = "black", ylim = c(0, N/12),
                   xlab = "km from road", main = paste(N,"all carcasses,",sims[jj,"n.in.eff"],"within",btrack,"of surv eff"))
              if(gamma.y) curve(d(temp.gamm)(x)*N/(length(breaks)-1), add = T, col = "purple", lwd = 2)
              if(!gamma.y) segments(0, 1/max.dist, max.dist, 1/max.dist,col = "purple", lwd = 2)
              ymax <- par("usr")[4]
            }else{
              if(gamma.y)
                {
                  curve(d(temp.gamm)(x)*N/(length(breaks)-1), add = F, col = "purple", lwd = 2, xlab = "km",
                        ylab = "[y]")
                }else{
                  plot(0,0, type = "n", xlim = c(0, max.dist), ylim = c(0,1.5), xlab = "km",
                       ylab = "[y]")
                  segments(0, 1/max.dist, max.dist, 1/max.dist,col = "purple", lwd = 2)
                }
            }
          if(real.dat) ymax <- nrow(out$carc)/5
          hist(out$carc$roaddist, breaks = breaks, freq = T, bty = "n", col = "yellow", ylim = c(0, ymax),
               xlab = "km from road", main = paste(nrow(out$carc),"detected carcasses"))
          mtext("likelihood", 3, outer = T, adj = .25)
          mtext("estimated detection prob", 3, outer = T, adj = .75)
          mtext(paste(c(round(pars.frame[jj,],2), round(sims[jj, c("ml.est","ml.l","ml.u")],0)),collapse = ", "), 1, outer = T)
        }
      for(jj in 1:nrow(pars.frame))
        {
          if(!real.dat)
            {
              print(table(outs[[jj]]$carc.all$roadwhich))
            }else{
              print(table(outs[[jj]]$carc$roadwhich))
            }
        }
      dev.off()
    }
  sims[,c(4:9,14)] <- round(sims[,c(4:9,14)],0)
  to.return <- list(sims = sims, pars = pars.frame, outs = outs, inits.frame = inits.frame, par.boot = par.boot)
  save(to.return, file = save.file)
  return(to.return)
}

## wrapper for multicore
optnll.mcl <- function(nc = 1, ...)
  {
    temp <- mclapply(nc, optnll, ...)
    for(ii in 1:length(nc))
      {
        if(ii==1)
          {
            out <- temp[[ii]]
          }else{
            out$sims <- rbind(out$sims, temp[[ii]]$sims)
            out$pars <- rbind(out$pars, temp[[ii]]$pars)
            out$outs <- rbind(out$outs, temp[[ii]]$outs)
            out$inits.frame <- rbind(out$inits.frame, temp[[ii]]$inits.frame)
            out$par.boot <- rbind(out$par.boot, temp[[ii]]$par.boot)
          }
      }
   return(out)
  }


## Get confint using true values of parameters for a bunch of
## simulations and then plot them
confinter <- function(out, cdr.output, ddrr.output, pars, fix.pars=NULL, covmat,  do.confint = T, numsim = 101,
                      gamma.y, temp.gamm, vectorized.quad = T, only.obs = F, st.est = T, max.dist, cue.trips,
                      lower = .025, upper = .975, print.int = 20, browse = F)
{
    if(browse) browser()
    ## concatenate parameters
    all.pars <- c(pars, fix.pars)
    n.obs <- nrow(out$carc)
    ## draw random parameters centered on mle with covmat frm hessian
    rand.pars <- rmnorm(numsim, mean = pars , varcov = covmat)
    if(length(fix.pars)>0)
    {             # if some parameters were fixed, keep them fixed
        for(ii in 1:length(fix.pars))
        {
            rand.pars <- cbind(rand.pars, rep(fix.pars[ii], numsim))
            colnames(rand.pars)[ncol(rand.pars)] <- names(fix.pars)[ii]
        }
    }
    ## calculate the quadrature at the fitted parameters
    if(st.est)
    {                                   # calculate the quadrature for the mle
        quads.mlest <- g.dk(ddrr.output = ddrr.output, out = out, yys = yy.seq,
                            only.obs = only.obs,
                            sig.v = exp(all.pars["lsig.v"]), sig.m = exp(all.pars["lsig.m"]),
                            sig25 = exp(all.pars["lsig25"]), d.shape= exp(all.pars["ld.shape"]),
                            gamma.y = gamma.y, temp.gamm = temp.gamm, max.dist = max.dist,
                            cue.trips = cue.trips,
                            verbose = F, browse = F)
        n.mlest <- n.obs/quads.mlest
        if(vectorized.quad)
        {
            ## use vectorized version of g.dk over parameters to calculate all the quadratures at once
            quads <- g.dk.array(ddrr.output = ddrr.output, out = out, yys = yy.seq,
                                only.obs = only.obs,
                                sig.v = exp(rand.pars[,1]), sig.m = exp(rand.pars[,2]),
                                sig25 = exp(rand.pars[,3]), d.shape= exp(rand.pars[,4]),
                                gamma.y = gamma.y, temp.gamm = temp.gamm, max.dist = max.dist,
                                cue.trips = cue.trips,
                                verbose = F, browse = F)
        }else{
            ## older version using loop, may not actually be any slower than above
            quads <- rep(NA, numsim)
            for(ii in 1:numsim)
            {
                if(ii%%print.int == 1) print(paste("Quadrature calculation", ii, "of", numsim))
                quads[ii] <- g.dk(ddrr.output = ddrr.output, out = out, yys = yy.seq,
                                  only.obs = only.obs,
                                  sig.v = exp(rand.pars[ii,1]), sig.m = exp(rand.pars[ii,2]),
                                  sig25 = exp(rand.pars[ii, 3]), d.shape= exp(rand.pars[ii,4]),
                                  gamma.y = gamma.y, temp.gamm = temp.gamm, max.dist = max.dist,
                                  cue.trips = cue.trips,
                                  verbose = F, browse = F)
            }
        }
        CI <- n.obs/quantile(quads, c(upper,lower))
        ## Use non st estimator
    }else{                              # calculate the horvitz-thompson indiv estimator at mle
        mcprs <- pr.fun.yspec(out, cdr.output = cdr.output, ddrr.output = ddrr.output,
                              sig.v = exp(all.pars["lsig.v"]), sig.m = exp(all.pars["lsig.m"]),
                              sig25 = exp(all.pars["lsig25"]), d.shape= exp(all.pars["ld.shape"]),
                              yys = yy.seq,       # quadrature subdivisions
                              gamma.y = gamma.y,
                              cue.trips = cue.trips,
                              temp.gamm = temp.gamm, browse = F)
        n.mlest <-  sum(1/mcprs)
        mcpr.rand <- rep(NA, numsim)
        for(ii in 1:numsim)
        {
            if(ii%%print.int == 1) print(paste("Quadrature calculation", ii, "of", numsim))
            mcpr.rand[ii] <- sum(1/ pr.fun.yspec(out, cdr.output, ddrr.output,
                                                 sig.v = exp(rand.pars[ii,1]), sig.m = exp(rand.pars[ii,2]),
                                                 sig25 = exp(rand.pars[ii, 3]), d.shape= exp(rand.pars[ii,4]),
                                                 yys = yy.seq,       # quadrature subdivisions
                                                 gamma.y = gamma.y,
                                                 cue.trips = cue.trips,
                                                 temp.gamm = temp.gamm, browse = F))
        }
                    CI <-  quantile(mcpr.rand, c(lower, upper))
    }
    output <- c(est = n.mlest, lower = CI[1], upper = CI[2])
    return(output)
}


library(abind)
## same as g.dk but can take data.frames of parameter values to do multiple quadrature calculations rapdily
g.dk.array <- function(ddrr.output,                 # output of ddrr
                       out,                         #
                       yys = yy.seq,       # quadrature subdivisions
                       sig.v.vec, sig.m.vec, sig25.vec, d.shape.vec, # must be the same length
                       only.obs = F,
                       gamma.y = F, temp.gamm, max.dist, cue.trips,
                       verbose = F, browse = F)
{
    if(browse) browser()
    depth <- length(sig.v.vec)
    pdet.allcues <- h.fun.v(xx = rep(rep(yys, each = length(allcues[1:3])),depth), cc = rep(rep(allcues[1:3], length(yys)), depth),
                            sig.v = rep(sig.v.vec, each = length(allcues[1:3])*length(yys)),
                            sig.m = rep(sig.m.vec, each = length(allcues[1:3])*length(yys)),
                            sig25 = rep(sig25.vec, each = length(allcues[1:3])*length(yys)),
                            d.shape = rep(d.shape.vec, each = length(allcues[1:3])*length(yys)))
    pdet.allcues <- array(pdet.allcues, c(length(allcues[1:3]), length(yys), depth))
    ## Calculate marginal probability of finding carcass for each trip type (t = 0, 1 , etc, or by time later on)
    pmarg.bytrip <- abind(lapply(1:depth, function(kkk) cue.trips %*% pdet.allcues[,,kkk]), along = 3)
    ## pmarg.bytrip <- cue.trips %*% pdet.allcues
    ddrr <- ddrr.output$ddrr
    tripmat <- ddrr.output$tripmat
    ## initialize gdrs (row is for each dd-rr combo, col is for subdivs)
    gdrs <- array(0, c(nrow(ddrr), length(yys), depth))
    left.trunc <- array(0, c(nrow(ddrr), length(yys), depth))
    for(ii in 1:ncol(tripmat))
    {
        if(verbose) print(ii)
        is.trip <- ddrr$numtrips >= ii
        which.trips <- tripmat[is.trip,ii] + 1
                                        # + 1 is because trip 0's correspond to 1st row of
                                        # pmarg.bytrip, need to generalize if add tod
        gdrs[is.trip,,] <- gdrs[is.trip,,] + (1-gdrs[is.trip,,]) * pmarg.bytrip[which.trips,,]
        ## Need to left truncate effort outside time window. To do
        ## this have a separate matrix that will save a snapshot of
        ## gdrs at the last trip before d.min (the cumulative
        ## probability of having detected a carcass between d.min)
        ## which will then be subtracted from gdrs at the end.
        befores <- as.numeric(ddrr$dd[is.trip]) + tripmat[is.trip,ii] < out$d.min #trips being considered are before
        if(ii<ncol(tripmat)) #if not at the max numtrips yet see if the next trip is after d.min
        {
            next.after <- as.numeric(ddrr$dd[is.trip]) + tripmat[is.trip,ii + 1] >= out$d.min #& it is the last trip before d.min
        }else{ #if at max numtrips, then nextafter is TRUE because we want to subtract all previous prob
            next.after <- TRUE # this way, if all trips were before d.min, left.trunc becomes the full marg prob
        }
        next.after[is.na(next.after)] <- TRUE
        if(sum(is.trip)==1)             #separate if statement if it is a vector not a matrix
        {
            left.trunc[is.trip,,][befores & next.after,] <- gdrs[is.trip,,][befores & next.after,]
        }else{
            left.trunc[is.trip,,][befores & next.after,,] <- gdrs[is.trip,,][befores & next.after,,]
        }
    }
    gdrs <- gdrs - left.trunc
    ## multiply by [y] and sum rows so we now have sum across dy for each ddrr combo
    ## using rectangular quadrature built around interval midpoints, so find pdf in each interval
    if(gamma.y)
    {                                   # use truncated gamma
        y.pdf <- d(temp.gamm)(yy.seq)
        ## y.pdf <- p(temp.gamm)(yy.seq + int/2) - p(temp.gamm)(yy.seq - int/2)
    }else{                          # use uniform
        y.pdf <- dunif(yy.seq, 0, max.dist)
        ## y.pdf <- punif(yy.seq + int/2, 0, max.dist) - punif(yy.seq - int/2, 0, max.dist)
    }
    gdrs <- abind(lapply(1:depth, function(kkk) gdrs[,,kkk] %*% y.pdf ), along = 3)
    ## for each dd-rr combo, multiply by [dd][rr]
    dd.p <- ddrr$dd.p
    rr.p <- ddrr$rr.p
    ddrr.p <- dd.p*rr.p
    if(only.obs)                        # normalize by rr-dd in window
    {
        ddrr.p <- ddrr.p / sum(ddrr.p[ddrr$numtrips>0])
    }
    # i can only do this because gdrs dd-rr organization matches that
    # of ddrr. and where gdrs = 0, it doesn't matter that we're
    # normalizing. so sum(ddrr.p) > 0 when there are trips included
    # that didn't have any effort in ddrr, but many are already
    # deleted since empty roads don't show up from non-factor usage of
    # xtabs above
    gdrs <- abind(lapply(1:depth, function(kkk) gdrs[,,kkk] * ddrr.p), along = 2)
    quadsum <- colSums(gdrs)*int
    ## quadsum <- sum(gdrs)*int
}

