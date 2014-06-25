###################################################################### 
## Bellan et al. 2012 - A hierarchical distance sampling approach to
## estimating mortality rates from opportunistic carcass surveillance
## data
######################################################################
## R scripts for analysis

## Set local directory
## setwd("local dir")
rm(list=ls())
library(abind)
library(multicore)                      # for parallelization if working on the amazon cloud
library(MASS)
library(mnormt)
source("Likelihood and Simulation Functions.R")                    # source likelihood functions
on.cloud <- T


if(on.cloud)                            #on amazon cloud?
  {
    args=(commandArgs(TRUE))
    ## args is now a list of character vectors
    ## First check to see if arguments are passed.
    ## Then cycle through each element of the list and evaluate the expressions.
    if(length(args)>0)
      for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
      }
    ## if arguments not supplied, set defaults
    ##  which country
    aaaa <- ""               #need at least one object so ls() isn't empty
    if(sum(ls() %in% "scen")==0)  scen <- 1 # default scenario 1
    if(sum(ls() %in% "term.on.finish")==0)  term.on.finish <- F # default scenario 1
    if(sum(ls() %in% "nboot")==0) nboot <- 1000
    if(sum(ls() %in% "quantsims")==0) quantsims <- 1000
    if(sum(ls() %in% "num.cores")==0) num.cores <- 8 # on the cloud
    if(sum(ls() %in% "nn")==0) nn <- 1
    if(sum(ls() %in% "verbose")==0) verbose <- T
  }

## to run from batch on Amazon Cloud
## nohup R CMD BATCH '--args scen=1 term.on.finish=T nboot=1000 quantsims=1000 num.cores=8 nn=13 verbose=T' runAnalysis.R &

## Initialize simulation parameters
max.dist <- .8                          # maximum half strip width
subdiv <- 100                           # number of subdivisions in rectangular quadrature
int <- max.dist / subdiv                # quadrature interal size
yy.seq <- seq(0 + int/2, max.dist - int/2, by = int) # midpoints of quadrature interval
rrs <- unique(out.real$carc$roadwhich)   #all roads with a carcass
d.min <- min(out.real$carc$nday)         #first day of time period
d.max <- max(out.real$carc$nday)         #last day of time period
## find all road-days with nonzero seff
dr.out <- ddrr.select(out.real, btrack = btrack, truncate.end = T)
ddrr.str <- paste(dr.out[["ddrr"]]$rr, dr.out[["ddrr"]]$dd) #for matching to n.in.eff below
ddrr.str <- ddrr.str[dr.out[["ddrr"]]$numtrips>0]
## folder to save all output
save.folder <- ""
true.vals <- c(sig.v = .4, sig.m = .15, sig25 = .1, d.shape = 2) # det fnx pars
start.time1 <- Sys.time()

## ######################################################################
## ## TEST Run: Only 1 simulation, 2 non-parametric bootstrap resamples,
## ## and 100 parametric bootstrap resamples, shouldn't take long. Look
## ## at output pdf file. All arguments are described in the "freqQuad.R"
## ## file
## ######################################################################
## stuff <- optnll.mcl(nc = 1,             # cores to run on
##                     true.vals = true.vals, which.pars = 1:3,
##                     nn = 1, only.obs = TRUE, st.like = FALSE, st.est = FALSE, only.in.eff = TRUE, ddrr.str = ddrr.str,
##                     real.dat = FALSE,
##                     hess.conf = T, boot = T, nboot = 2, boot.y = T, boot.c = T,
##                     quantsims = 2, like.pts = TRUE, max.dist = max.dist,
##                     gamma.y = FALSE, g.shape = 1.2, g.scale = .8, vectorized.quad = FALSE,
##                     make.pdf = TRUE, d.min = d.min, d.max = d.max, rrs = rrs, N = 300,
##                     rand.dd.p = FALSE, rand.rr.p = FALSE,
##                     file.name = paste(save.folder, "test", sep = ""),
##                     sd.pert = .5,  do.sann = T, hess = T, nm.its = 500, tr = 0, sann.its = 40, reltol = 10^-7,
##                     show.pars = F, verbose = T, timer = F, browse =F)



######################################################################
######################################################################
######################################################################


######################################################################
## Simulations for Biometrics Manuscript
######################################################################

if(on.cloud)
  {
    if(scen == 1)
      {
######################################################################
        ## Scenario 1
        ## Uniform [y], No ST het
######################################################################
        true.vals <- c(sig.v = .4, sig.m = .15, sig25 = .1, d.shape = 2) # det fnx pars
        stuff <- optnll.mcl(nc = 1:num.cores, true.vals = true.vals, which.pars = 1:3,
                        nn = nn, only.obs = TRUE, st.like = FALSE, st.est = FALSE, only.in.eff = TRUE, ddrr.str = ddrr.str,
                        real.dat = FALSE,
                        hess.conf = T, boot = T, nboot = nboot,
                        quantsims = quantsims, like.pts = TRUE, max.dist = max.dist,
                        gamma.y = FALSE, g.shape = 1.2, g.scale = .8, vectorized.quad = FALSE,
                        make.pdf = TRUE, d.min = d.min, d.max = d.max, rrs = rrs, N = 300,
                        rand.dd.p = FALSE, rand.rr.p = FALSE,
                        file.name = paste(save.folder, "biom yuni nst", sep = ""),
                        sd.pert = .5,  do.sann = T, hess = T, nm.its = 500, tr = 0, sann.its = 40, reltol = 10^-7,
                        show.pars = F, verbose = verbose, timer = F, browse =F)
      }
    if(scen == 2)
      {
######################################################################
        ## Scenario 2
        ## Uniform [y], ST het
######################################################################
        true.vals <- c(sig.v = .4, sig.m = .15, sig25 = .1, d.shape = 2) # det fnx pars
        stuff <- optnll.mcl(nc = 1:num.cores, true.vals = true.vals, which.pars = 1:3,
                        nn = nn, only.obs = TRUE, st.like = FALSE, st.est = FALSE, only.in.eff = TRUE, ddrr.str = ddrr.str,
                        real.dat = FALSE,
                        hess.conf = T, boot = TRUE, nboot = nboot,
                        quantsims = quantsims, like.pts = TRUE, max.dist = max.dist,
                        gamma.y = FALSE, g.shape = .942, g.scale = 1.89, vectorized.quad = FALSE,
                        make.pdf = TRUE, d.min = d.min, d.max = d.max, rrs = rrs, N = 300,
                        rand.dd.p = TRUE, rand.rr.p = TRUE,
                        file.name = paste(save.folder, "biom yuni st", sep = ""),
                        sd.pert = .5,  do.sann = T, hess = T, nm.its = 500, tr = 0, sann.its = 40, reltol = 10^-7,
                        show.pars = F, verbose = verbose, timer = F, browse =F)
      }
    if(scen == 3)
      {
######################################################################
        ## Scenario 3
        ## GAMMA [y], No ST het
######################################################################
        true.vals <- c(sig.v = .4, sig.m = .15, sig25 = .1, d.shape = 2) # det fnx pars
        stuff <- optnll.mcl(nc = 1:num.cores, true.vals = true.vals, which.pars = 1:3,
                        nn = nn, only.obs = TRUE, st.like = FALSE, st.est = FALSE, only.in.eff = TRUE, ddrr.str = ddrr.str,
                        real.dat = FALSE,
                        hess.conf = T, boot = TRUE, nboot = nboot,
                        quantsims = quantsims, like.pts = TRUE, max.dist = max.dist,
                        gamma.y = TRUE, g.shape = .942, g.scale = 1.89, vectorized.quad = FALSE,
                        make.pdf = TRUE, d.min = d.min, d.max = d.max, rrs = rrs, N = 300,
                        rand.dd.p = FALSE, rand.rr.p = FALSE,
                        file.name = paste(save.folder, "biom ygam nst", sep = ""),
                        sd.pert = .5,  do.sann = T, hess = T, nm.its = 500, tr = 0, sann.its = 40, reltol = 10^-7,
                        show.pars = F, verbose = verbose, timer = F, browse =F)
      }
    if(scen == 4)
      {
######################################################################
        ## Scenario 4
        ## GAMMA [y], ST het
######################################################################
        true.vals <- c(sig.v = .4, sig.m = .15, sig25 = .1, d.shape = 2) # det fnx pars
        stuff <- optnll.mcl(nc = 1:num.cores, true.vals = true.vals, which.pars = 1:3,
                        nn = nn, only.obs = TRUE, st.like = FALSE, st.est = FALSE, only.in.eff = TRUE, ddrr.str = ddrr.str,
                        real.dat = FALSE,
                        hess.conf = T, boot = TRUE, nboot = nboot,
                        quantsims = quantsims, like.pts = TRUE, max.dist = max.dist,
                        gamma.y = TRUE, g.shape = .942, g.scale = 1.89, vectorized.quad = FALSE,
                        make.pdf = TRUE, d.min = d.min, d.max = d.max, rrs = rrs, N = 300,
                        rand.dd.p = TRUE, rand.rr.p = TRUE,
                        file.name = paste(save.folder, "biom ygam st", sep = ""),
                        sd.pert = .5,  do.sann = T, hess = T, nm.its = 500, tr = 0, sann.its = 40, reltol = 10^-7,
                        show.pars = F, verbose = verbose, timer = F, browse =F)
      }
    hours <- round(as.numeric(difftime(Sys.time(), start.time1, unit = "hour")),3)
    save(hours, file = paste("took",hours,"hrs.Rdata"))
    save.image(file="workspace.Rdata")      #resave workspace
    if(getwd()=="/home/ubuntu/distsamp")       # if on cloud
      {
        ## Upload all output to S3 Bucket
        dirnm <- paste("Scen", scen, "-", format(Sys.time(), "%Y%m%d-%H:%M"), sep = "")
        comnd <- paste("s3cmd put /home/ubuntu/distsamp s3://distsamp/",dirnm,"/ --recursive", sep ="")
        system(comnd)
        if(term.on.finish)
          {
            ## Shutdown instance in 15 minutes, give time for copying and for email alert
            system("sudo shutdown -h 15")
          }
      }
  }

if(!on.cloud)
  {
######################################################################
    ## Table 1: Results of Simulations
######################################################################
    load.folder <- ""
    nms <- c("biom yuni nst.Rdata", "biom yuni st.Rdata", "biom ygam nst.Rdata", "biom ygam st.Rdata")
    tab1 <- data.frame(uniSt = rep(NA, 13), uniNst=NA, gamNst=NA, gamSt=NA)
    sef <- function(x) {sd(x)/sqrt(length(x))}
    for(ii in 1:4)
      {
        load(paste(load.folder, nms[ii], sep = ""))
        sims <- to.return$sims
        scen1 <- paste(round(apply(sims[,c(3:6,8:9)], 2, mean)), " (",
                       round(apply(sims[,c(3:6,8:9)], 2, sef)), ")", sep = "")
        scen1 <- c(scen1[1:2], paste(scen1[5],"-", scen1[6]), paste(scen1[3],"-", scen1[4]))
        scen1 <-  c(scen1, paste(signif(apply(to.return$pars[,c(1:3)], 2, mean), 2), " (",
                                 signif(apply(to.return$pars[,c(1:3)], 2, sef), 2), ")", sep = ""))
        mse <- mean((sims$n.in.eff - sims$ml.est)^2)
        bias <- mean(sims$ml.est - sims$n.in.eff)
        pval.hess <- mean(sims$ml.l <= sims$n.in.eff & sims$ml.u >= sims$n.in.eff)
        pval.boot <- mean(sims$ml.l.b <= sims$n.in.eff & sims$ml.u.b >= sims$n.in.eff)
        ## look at pval with parametric lower bound & nonparametric upper bound
                                        #    pval.boothess <- mean(sims$ml.l <= sims$n.in.eff & sims$ml.u.b >= sims$n.in.eff)    
        scen1 <- c(scen1, bias, mse, pval.hess, pval.boot, pval.boothess)
        tab1[1:length(scen1),ii] <- scen1
      }
    tab1 <- tab1[c(1:2,4,3,5:nrow(tab1)),]
    tab1
    write.csv(tab1, paste(save.folder, "simresults.csv", sep = ""))
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
    ## Estimate real data
######################################################################
    load("realcarc120105.Rdata")
    true.vals <- c(sig.v = .4, sig.m = .12, sig25 = .1, d.shape = 2) # det fnx pars
    quantsims <- 1
    nboot <- 1000
    boot <- T

    ## [y] unif, 10,000 parametric bootstraps
    stuff1 <- optnll(true.vals = true.vals, which.pars = 1:3,
                     nn = 1, only.obs = TRUE, st.like = FALSE, st.est = FALSE,
                     real.dat = TRUE, out.real = out.real,
                     hess.conf = T, boot = boot, nboot = 2, boot.y = TRUE, #boot over y
                     quantsims = 10000, like.pts = TRUE, max.dist = max.dist,
                     gamma.y = F, g.shape = .942, g.scale = 1.89, vectorized.quad = FALSE,
                     make.pdf = TRUE, d.min = d.min, d.max = d.max, rrs = rrs, N = 300,
                     rand.dd.p = FALSE, rand.rr.p = FALSE, #st het!!
                     file.name = paste(save.folder, "biom real ygam quants", sep = ""),
                     sd.pert = .5,  do.sann = T, hess = T, nm.its = 500, tr = 0, sann.its = 40, reltol = 10^-7,
                     show.pars = F, verbose = T, timer = F, browse =F)

    ## [y] uniform, 1000 non-parametric bootstraps
    stuff1 <- optnll(true.vals = true.vals, which.pars = 1:3,
                     nn = 1, only.obs = TRUE, st.like = FALSE, st.est = FALSE,
                     real.dat = TRUE, out.real = out.real,
                     hess.conf = T, boot = boot, nboot = nboot, boot.y = FALSE, boot.c = F, #normal boot
                     quantsims = quantsims, like.pts = TRUE, max.dist = max.dist,
                     gamma.y = FALSE, g.shape = .942, g.scale = 1.89, vectorized.quad = FALSE,
                     make.pdf = TRUE, d.min = d.min, d.max = d.max, rrs = rrs, N = 300,
                     rand.dd.p = FALSE, rand.rr.p = FALSE, #st het!!
                     file.name = paste(save.folder, "biom real yuni boot1", sep = ""),
                     sd.pert = .5,  do.sann = T, hess = T, nm.its = 500, tr = 0, sann.its = 40, reltol = 10^-7,
                     show.pars = F, verbose = verbose, timer = F, browse =F)


    ## [y] uniform, 1000 non-parametric bootstraps, resampling [c|t]
    stuff1 <- optnll(true.vals = true.vals, which.pars = 1:3,
                     nn = 1, only.obs = TRUE, st.like = FALSE, st.est = FALSE,
                     real.dat = TRUE, out.real = out.real,
                     hess.conf = T, boot = boot, nboot = nboot, boot.y = TRUE, boot.c = T, #boot over y & c
                     quantsims = quantsims, like.pts = TRUE, max.dist = max.dist,
                     gamma.y = FALSE, g.shape = .942, g.scale = 1.89, vectorized.quad = FALSE,
                     make.pdf = TRUE, d.min = d.min, d.max = d.max, rrs = rrs, N = 300,
                     rand.dd.p = FALSE, rand.rr.p = FALSE, #st het!!
                     file.name = paste(save.folder, "biom real yuni bootc", sep = ""),
                     sd.pert = .5,  do.sann = T, hess = T, nm.its = 500, tr = 0, sann.its = 40, reltol = 10^-7,
                     show.pars = F, verbose = verbose, timer = F, browse =F)
    stuff1$s

    ## [y] gamma, 10,000 parametric bootstraps
    stuff1 <- optnll(true.vals = true.vals, which.pars = 1:3,
                     nn = 1, only.obs = TRUE, st.like = FALSE, st.est = FALSE,
                     real.dat = TRUE, out.real = out.real,
                     hess.conf = T, boot = boot, nboot = 2, boot.y = TRUE, #boot over y
                     quantsims = 10000, like.pts = TRUE, max.dist = max.dist,
                     gamma.y = TRUE, g.shape = .942, g.scale = 1.89, vectorized.quad = FALSE,
                     make.pdf = TRUE, d.min = d.min, d.max = d.max, rrs = rrs, N = 300,
                     rand.dd.p = FALSE, rand.rr.p = FALSE, #st het!!
                     file.name = paste(save.folder, "biom real ygam quants", sep = ""),
                     sd.pert = .5,  do.sann = T, hess = T, nm.its = 500, tr = 0, sann.its = 40, reltol = 10^-7,
                     show.pars = F, verbose = verbose, timer = F, browse =F)

    ## [y] gamma, 1000 non-parametric bootstraps
    stuff1 <- optnll(true.vals = true.vals, which.pars = 1:3,
                     nn = 1, only.obs = TRUE, st.like = FALSE, st.est = FALSE,
                     real.dat = TRUE, out.real = out.real,
                     hess.conf = T, boot = boot, nboot = nboot, boot.y = FALSE, #normal boot
                     quantsims = quantsims, like.pts = TRUE, max.dist = max.dist,
                     gamma.y = TRUE, g.shape = .942, g.scale = 1.89, vectorized.quad = FALSE,
                     make.pdf = TRUE, d.min = d.min, d.max = d.max, rrs = rrs, N = 300,
                     rand.dd.p = FALSE, rand.rr.p = FALSE, #st het!!
                     file.name = paste(save.folder, "biom real ygam boot1", sep = ""),
                     sd.pert = .5,  do.sann = T, hess = T, nm.its = 500, tr = 0, sann.its = 40, reltol = 10^-7,
                     show.pars = F, verbose = verbose, timer = F, browse =F)

    ## [y] gamma, 1000 non-parametric bootstraps, resampling over [y]
    stuff1 <- optnll(true.vals = true.vals, which.pars = 1:3,
                     nn = 1, only.obs = TRUE, st.like = FALSE, st.est = FALSE,
                     real.dat = TRUE, out.real = out.real,
                     hess.conf = T, boot = boot, nboot = nboot, boot.y = TRUE, #boot over y
                     quantsims = quantsims, like.pts = TRUE, max.dist = max.dist,
                     gamma.y = TRUE, g.shape = .942, g.scale = 1.89, vectorized.quad = FALSE,
                     make.pdf = TRUE, d.min = d.min, d.max = d.max, rrs = rrs, N = 300,
                     rand.dd.p = FALSE, rand.rr.p = FALSE, #st het!!
                     file.name = paste(save.folder, "biom real ygam booty", sep = ""),
                     sd.pert = .5,  do.sann = T, hess = T, nm.its = 500, tr = 0, sann.its = 40, reltol = 10^-7,
                     show.pars = F, verbose = verbose, timer = F, browse =F)


    ## [y] gamma, 1000 non-parametric bootstraps, resampling over [c|t]
    stuff1 <- optnll(true.vals = true.vals, which.pars = 1:3,
                     nn = 1, only.obs = TRUE, st.like = FALSE, st.est = FALSE,
                     real.dat = TRUE, out.real = out.real,
                     hess.conf = T, boot = boot, nboot = nboot, boot.y = F, boot.c = T, #boot over y
                     quantsims = quantsims, like.pts = TRUE, max.dist = max.dist,
                     gamma.y = TRUE, g.shape = .942, g.scale = 1.89, vectorized.quad = FALSE,
                     make.pdf = TRUE, d.min = d.min, d.max = d.max, rrs = rrs, N = 300,
                     rand.dd.p = FALSE, rand.rr.p = FALSE, #st het!!
                     file.name = paste(save.folder, "biom real ygam bootc", sep = ""),
                     sd.pert = .5,  do.sann = T, hess = T, nm.its = 500, tr = 0, sann.its = 40, reltol = 10^-7,
                     show.pars = F, verbose = verbose, timer = F, browse =F)
    stuff1$s

    ## [y] gamma, 1000 non-parametric bootstraps, resampling over [y] and [c|t]
    stuff1 <- optnll(true.vals = true.vals, which.pars = 1:3,
                     nn = 1, only.obs = TRUE, st.like = FALSE, st.est = FALSE,
                     real.dat = TRUE, out.real = out.real,
                     hess.conf = T, boot = boot, nboot = nboot, boot.y = TRUE, boot.c = T, #boot over y
                     quantsims = quantsims, like.pts = TRUE, max.dist = max.dist,
                     gamma.y = TRUE, g.shape = .942, g.scale = 1.89, vectorized.quad = FALSE,
                     make.pdf = TRUE, d.min = d.min, d.max = d.max, rrs = rrs, N = 300,
                     rand.dd.p = FALSE, rand.rr.p = FALSE, #st het!!
                     file.name = paste(save.folder, "biom real ygam bootyc", sep = ""),
                     sd.pert = .5,  do.sann = T, hess = T, nm.its = 500, tr = 0, sann.its = 40, reltol = 10^-7,
                     show.pars = F, verbose = verbose, timer = F, browse =F)

######################################################################
    ## Table 1: Analyses of Real Data
######################################################################
    load.folder <- save.folder
    all.fls <- list.files(load.folder)
    nms <- c("biom real yuni quants","biom real yuni boot1","biom real yuni booty","biom real yuni bootyc",
             "biom real ygam quants","biom real ygam boot1","biom real ygam booty","biom real ygam bootyc")
    for(ii in 1:length(nms))
      {
        load(paste(load.folder,nms[ii],".Rdata",sep=""))
        if(ii==1)
          {
            tab2 <- to.return$sims
            tab3 <- to.return$pars
          }else{
            tab2 <- rbind(tab2, to.return$sims)
            tab3 <- rbind(tab3, to.return$pars)
          }
      }
    tab2[,-c(1:2,14)]
    tab3

######################################################################
    ## Fig 2: Plot cue temporal distribution
######################################################################
    pdf(paste(save.folder, "cueplot.pdf", sep = ""), width = 3, height = 3)
    load(file = "cue.trips.Rdata")
    par(mar = c(4,4,.5,.5))
    plot(0,0, type = "n", xlim = c(0,5), ylim = c(0,.6), xlab="days since death", ylab = "proportion of daytime", bty = "n")
    for(ii in 2:4) lines(cue.trips[1:6,c(1,ii)], lty = ii-1)
    legend("topright", leg = c("avian", "mammalian", "fresh carcass"), lty = 1:3, bty = "n")
    dev.off()
######################################################################

######################################################################
    ## Figure 4: Sighting distance data & fitted det fxn (using gamma [y])
######################################################################
    ## [y] gamma, no interval estimation, just want sigma pars for plotting
    stuff1 <- optnll(true.vals = true.vals, which.pars = 1:3,
                     nn = 1, only.obs = TRUE, st.like = FALSE, st.est = FALSE,
                     real.dat = TRUE, out.real = out.real,
                     hess.conf = T, boot = boot, nboot = 2, boot.y = TRUE, #boot over y
                     quantsims = 2, like.pts = TRUE, max.dist = max.dist,
                     gamma.y = TRUE, g.shape = .942, g.scale = 1.89, vectorized.quad = FALSE,
                     make.pdf = TRUE, d.min = d.min, d.max = d.max, rrs = rrs, N = 300,
                     rand.dd.p = FALSE, rand.rr.p = FALSE, #st het!!
                     file.name = paste(save.folder, "biom real ygam quants", sep = ""),
                     sd.pert = .5,  do.sann = T, hess = T, nm.its = 500, tr = 0, sann.its = 40, reltol = 10^-7,
                     show.pars = F, verbose = verbose, timer = F, browse =F)
    pdf(paste(save.folder, "hist and det fxn.pdf", sep = ""), width = 7.5, height = 3)
    par(mfrow=c(1,3), mar = c(3, 3, 1.5, 3.5), oma = c(2,2,0,0))
    breaks <- seq(0,.8, by = .1)*1000
    sel <- out.real$carc$roaddist < .8
    ymax <- 21
    mains <- c("avian scavengers", "mammalian scavengers", "fresh carcass")
    for(ii in 1:3)
      {
        cc <- allcues[ii]
        h1 <- hist(1000*out.real$carc$roaddist[sel & out.real$carc$cue==cc], col = "gray", main = mains[ii],
                   breaks = breaks, ylim = c(0, ymax), border = NA, freq = T, xlab = "", ylab = "", yaxt = "n")
        area1 <- sum(h1$counts)*.1
        d.funx <- function(x) {sapply(x, function(x) as.numeric(exp( -(x/ (stuff1$pars[ii]) )^2) ))}
        area2 <- integrate(d.funx, 0, max.dist)$value
        scalar <- area1/area2
        if(ii==1) axis(2, at = seq(0,20,5))
                                        #        if(ii==3) axis(4, at = c(0,5,10), labels = c(0,.5,1), col = "red")
        print(max(h1$counts))
        xseq <- c(1:800)/1000
        yseq <- sapply(xseq, d.funx)
        lines(1000*xseq, scalar*yseq, col = "black")
      }
    mtext("meters from the road", side = 1, line = .5, outer = T)
    mtext("# carcasses detected", side = 2, line = .5, outer = T)
    dev.off()
    ## look at distribution of parameters
    parsf <- to.return$pars
    true.vals
    breaks <- seq(0,1,by=.01)
    par(mfrow = c(3,1))
    for(ii in 1:3)
      {
        hist(parsf[,ii], xlim = c(0,1), breaks = breaks, col = "black", xlab = names(true.vals)[ii])
        abline(v = true.vals[ii], col="red")
      }
  }                                     # analysis if not on cloud
