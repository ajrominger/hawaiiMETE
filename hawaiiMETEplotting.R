## ===========================================
## set-up
## ===========================================

devtools::load_all('~/Dropbox/Research/meteR')

setwd('~/Dropbox/Research/hawaiiMETE')

site.info <- read.csv('~/Research/data/Gruner/guner_site.csv')

load('mete.byTS.RData')
source('~/R_functions/logAxis.R')



## ===========================================
## organize all info into byTrophBySite
## ===========================================

## state vars
byTrophBySite <- cbind(byTrophBySite,
                       t(sapply(mete.byTS, function(x) {
                           c(x$sad$state.var, maxN=max(x$sad$data), maxE=max(x$ipd$data))
                       })))

## ages
byTrophBySite$siteAge <- site.info$age[match(byTrophBySite$Site, site.info$code)]
byTrophBySite$islandAge <- site.info$island.age[match(byTrophBySite$Site, site.info$code)]

## proportion non-native
byTrophBySite <- cbind(byTrophBySite, t(apply(byTrophBySite[, 1:2], 1, function(m) {
    this.match <- x[paste(m['trophic'], m['Site']) == paste(x$trophic, x$Site), , drop=FALSE]
    Nnon <- 1 - sum((this.match$Origin == 'endemic' | this.match$Origin == 'indigenous')) / nrow(this.match)
    Snon <- 1 - length(unique(this.match$SpeciesCode[this.match$Origin == 'endemic' | this.match$Origin == 'indigenous'])) / 
        length(unique(this.match$SpeciesCode))
    Enon <- 1 - sum(this.match$IND_BIOM[this.match$Origin == 'endemic' | this.match$Origin == 'indigenous']^0.75) /
        sum(this.match$IND_BIOM^0.75)
    
    return(c(Nnon, Snon, Enon))
})))

names(byTrophBySite)[10:12] <- c('N.none', 'S.none', 'E.none')

## logLik and associated Z for SAD
sad.z <- lapply(mete.byTS, function(x) logLikZ(x$sad, nrep=999, return.sim=TRUE, type='cumulative'))
byTrophBySite$sad.z <- sapply(sad.z, function(x) x$z)
byTrophBySite$sad.ll <- sapply(sad.z, function(x) x$obs)


## logLik and associated Z for IPD
ipd.z <- lapply(mete.byTS, function(x) logLikZ(x$ipd, nrep=999, return.sim=TRUE, type='cumulative'))
byTrophBySite$ipd.z <- sapply(ipd.z, function(x) x$z)
byTrophBySite$ipd.ll <- sapply(ipd.z, function(x) x$obs)

## ===========================================
## plot SAD
## ===========================================

## SAD rank plots
pdf(file='fig_allSAD.pdf', width=6, height=4)

par(mar=rep(0.1, 4), oma=c(4, 4, 2, 2)+0.1)
layout(matrix(1:(3*5), nrow=3, byrow=TRUE))

for(troph in c('D', 'H', 'P')) {
    for(site in c('VO', 'LA', 'KH', 'MO', 'KA')) {
        plot(mete.byTS[[which(byTrophBySite$trophic==troph & 
                                  byTrophBySite$Site==site)]]$sad, 
             ptype='cdf', add.legend=FALSE,
             yaxt=ifelse(site=='VO', 's', 'n'), xaxt=ifelse(troph=='P', 's', 'n'), 
#              log='y', 
             log='x', 
#             xlim=c(1, max(byTrophBySite$S0)),
             xlim=c(1, max(byTrophBySite$maxN)),
#             ylim=c(1, max(byTrophBySite$maxN))
             ylim=0:1
             )
        
        if(site=='KA') mtext(troph, side=4, line=1)
        if(troph=='D') mtext(site, side=3, line=1)
    }
}

dev.off()


## plot z-score for SAD

pdf('fig_sad_zscore.pdf', width=5, height=5)  # z-score by age

layout(matrix(1:3, nrow=3))
par(mar=rep(0.1, 4), oma=c(4, 4, 0, 2)+0.1)

for(troph in c('D', 'H', 'P')) {
    dat <- cbind(byTrophBySite$siteAge[byTrophBySite$trophic==troph], 
                 byTrophBySite$sad.z[byTrophBySite$trophic==troph])
    dat <- dat[order(dat[, 1]), ]
    plot(dat, type='b', log='x', 
         ylim=c(0, 2.25), 
         xaxt='n')
    abline(h=2, lty=2, col='gray')
    mtext(troph, side=4, line=1)
}
logAxis(1)
axis(1, at=5)

dev.off()

## z by invaded 

palette(c(D = 'black', H = 'green', P = 'blue'))

with(byTrophBySite, plot(N.none, sad.z, col=trophic))

with(byTrophBySite, plot(E.none, ipd.z, col=trophic))


## ===========================================
## plot IPD
## ===========================================

## IPD rank plots
pdf(file='fig_allIPD.pdf', width=6, height=4)

par(mar=rep(0.1, 4), oma=c(4, 4, 2, 2)+0.1)
layout(matrix(1:(3*5), nrow=3, byrow=TRUE))

for(troph in c('D', 'H', 'P')) {
    for(site in c('VO', 'LA', 'KH', 'MO', 'KA')) {
        plot(mete.byTS[[which(byTrophBySite$trophic==troph & 
                                  byTrophBySite$Site==site)]]$ipd, 
             ptype='rad', add.legend=FALSE,
             log='xy', yaxt=ifelse(site=='VO', 's', 'n'), xaxt=ifelse(troph=='P', 's', 'n'), 
             xlim=c(1, max(byTrophBySite$N0)), ylim=c(1, max(byTrophBySite$maxN)))
        
        if(site=='KA') mtext(troph, side=4, line=1)
        if(troph=='D') mtext(site, side=3, line=1)
    }
}

dev.off()


