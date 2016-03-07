## ===========================================
## set-up
## ===========================================

devtools::load_all('~/Dropbox/Research/meteR')
library(vegan)

setwd('~/Dropbox/Research/hawaiiMETE')

site.info <- read.csv('~/Research/data/Gruner/guner_site.csv')

load('mete.byTS.RData')
source('~/R_functions/logAxis.R')



## ===========================================
## function organize all info into byTrophBySite or similar data.frame
## ===========================================

composeData <- function(dat, mete) {
    ## state vars
    dat <- cbind(dat, 
                 t(sapply(mete, function(x) {
                     c(x$sad$state.var, maxN=max(x$sad$data), maxE=max(x$ipd$data))
                 })))
    
    ## ages
    dat$siteAge <- site.info$age[match(dat$Site, site.info$code)]
    dat$islandAge <- site.info$island.age[match(dat$Site, site.info$code)]
    
    ## prop non-native
    if('trophic' %in% names(dat)) {
        matchString <- function(m) paste(m['trophic'], m['Site'])
        matchTable <- paste(x$trophic, x$Site)
    } else {
        matchString <- function(m) paste(m['Site'])
        matchTable <- x$Site
    }
    
    if('Tree' %in% names(dat)) {
        matchStringFun <- function(m) paste(matchString(m), m['Tree'])
        matchTable <- paste(matchTable, x$Tree)
    } else {
        matchStringFun <- matchString
    }
    
    dat <- cbind(dat, t(apply(dat, 1, function(m) {
        this.match <- x[matchStringFun(m) == matchTable, , drop=FALSE]
        Nnon <- 1 - sum((this.match$Origin == 'endemic' | this.match$Origin == 'indigenous')) / nrow(this.match)
        Snon <- 1 - length(unique(this.match$SpeciesCode[this.match$Origin == 'endemic' | this.match$Origin == 'indigenous'])) / 
            length(unique(this.match$SpeciesCode))
        Enon <- 1 - sum(this.match$IND_BIOM[this.match$Origin == 'endemic' | this.match$Origin == 'indigenous']^0.75) /
            sum(this.match$IND_BIOM^0.75)
        
        return(c(N.none=Nnon, S.none=Snon, E.none=Enon))
    })))
    
    
    ## normalized beta-diversity
    if('trophic' %in% names(dat)) {
        matchString <- function(m) paste(m['trophic'], m['Site'])
        matchTable <- paste(x$trophic, x$Site)
    } else {
        matchString <- function(m) paste(m['Site'])
        matchTable <- x$Site
    }
    
    dat$bdiv <- apply(dat, 1, function(m) {
        this.match <- x[matchString(m) == matchTable, , drop=FALSE]
        mat <- samp2sitespp(this.match$Tree, this.match$SpeciesCode, this.match$Abundance)
        
        if('Tree' %in% names(dat)) {
            obs <- NA
        } else {
            obs <- mean(vegdist(mat, 'chao'))
        }
        
        
        sim <- replicate(999, {
            if('Tree' %in% names(dat)) {
                return(NA)
            } else {
                newdat <- this.match
                newdat$Tree <- sample(newdat$Tree)
                mat <- samp2sitespp(newdat$Tree, newdat$SpeciesCode, newdat$Abundance)
                return(mean(vegdist(mat, 'chao')))
            }
        })
        sim <- c(sim, obs)
        
        return((obs - mean(sim))/sd(sim))
    })
    
    ## logLik and associated Z for SAD
    sad.z <- lapply(mete, function(x) logLikZ(x$sad, nrep=999, return.sim=FALSE, type='cumulative'))
    dat$sad.z <- sapply(sad.z, function(x) x$z)
    dat$sad.ll <- sapply(sad.z, function(x) x$obs)
    
    
    ## logLik and associated Z for IPD
    ipd.z <- lapply(mete, function(x) logLikZ(x$ipd, nrep=999, return.sim=FALSE, type='cumulative'))
    dat$ipd.z <- sapply(ipd.z, function(x) x$z)
    dat$ipd.ll <- sapply(ipd.z, function(x) x$obs)
    
    return(dat)
}

## convenience function for making site by spp matrix
samp2sitespp <- function(site, spp, abund) {
    x <- tapply(abund, list(site=site, spp=spp), sum)
    x[is.na(x)] <- 0
    
    return(x)
}


## orgainize byTrophBySite
byTrophBySite <- composeData(byTrophBySite, mete.byTS)

## explore beta-diversity
with(byTrophBySite, plot(bdiv, sad.z, col=trophic))


## same for bySite
bySite <- composeData(bySite, mete.byS)
with(bySite, plot(bdiv, sad.z))
with(bySite, plot(siteAge, sad.z, log='x'))

plot(mete.byS[[which(bySite$Site=='KA')]]$sad, ptype='rad', log='y')


## ===========================================
## byTrophBySiteByTree shows pretty much everything 
## is METEish and doesn't change with age...frankly 
## going to ignore cause can't rule out small sample 
## sizes
## ===========================================


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


