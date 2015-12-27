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




## ===========================================
## use byTrophBySiteByTree to get beta-div
## ===========================================

