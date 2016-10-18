## load needed packages
library(meteR)
library(parallel)

setwd('~/Dropbox/Research/hawaiiMETE')

## read cleaned data from Gruner 2007 Syst Biol
grun <- read.csv('~/Research/data/Gruner/gruner_clean.csv')


## data.frames holding all unique combinations of trophic, site and tree factors
byTrophBySite <- unique(grun[, c('trophic', 'Site')])
byTrophBySite <- byTrophBySite[!(byTrophBySite$trophic %in% c('U', 'T')), ]
row.names(byTrophBySite) <- NULL


## mete objects divided by trophic and site
mete.byTS <- mclapply(1:nrow(byTrophBySite), mc.cores = 6, FUN = function(i) {
    dat <- grun[grun$trophic == byTrophBySite[i, 1] & grun$Site == byTrophBySite[i, 2], ]
    esf <- meteESF(dat$SpeciesCode, dat$Abundance, dat$IND_BIOM^0.75)
    out <- list(sad = sad(esf), ipd = ipd(esf))
    return(out)
})


## loop over METE objects and calculate z scores
zz <- mclapply(mete.byTS, mc.cores = 6, FUN = function(x) {
    sad.z <- logLikZ(x$sad, nrep = 999)$z
    ipd.z <- logLikZ(x$ipd, nrep = 999)$z
    
    return(c(sad.z = sad.z, ipd.z = ipd.z))
})

## combine z scores with summary table
byTrophBySite <- cbind(byTrophBySite, do.call(rbind, zz))

## save mete and summary objects
save(mete.byTS, byTrophBySite, grun, file = 'mete.byTS.RData')
