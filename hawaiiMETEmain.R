## load needed packages
devtools::load_all('~/Dropbox/Research/meteR')

setwd('~/Dropbox/Research/hawaiiMETE')

## read cleaned data from Gruner 2007 Syst Biol
x <- read.csv('~/Research/data/Gruner/gruner_clean.csv')


## data.frames holding all unique combinations of trophic, site and tree factors
byTrophBySite <- unique(x[, c('trophic', 'Site')])
byTrophBySite <- byTrophBySite[!(byTrophBySite$trophic %in% c('U', 'T')), ]

byTrophBySiteByTree <- unique(x[, c('trophic','Site', 'Tree')])
byTrophBySiteByTree <- byTrophBySiteByTree[!(byTrophBySiteByTree$trophic %in% c('U', 'T')), ]

bySite <- unique(x$Site)

bySiteByTree <- unique(x[, c('Site', 'Tree')])

## mete objects divided by trophic and site
mete.byTS <- apply(byTrophBySite, 1, function(s) {
    dat <- x[x$trophic==s[1] & x$Site==s[2], ]
    esf <- meteESF(dat$SpeciesCode, dat$Abundance, dat$IND_BIOM^0.75)
    out <- list(sad=sad(esf), ipd=ipd(esf))
    return(out)
})

## mete objects divided by trophic, site and tree
mete.byTST <- apply(byTrophBySiteByTree, 1, function(s) {
    dat <- x[x$trophic==s[1] & x$Site==s[2] & x$Tree==as.numeric(s[3]), ]
    this.n <- sum(dat$Abundance)
    if(sum(dat$Abundance) > 50) {
        esf <- meteESF(dat$SpeciesCode, dat$Abundance, dat$IND_BIOM^0.75)
        out <- list(sad=sad(esf), ipd=ipd(esf))  
    } else {
        out <- NULL
    }  
    return(out)
})

## save mete objects
save(mete.byTS, byTrophBySite, byTrophBySiteByTree, x, file='mete.byTS.RData')
