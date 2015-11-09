## load needed packages
devtools::load_all('~/Dropbox/Research/meteR')

setwd('~/Dropbox/Research/hawaiiMETE')

## read cleaned data from Gruner 2007 Syst Biol
x <- read.csv('~/Research/data/Gruner/gruner_clean.csv')
site.info <- read.csv('~/Research/data/Gruner/guner_site.csv')

## vectors holding all unique combinations of trophic, site and tree factors
byTrophBySite <- unique(x[, c('trophic', 'Site')])
byTrophBySite <- byTrophBySite[!(byTrophBySite$trophic %in% c('U', 'T')), ]

byTrophBySiteByTree <- unique(x[, c('trophic','Site', 'Tree')])
byTrophBySiteByTree <- byTrophBySiteByTree[!(byTrophBySiteByTree$trophic %in% c('U', 'T')), ]

## mete objects divided by trophic and site
mete.byTS <- apply(byTrophBySite, 1, function(s) {
    dat <- x[x$trophic==s[1] & x$Site==s[2], ]
    esf <- meteESF(dat$SpeciesCode, dat$Abundance, dat$IND_BIOM^0.75)
    out <- list(sad=sad(esf), ipd=ipd(esf))
    return(out)
})


## plot SAD
pdf(file='fig_allSAD.pdf', width=5, height=5)
par(mar=rep(0.1, 4), oma=c(4, 4, 0, 0)+0.1)
layout(matrix(1:(3*5), nrow=3, byrow=TRUE))

for(troph in c('D', 'H', 'P')) {
    for(site in c('VO', 'LA', 'KH', 'MO', 'KA')) {
        plot(mete.byTS[[which(byTrophBySite$trophic==troph & 
                                  byTrophBySite$Site==site)]]$sad, 
             ptype='rad', log='y', xaxt='n', yaxt='n', add.legend=FALSE)
        legend('topright', legend=paste(site, troph))
    }
}

dev.off()


## plot IPD
pdf(file='fig_allIPD.pdf', width=5, height=5)
par(mar=rep(0.1, 4), oma=c(4, 4, 0, 0)+0.1)
layout(matrix(1:(3*5), nrow=3, byrow=TRUE))

for(troph in c('D', 'H', 'P')) {
    for(site in c('VO', 'LA', 'KH', 'MO', 'KA')) {
        plot(mete.byTS[[which(byTrophBySite$trophic==troph & 
                                  byTrophBySite$Site==site)]]$ipd, 
             ptype='rad', log='y', xaxt='n', yaxt='n', add.legend=FALSE)
        legend('topright', legend=paste(site, troph))
    }
}

dev.off()

