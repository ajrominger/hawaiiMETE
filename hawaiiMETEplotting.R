## load needed packages
devtools::load_all('~/Dropbox/Research/meteR')

setwd('~/Dropbox/Research/hawaiiMETE')

source('mete.byTS.RData')

byTrophBySite <- cbind(byTrophBySite,
                       t(sapply(mete.byTS, function(x) {
                           c(x$sad$state.var, maxN=max(x$sad$data), maxE=max(x$ipd$data))
                       })))


## plot SAD
pdf(file='fig_allSAD.pdf', width=6, height=4)

par(mar=rep(0.1, 4), oma=c(4, 4, 2, 2)+0.1)
layout(matrix(1:(3*5), nrow=3, byrow=TRUE))

for(troph in c('D', 'H', 'P')) {
    for(site in c('VO', 'LA', 'KH', 'MO', 'KA')) {
        plot(mete.byTS[[which(byTrophBySite$trophic==troph & 
                                  byTrophBySite$Site==site)]]$sad, 
             ptype='rad', add.legend=FALSE,
             log='y', yaxt=ifelse(site=='VO', 's', 'n'), xaxt=ifelse(troph=='P', 's', 'n'), 
             xlim=c(1, max(byTrophBySite$S0)), ylim=c(1, max(byTrophBySite$maxN)))
        
        if(site=='KA') mtext(troph, side=4, line=1)
        if(troph=='D') mtext(site, side=3, line=1)
    }
}

dev.off()


## plot IPD
pdf(file='fig_allIPD.pdf', width=6, height=4)

par(mar=rep(0.1, 4), oma=c(4, 4, 2, 2)+0.1)
layout(matrix(1:(3*5), nrow=3, byrow=TRUE))

for(troph in c('D', 'H', 'P')) {
    for(site in c('VO', 'LA', 'KH', 'MO', 'KA')) {
        plot(mete.byTS[[which(byTrophBySite$trophic==troph & 
                                  byTrophBySite$Site==site)]]$ipd, 
             ptype='rad', add.legend=FALSE,
             log='y', yaxt=ifelse(site=='VO', 's', 'n'), xaxt=ifelse(troph=='P', 's', 'n'), 
             xlim=c(1, max(byTrophBySite$S0)), ylim=c(1, max(byTrophBySite$maxN)))
        
        if(site=='KA') mtext(troph, side=4, line=1)
        if(troph=='D') mtext(site, side=3, line=1)
    }
}

dev.off()

