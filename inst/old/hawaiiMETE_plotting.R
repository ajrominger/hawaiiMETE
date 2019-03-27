library(meteR)

setwd('~/Dropbox/Research/hawaiiMETE')

load('mete.byTS.RData')

with(byTrophBySite[byTrophBySite$trophic == 'H', ], 
     plot(factor(as.character(Site), levels = c('VO', 'LA', 'KH', 'MO', 'KA')), 
          sad.z))

plot(byTrophBySite$S0, byTrophBySite$sad.z, col = factor(byTrophBySite$trophic))
