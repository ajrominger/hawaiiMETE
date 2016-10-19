library(meteR)

setwd('~/Dropbox/Research/hawaiiMETE')

load('mete.byTS.RData')

with(byTrophBySite[byTrophBySite$trophic == 'H', ], 
     plot(factor(as.character(Site), levels = c('VO', 'LA', 'KH', 'MO', 'KA')), 
          S0))


with(byTrophBySite[byTrophBySite$trophic == 'H', ], plot(S0, sad.z, xlim = range(byTrophBySite$S0)))

plot(byTrophBySite$S0, byTrophBySite$sad.z, col = factor(byTrophBySite$trophic))
