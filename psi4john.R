library(meteR)
library(socorro)
load('mete.byTS.RData')

pdf('fig_psi4john.pdf', width = 5, height = 5)
par(mar = c(3, 3, 2, 0) + 0.1, mgp = c(2, 0.75, 0))
i <- 6
plot(mete.byTS[[i]]$ipd, ptype = 'rad')
mtext('Kohala predatory arthropods', line = 1)
dev.off()
