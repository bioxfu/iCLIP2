library(RColorBrewer)
library(scales)

col_set <- brewer.pal(n = 8, 'Set1')

smooth_occu <- function(x, step=3) {
  y <- x
  x <- c(x[1:step], x)
  for (i in 1:length(y)) {
    y[i] <- mean(x[i:(i+step)])
  }
  return(y)
}

fix_short <- function(x) {
  full <- data.frame(V1 = 0:149, V2 = 0)
  full[full$V1 %in% x$V1, 2] <- x$V2
  return(full)
}

ko_si_5ss <- fix_short(read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_KO_dist2occu'))
wt_si_5ss <- fix_short(read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_WT_dist2occu'))
ko_en_5ss <- fix_short(read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_KO_dist2occu'))
wt_en_5ss <- fix_short(read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_WT_dist2occu'))
ko_si_3ss <- fix_short(read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_KO_dist2occu'))
wt_si_3ss <- fix_short(read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_WT_dist2occu'))
ko_en_3ss <- fix_short(read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_KO_dist2occu'))
wt_en_3ss <- fix_short(read.table('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_WT_dist2occu'))

fc_cutoff <- 2
p_cutoff <- 0.01

#pdf('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_RNAmap_v1.pdf', hei=5) # p.adjust(p)
pdf('occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_RNAmap_v2.pdf', hei=5) # p.adjust(p, method = 'bonferroni') 
par(mfrow=c(2, 2))
#par(mar=c(5,4,4,1))
par(mar=c(2,4,1,1))
# enhanced SE events
plot(wt_en_3ss$V1, smooth_occu(wt_en_3ss$V2, step=10)+0.05, ylim=c(0,0.6), yaxt='n', type='l', bty='l', las=2, lwd=1.5, col=col_set[1], ylab='Normalized crosslink events', xlab='', xaxt='n')
lines(ko_en_3ss$V1, smooth_occu(ko_en_3ss$V2, step=10)+0.05, col=col_set[2], lwd=1.5)
axis(1, at=c(0,50,100,150), lab=c(-100, -50, '3SS', 50))
axis(2, at=c(0, 0.1, 0.2, 0.3, 0.4, 0.5)+0.05, lab=c(0, 0.1, 0.2, 0.3, 0.4, 0.5), las=2)
legend('topleft', c('WT', 'KO'), col = col_set[1:2], lwd = 3, bty='n')
abline(h=0)
rect(100, -0.03, 200, 0.03, col='gray50', border = 'gray50')
x <- smooth_occu(wt_en_3ss$V2, step=10)
y <- smooth_occu(ko_en_3ss$V2, step=10)
x <- c(x, x[-(1:120)])
y <- c(y, y[-(1:120)])
p <- c()
fc <- c()
for (i in 1:150) {
  p[i] <- t.test(x[i:(i+30)], y[i:(i+30)])$p.value
  fc[i] <- log2(mean(x[i:(i+30)])/mean(y[i:(i+30)]))
}
fc[is.na(fc)] <- 0
start <- which(p.adjust(p, method = 'bonferroni') < p_cutoff & fc > log2(fc_cutoff))
end <- which(p.adjust(p, method = 'bonferroni') < p_cutoff & fc > log2(fc_cutoff))+30
rect(min(start), 0, max(end), 10, col = alpha(col=col_set[1], 0.1), border = NA)

#par(mar=c(5,1,4,4))
par(mar=c(2,1,1,4))
plot(wt_en_5ss$V1, smooth_occu(wt_en_5ss$V2, step=10)+0.05, ylim=c(0,0.6), yaxt='n', type='l', bty='n', las=2, lwd=1.5, col=col_set[1], ylab='', xlab='', main='', xaxt='n')
lines(ko_en_5ss$V1, smooth_occu(ko_en_5ss$V2, step=10)+0.05, col=col_set[2], lwd=1.5)
axis(1, at=c(-100,0,50,100,150,200), lab=c('', -50, '5SS', 50, 100, ''))
abline(h=0)
rect(-50, -0.03, 50, 0.03, col='gray50', border = 'gray50')
x <- smooth_occu(wt_en_5ss$V2, step=10)
y <- smooth_occu(ko_en_5ss$V2, step=10)
x <- c(x, x[-(1:120)])
y <- c(y, y[-(1:120)])
p <- c()
fc <- c()
for (i in 1:150) {
  p[i] <- t.test(x[i:(i+30)], y[i:(i+30)])$p.value
  fc[i] <- log2(mean(x[i:(i+30)])/mean(y[i:(i+30)]))
}
start <- which(p.adjust(p, method = 'bonferroni') < p_cutoff & fc > log2(fc_cutoff))
end <- which(p.adjust(p, method = 'bonferroni') < p_cutoff & fc > log2(fc_cutoff))+30
rect(min(start), 0, max(end), 10, col = alpha(col=col_set[1], 0.1), border = NA)

#par(mar=c(5,4,4,1))
par(mar=c(2,4,1,1))
# silenced SE events
plot(wt_si_3ss$V1, smooth_occu(wt_si_3ss$V2, step=10)+0.05, ylim=c(0,0.6), yaxt='n', type='l', bty='l', las=2, lwd=1.5, col=col_set[3], ylab='Normalized crosslink events', xlab='', xaxt='n')
lines(ko_si_3ss$V1, smooth_occu(ko_si_3ss$V2, step=10)+0.05, col=col_set[4], lwd=1.5)
axis(1, at=c(0,50,100,150), lab=c(-100, -50, '3SS', 50))
axis(2, at=c(0, 0.1, 0.2, 0.3, 0.4, 0.5)+0.05, lab=c(0, 0.1, 0.2, 0.3, 0.4, 0.5), las=2)
legend('topleft', c('WT', 'KO'), col = col_set[3:4], lwd = 3, bty='n')
abline(h=0)
rect(100, -0.03, 200, 0.03, col='gray50', border = 'gray50')
x <- smooth_occu(wt_si_3ss$V2, step=10)
y <- smooth_occu(ko_si_3ss$V2, step=10)
x <- c(x, x[-(1:120)])
y <- c(y, y[-(1:120)])
p <- c()
fc <- c()
for (i in 1:150) {
  p[i] <- t.test(x[i:(i+30)], y[i:(i+30)])$p.value
  fc[i] <- log2(mean(x[i:(i+30)])/mean(y[i:(i+30)]))
}
start <- which(p.adjust(p, method = 'bonferroni') < p_cutoff & fc > log2(fc_cutoff))
end <- which(p.adjust(p, method = 'bonferroni') < p_cutoff & fc > log2(fc_cutoff))+30
rect(min(start), 0, max(end), 10, col = alpha(col=col_set[1], 0.1), border = NA)

#par(mar=c(5,1,4,4))
par(mar=c(2,1,1,4))
plot(wt_si_5ss$V1, smooth_occu(wt_si_5ss$V2, step=10)+0.05, ylim=c(0,0.6), type='l', bty='n', yaxt='n', las=2, lwd=1.5, col=col_set[3], ylab='', xlab='', main='', xaxt='n')
lines(ko_si_5ss$V1, smooth_occu(ko_si_5ss$V2, step=10)+0.05, col=col_set[4], lwd=1.5)
axis(1, at=c(-100,0,50,100,150,200), lab=c('', -50, '5SS', 50, 100, ''))
abline(h=0)
rect(-50, -0.03, 50, 0.03, col='gray50', border = 'gray50')
x <- smooth_occu(wt_si_5ss$V2, step=10)
y <- smooth_occu(ko_si_5ss$V2, step=10)
x <- c(x, x[-(1:120)])
y <- c(y, y[-(1:120)])
p <- c()
fc <- c()
for (i in 1:150) {
  p[i] <- t.test(x[i:(i+30)], y[i:(i+30)])$p.value
  fc[i] <- log2(mean(x[i:(i+30)])/mean(y[i:(i+30)]))
}
start <- which(p.adjust(p, method = 'bonferroni') < p_cutoff & fc > log2(fc_cutoff))
end <- which(p.adjust(p, method = 'bonferroni') < p_cutoff & fc > log2(fc_cutoff))+30
rect(min(start), 0, max(end), 10, col = alpha(col=col_set[1], 0.1), border = NA)

dev.off()


