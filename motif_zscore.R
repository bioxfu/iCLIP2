library(RColorBrewer)
library(scales)
cols <- brewer.pal(n = 8, 'Set1')

cal_zscore <- function(dumps, N, K) {
  kmer <- read.table(paste0('zscore/crosslink.window.', K, '.list'), stringsAsFactors = F)[,1]
  x <- matrix(0, nrow=length(kmer), ncol=101)
  rownames(x) <- kmer
  colnames(x) <- 0:100
  
  input <- read.table(dumps, stringsAsFactors=F)
  x[input[, 1], 1] <- input[, 2]
  
  files <- list.files(paste0('zscore/', K, '/random', N), '*.dump', full.names = T)
  for (i in 1:length(files)) {
    input <- read.table(files[i], stringsAsFactors=F)
    x[input[, 1], i+1] <- input[, 2]
  }
  
  z_score <- (x[, 1]- apply(x[, -1], 1, mean)) / apply(x[, -1], 1, sd)
  names(z_score) <- gsub('T', 'U', names(z_score))
  #z_score <- data.frame(score=sort(z_score, decreasing=T), rank=1:length(z_score))
  return(z_score)
}

K <- '5mer'

WT_zscore <- NULL
KO_zscore <- NULL

for (N in 1:10) {
  zscore1 <- cal_zscore(paste0('zscore/WT_binding_region_specific_intron_7K_', K, '.dumps'), N, K)
  zscore2 <- cal_zscore(paste0('zscore/KO_binding_region_specific_intron_7K_', K, '.dumps'), N, K)
  WT_zscore <- cbind(WT_zscore, zscore1)
  KO_zscore <- cbind(KO_zscore, zscore2)
}

all <- data.frame(WT_zscore = rowMeans(WT_zscore, 1),
                  KO_zscore = rowMeans(KO_zscore, 1))
all <- all[order(-all$WT_zscore),]
all$WT_zscore_rank <- 1:nrow(all)
all <- all[order(-all$KO_zscore),]
all$KO_zscore_rank <- 1:nrow(all)

motifs <- read.table('motif_hnRNP', stringsAsFactors = F)
colnames(motifs) <- c('motif', 'class')
motifs <- merge(motifs, all, by.x=1, by.y=0, all.x = T)

pdf(paste0('zscore/', K, '_zscore_motif_hnRNP.pdf'), wid=9, hei=9)
layout(matrix(c(1,2,3,3), nrow = 2, byrow = T))
#par(mfrow=c(1, 3))
motifs_color <- motifs$class
motifs_color[motifs_color=='Fox'] <- cols[1]   
motifs_color[motifs_color=='green'] <- cols[3]
motifs_color[motifs_color=='blue'] <- cols[2]
motifs_color[motifs_color=='orange'] <- cols[5]

hist(all$WT_zscore, ylim=c(-40, 550), breaks = 30, ylab='Number of pentamers', xlab='pentamer z-score', main='WT', col='gray90')
axis(1, at=c(-1000, 1000))
points(all[motifs$motif[motifs$class=='Fox'],'WT_zscore'], rep(-10, sum(motifs$class=='Fox')), col=motifs_color[motifs$class=='Fox'], pch=16, cex=0.9)
points(all[motifs$motif[motifs$class=='orange'],'WT_zscore'], rep(-10, sum(motifs$class=='orange')), col=motifs_color[motifs$class=='orange'], pch=16, cex=0.9)
points(all[motifs$motif[motifs$class=='blue'],'WT_zscore'], rep(-25, sum(motifs$class=='blue')), col=motifs_color[motifs$class=='blue'], pch=16, cex=0.9)
points(all[motifs$motif[motifs$class=='green'],'WT_zscore'], rep(-40, sum(motifs$class=='green')), col=motifs_color[motifs$class=='green'], pch=16, cex=0.9)
points(all[motifs$motif[motifs$class=='gray'],'WT_zscore'], rep(-55, sum(motifs$class=='gray')), col=motifs_color[motifs$class=='gray'], pch=16, cex=0.9)

hist(all$KO_zscore, ylim=c(-40, 550), breaks = 30, ylab='Number of pentamers', xlab='pentamer z-score', main='KO', col='gray90')
axis(1, at=c(-1000, 1000))
points(all[motifs$motif[motifs$class=='Fox'],'KO_zscore'], rep(-10, sum(motifs$class=='Fox')), col=motifs_color[motifs$class=='Fox'], pch=16, cex=0.9)
points(all[motifs$motif[motifs$class=='orange'],'KO_zscore'], rep(-10, sum(motifs$class=='orange')), col=motifs_color[motifs$class=='orange'], pch=16, cex=0.9)
points(all[motifs$motif[motifs$class=='blue'],'KO_zscore'], rep(-25, sum(motifs$class=='blue')), col=motifs_color[motifs$class=='blue'], pch=16, cex=0.9)
points(all[motifs$motif[motifs$class=='green'],'KO_zscore'], rep(-40, sum(motifs$class=='green')), col=motifs_color[motifs$class=='green'], pch=16, cex=0.9)
points(all[motifs$motif[motifs$class=='gray'],'KO_zscore'], rep(1-55, sum(motifs$class=='gray')), col=motifs_color[motifs$class=='gray'], pch=16, cex=0.9)
legend('topright', c('RBFOX2 motifs', 'RBFOX2 motifs', 'hnRNP M motifs', 'hnRNP C motifs', 'hnRNP H motifs'), pch=16,col = c(cols[1], cols[5], cols[2], cols[3], 'gray'), bty='n', border = NA)

# plot(motifs$WT_zscore_rank, motifs$KO_zscore_rank, xlim = c(0,1024), ylim=c(0,1024), col=motifs_color, pch=16, cex=2, lwd=2, xlab = 'Z-score rank in WT', ylab = 'Z-score rank in KO')
# abline(a=0, b=1, lty=2)

write.table(motifs[order(-motifs$WT_zscore), ], paste0('zscore/', K, '_zscore_motif_hnRNP.tsv'), row.names = F, sep = '\t')

motifs <- motifs[motifs$class != 'Fox',]
motifs$rank_diff <- motifs$WT_zscore_rank - motifs$KO_zscore_rank
x <- sort(motifs$rank_diff)
n <- motifs$motif[order(motifs$rank_diff)]
y <- as.character(motifs$class[order(motifs$rank_diff)])
y[y=='green'] <- cols[3]
y[y=='blue'] <- cols[2]

barplot(x, names.arg = n, las=2, col=y, border = NA, ylab='Z-score rank difference (WT - KO)', space = 0.5, ylim=c(-250, 400), yaxt='n')
axis(2, at=seq(-250, 400, 50), lab=seq(-250, 400, 50))
legend('bottomright', c('hnRNP M motifs', 'hnRNP C motifs', 'hnRNP H motifs'), fill = c(cols[2], cols[3], 'gray'), bty='n', border = NA)
dev.off()

