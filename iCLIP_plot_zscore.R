#! /usr/bin/env Rscript

plot_z_score <- function(z) {
  z_score <- z[,1]
  names(z_score) <- rownames(z)
  top_motif = rev(sort(z_score,decreasing=T)[1:5])
  ht = hist(z_score, col='darkblue', xlab='Z-score',main='',breaks=100)
  axis(1,at=c(-10000,10000))
  axis(2,at=c(-10000,10000))
  txt_x1 = diff(range(ht$breaks))/3
  txt_y1 = seq(max(ht$counts)/3*2,max(ht$counts),len=5)
  text(txt_x1,txt_y1, paste(names(top_motif)[1:5],'   ',sprintf('%-6.2f',top_motif)))
  return(names(top_motif))
}

args <- commandArgs(T)
input <- args[1]
output <- args[2]

pdf(output,width=5,height=5)
x <- read.table(input,row.names=1)
plot_z_score(x)
dev.off()
