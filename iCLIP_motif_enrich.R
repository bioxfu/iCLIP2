#! /usr/bin/env Rscript
library(RWebLogo)

args <- commandArgs(trailingOnly = T)
input_path <- args[1]
input_kmer <- args[2]
input_file0 <- args[3]
input_pattern <- args[4]
output <- args[5]

kmer <- read.table(paste0(input_path,'/',input_kmer))[,1]
x <- matrix(0,nrow=length(kmer),ncol=101)
rownames(x) <- kmer
colnames(x) <- 0:100
input <- read.table(paste0(input_path,'/',input_file0), stringsAsFactors=F)
x[input[,1],1] <- input[,2]

files = paste0(input_path,'/',list.files(input_path,input_pattern))
for (i in 1:length(files)) {
  input = read.table(files[i], stringsAsFactors=F)
  x[input[,1],i+1] <- input[,2]
}
z_score = (x[,1]-apply(x[,-1],1,mean))/apply(x[,-1],1,sd)
names(z_score) = gsub('T','U',names(z_score))
z_score = sort(z_score,decreasing=T)
write.table(x, paste0(output,'.kmer.tsv'), quote=F, sep='\t', col.names=NA)
write.table(z_score, paste0(output,'.zscore.tsv'), quote=F, sep='\t', col.names=F)

plot_z_score = function(z) {
  top_motif = rev(sort(z,decreasing=T)[1:5])
  ht = hist(z, col='gray', xlab='Z-score',main='')
  txt_x1 = diff(range(ht$breaks))/3
  txt_y1 = seq(max(ht$counts)/3*2,max(ht$counts),len=5)
  text(txt_x1,txt_y1, paste(names(top_motif),sprintf('%-6.2f',top_motif)))
  return(names(top_motif))
}

#pdf(paste0(output,'.zscore.pdf'),heigh=5,width=5)
#tpm = plot_z_score(z_score)
#dev.off()
#weblogo(tpm,errorbars=FALSE,annotate=c(1:nchar(names(z_score)[1])),file.out=paste0(output,'.logo.pdf'),color.scheme='classic',open=F)



