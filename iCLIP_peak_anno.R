suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(ChIPseeker))

argv <- commandArgs(T)
input_bed <- argv[1]
sqlite <- argv[2]
output_pie <- argv[3]

# input_bed <- 'peak/hs_PureCLIP.crosslink_sites_short.bed'
# sqlite <- '/home/xfu/Gmatic7/gene/human/txdb/GRCh38_v29_txdb.sqlite'
# output_pie <- 'figure/test_pie.pdf'

peak <- readPeakFile(input_bed)
txdb <- loadDb(sqlite)

peakAnno <- annotatePeak(peak, TxDb=txdb, tssRegion = c(-2000,0))
anno <- peakAnno@annoStat
anno2 <- data.frame(Feature = c('Intergenic', 'Exon', 'Intron', "5' UTR", "3' UTR"),
                    Frequency = c(sum(anno[c(grep('Promoter', anno$Feature), grep('Downstream', anno$Feature), grep('Intergenic', anno$Feature)), 2]),
                                  sum(anno[grep('Exon', anno$Feature), 2]),
                                  sum(anno[grep('Intron', anno$Feature), 2]),
                                  anno[grep("5' UTR", anno$Feature), 2],
                                  anno[grep("3' UTR", anno$Feature), 2]))

peakAnno@annoStat <- anno2 

pdf(output_pie, width=7, height=5)
plotAnnoPie(peakAnno)
dev.off()


