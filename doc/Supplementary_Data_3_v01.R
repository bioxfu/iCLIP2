# SUPPLEMENTARY DATA 3
# Busch et al, Methods, 2019

# This file provides the R code to reproducibility analyses and assignment of genes and transcript regions as described in Chapters 6.1 and 6.2.

# Running this code requires a bed file with binding sites (e.g. curated PureCLIP output),
# a gtf file with gene/transcript annotations,
# and bw files with crosslink events in the individual replicates (sampleX.strand.bw).

###############################

# The following files need to be specified to run the code:

### Path to .gtf file that holds gene annotations 
annotation.file = < /path/to/annotation.gtf >
  
### Path to folder with bw files
mypath = < /path/to/my/data/dir/ >

# Path to bed file with curated binding sites
bindingsites = < /path/to/bindingsites.bed >
      

### ============================================================================
###
###              Reproducibility & Genomic Annotation Script
###
### ============================================================================

### ============================================================================
### 0) Load packages & files
### ----------------------------------------------------------------------------
### 
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(AnnotationDbi)
library(dplyr)
library(reshape2)
library(UpSetR)
library(GenomicFeatures)

### Load clip tracks per replicate and strand 
myfiles = list.files(mypath, pattern = "*.bw$", full.names = T)
r1.m = import(myfiles[1], "BigWig", as = "Rle")
r1.p = import(myfiles[2], "BigWig", as = "Rle")
r2.m = import(myfiles[3], "BigWig", as = "Rle")
r2.p = import(myfiles[4], "BigWig", as = "Rle")
r3.m = import(myfiles[5], "BigWig", as = "Rle")
r3.p = import(myfiles[6], "BigWig", as = "Rle")
r4.m = import(myfiles[7], "BigWig", as = "Rle")
r4.p = import(myfiles[8], "BigWig", as = "Rle")

### Load binding sites
bindingsites = "./data/curated_peaks_file.bed"
bs = import(bindingsites, format = "BED")

### Object cleaning
bs = keepStandardChromosomes(bs, pruning.mode = "coarse")
mcols(bs) <- NULL



### ============================================================================
### Reproducibility filtering
### ----------------------------------------------------------------------------

### 
### Sum up xlinks per binding site 
###

### Positiv strand
bs.p = bs[strand(bs) == "+"]
bs.p$clp_rep1 = r1.p[bs.p] %>% sum
bs.p$clp_rep2 = r2.p[bs.p] %>% sum
bs.p$clp_rep3 = r3.p[bs.p] %>% sum
bs.p$clp_rep4 = r4.p[bs.p] %>% sum

### Negativ strand
bs.m = bs[strand(bs) == "-"]
bs.m$clp_rep1 = r1.m[bs.m] %>% sum
bs.m$clp_rep2 = r2.m[bs.m] %>% sum
bs.m$clp_rep3 = r3.m[bs.m] %>% sum
bs.m$clp_rep4 = r4.m[bs.m] %>% sum

### Combine
bs = c(bs.p, bs.m)

###
### Make per replicate reproducibility plots 
###

### Scale data
df = mcols(bs) %>% as.matrix
df[df > 40] = 40
df = as.data.frame(df) %>% melt

### Caclulate percentile based threshold
aa = data.frame(cutoff = c(quantile(bs$clp_rep1, probs = seq(0,1, by = 0.1))[2],
                           quantile(bs$clp_rep2, probs = seq(0,1, by = 0.1))[2],
                           quantile(bs$clp_rep3, probs = seq(0,1, by = 0.1))[2],
                           quantile(bs$clp_rep4, probs = seq(0,1, by = 0.1))[2]),
                variable = c("clp_rep1", "clp_rep2", "clp_rep3", "clp_rep4"))
### Add lower boundary 
aa$cutoff[aa$cutoff < 2] = 2

### Make per replicate distribution plot
ggplot(df, aes(x = value, group = variable)) + 
  geom_bar() + 
  facet_grid(~variable) +
  geom_vline(data = aa, aes(xintercept = cutoff), color = "red") +
  geom_text(data = aa, 
            mapping = aes(x = cutoff, y = Inf, label = cutoff, color = "red"),
            hjust = -0.1,
            vjust = 1.5) +
  theme(legend.position="none")

### Make summary plot
names(bs) = 1:length(bs)
NameList = list(rep1 = names(bs[bs$clp_rep1 >= aa$cutoff[1]]),
                rep2 = names(bs[bs$clp_rep2 >= aa$cutoff[2]]),
                rep3 = names(bs[bs$clp_rep3 >= aa$cutoff[3]]),
                rep4 = names(bs[bs$clp_rep3 >= aa$cutoff[4]]))
upset(fromList(NameList), order.by = "freq", nsets = 4)

### Annotate binding sites with reproducibility filter status
status = data.frame(r1 = ifelse(bs$clp_rep1 >= aa$cutoff[1], 1, 0), 
                    r2 = ifelse(bs$clp_rep2 >= aa$cutoff[2], 1, 0),
                    r3 = ifelse(bs$clp_rep3 >= aa$cutoff[3], 1, 0),
                    r4 = ifelse(bs$clp_rep4 >= aa$cutoff[4], 1, 0))
### apply filter to binding sites object
bs$filter = ifelse(rowSums(status) >= 3, T, F)
bs.filter = bs[bs$filter == T]
names(bs.filter) = 1:length(bs.filter)



### ============================================================================
### Genomic Annotation
### ----------------------------------------------------------------------------

### 
### Filter gene annotations 
###

### Import annotation file
anno = import(annotation.file, format = "GTF")

### Filter feature level annotation
anno = anno[anno$level != 3]

### Filter transcript level annotation
anno$transcript_support_level[is.na(anno$transcript_support_level)] = 0
anno$transcript_support_level[anno$transcript_support_level == "NA"] = 10
anno = anno[anno$transcript_support_level == 0
            | anno$transcript_support_level == 1
            | anno$transcript_support_level == 2
            | anno$transcript_support_level == 3 ]

### Create txdb databse from filtered annotations
anno.db = makeTxDbFromGRanges(anno)


### 
### Genomic targets overview
### 

### Select annotated genes
gns = genes(anno.db)
idx = match(gns$gene_id, anno$gene_id)
elementMetadata(gns) = cbind(elementMetadata(gns), elementMetadata(anno)[idx,])

### Gene perspective
df = data.frame(gene_type = subsetByOverlaps(gns, bs.filter)$gene_type) %>% 
  table %>% 
  as.data.frame

ggplot(df, aes(x = reorder(df[,1], -df[,2]), y = df[,2])) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Gene target fractions") + 
  scale_y_log10() + 
  xlab("labels") + 
  ylab("counts [log10]")


### Binding site perspective
bs.tg = subsetByOverlaps(bs.filter,gns)
bs.tg = bs.tg[countOverlaps(bs.tg,targets) == 1]
overlaps = findOverlaps(bs.tg, targets) %>% as.data.frame
df = data.frame(gene_type = targets$gene_type[overlaps$subjectHits]) %>%
  table %>%
  as.data.frame

ggplot(df, aes(x = reorder(df[,1], -df[,2]), y = df[,2])) + 
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Binding site fractions") +
  scale_y_log10() +
  xlab("labels") + 
  ylab("counts [log10]")


### 
### Binding site annotation overlaps
### 

### Select most frequent targets
targets = subsetByOverlaps(gns,bs.filter)
targets.prot = targets[targets$gene_type == "protein_coding"]
bs.prot = subsetByOverlaps(bs.filter,targets.prot)

## How many times does a binding site overlap with two different annotations
df = data.frame(ols = countOverlaps(bs.prot,targets.prot)) %>%
  table %>% as.data.frame
ggplot(df, aes(x = df[,1], y = df[,2])) + geom_bar(stat = "identity") +
  ggtitle("Binding site annoation overlaps") + 
  xlab("annotations") +
  ylab("count [log10]")

### Keep only binding sites that overlap with one distinct target
bs.clean = bs.prot[countOverlaps(bs.prot,targets.prot) == TRUE]


### 
### RBP exact location - Problem visualization
### 

### Count the overlap of each binidng site within each part of the gene
cdseq = cds(anno.db) %>% countOverlaps(bs.clean,.)
intrns = unlist(intronsByTranscript(anno.db)) %>% countOverlaps(bs.clean,.)
utrs3 = unlist(threeUTRsByTranscript(anno.db)) %>% countOverlaps(bs.clean,.)
utrs5 = unlist(fiveUTRsByTranscript(anno.db)) %>% countOverlaps(bs.clean,.)
count.df = data.frame(cds = cdseq, intron = intrns, utr3 = utrs3, utr5 = utrs5)

### Plot the number of different transcript annotaitons at the same position
df = apply(count.df,1,function(x) length(x[x == 0]))
df = data.frame(type = rev(names(table(df))), count = as.vector(table(df)))
ggplot(df, aes(x = df[,1], y = df[,2])) + 
  geom_bar(stat = "identity") + 
  ggtitle("Different transcript region overlaps") + 
  xlab("number of different annotation") +
  ylab("counts")


### 
### RBP exact location - Problem solution
### 

### Setting the hierarchical rule for ties
rule = c("intron", "cds", "utr3", "utr5")

### Applying the majority vote
count.df = count.df[, rule] %>% 
  as.matrix %>% 
  cbind.data.frame(., outside = ifelse(rowSums(count.df) == 0, 1, 0) )
names = colnames(count.df)
reg = apply(count.df, 1, function(x){ names[which.max(x)] })

### Add region annotation to binding sites object
mcols(bs.clean)$region = reg

### Make plot
df = data.frame(region = bs.clean$region)
ggplot(df, aes(x = df[,1])) +
  geom_bar() +
  xlab("region") +
  ylab("count")

### Remove binding sites outside of annotated regions
bs.clean = bs.clean[bs.clean$region != "outside"]



### ============================================================================
### Export filtered & annotated binding sites
### ----------------------------------------------------------------------------

### Write bed file
export(bs.clean, "./bs_clean.bed", format = "BED")
