wget https://github.com/seqan/flexbar/releases/download/v3.5.0/flexbar-3.5.0-linux.tar.gz
tar zxf flexbar-3.5.0-linux.tar.gz
conda create -n iCLIP2 python=3.5 fastqc fastx_toolkit seqtk star umi_tools bedtools samtools pureclip

DATA=/var/data/raw_data/HJY/HJY_iCLIP_JiangYan190722-X1B_L003/Undetermined_S0_L003_R1_001.fastq.gz
barcodeLength=15
minBaseQuality=10
minReadLength=15
maxReadLength=150
adapter=TGAGATCGGAAGAGCGGTTCAG
genomeMappingIndex=/home/xfu/Gmatic7/genome/human/STAR/
GTF=/home/xfu/Gmatic7/gene/human/GRCh38_v29.gtf
FAI=/home/xfu/Gmatic7/genome/human/GRCh38.fa.fai
#genomeMappingIndex=/home/xfu/Gmatic7/genome/mouse/STAR/
#GTF=/home/xfu/Gmatic7/gene/mouse/GRCm38_vM19.gtf
#FAI=/home/xfu/Gmatic7/genome/mouse/GRCm38.fa.fai

##====================
## Quality control
##====================

###--------------------------
### General quality check
###--------------------------

#### FastQC on the full data set
mkdir fastqc
fastqc --extract --nogroup --outdir fastqc $DATA


###------------------------------------------
### Quality filter on the barcode regions
###------------------------------------------

#### List of read IDs of reads with high quality barcode regions (using FASTX-Toolkit)
mkdir tmp
zcat $DATA | fastx_trimmer -l $barcodeLength | fastq_quality_filter -q $minBaseQuality -p 100 | awk 'FNR%4==1 { print $1 }' | sed 's/@//' > tmp/data.qualFilteredIDs.list

#### Extract reads of given read IDs (using seqtk) and remove problematic characters and whitespaces from read IDs

#seqtk subseq $DATA tmp/data.qualFilteredIDs.list | sed 's/ /#/g; s/\\//#/g' | gzip > filtered.fastq.gz
./seqtk_batch.sh $DATA

###------------------------
### Barcode frequencies 
###------------------------

#### Extract all detected experimental barcodes and their frequencies (x = length of UMI1, y = length of the experimental barcodes)

#zcat filtered.fastq.gz | awk -v umi1_len=5 -v exp_bc_len=6 '{ if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len) }' | sort | uniq -c | sort -k1,1rn > exp_barcodes.detected
mkdir results
./barcode_freq_batch.sh



##=================================================
## Demultiplexing, adapter and barcode trimming
##=================================================

#### check 2nt before adapter sequence if necessary
# zcat x00.fastq.gz|grep 'AGATCGGAAGAGCGGTTCAG'|sed -r 's/AGATCGGAAGAGCGGTTCAG.+//'|awk '{print substr($1, (length($1)-1), 2)}'|sort|uniq -c


#### Demultiplexing, adapter and barcode trimming using Flexbar
mkdir demultiplex
cd demultiplex
#../flexbar-3.5.0-linux/flexbar -n 3 -t x00 -r filtered.fastq.gz --zip-output GZ --barcodes ../barcodes.fasta --barcode-unassigned --barcode-trim-end LTAIL --barcode-error-rate 0 --adapter-seq $adapter --adapter-trim-end RIGHT --adapter-error-rate 0.1 --adapter-min-overlap 1 --min-read-length $minReadLength --umi-tags --length-dist
../demultiplex_batch.sh

cat x*_Ctrl.fastq.gz > Ctrl.fastq.gz
cat x*_KO_rep1.fastq.gz > KO_rep1.fastq.gz
cat x*_KO_rep2.fastq.gz > KO_rep2.fastq.gz
cat x*_KO_rep3.fastq.gz > KO_rep3.fastq.gz
cat x*_WT_rep1.fastq.gz > WT_rep1.fastq.gz
cat x*_WT_rep2.fastq.gz > WT_rep2.fastq.gz
cat x*_WT_rep3.fastq.gz > WT_rep3.fastq.gz

cd ..
cat demultiplex/*.lengthdist|sort -n|uniq|groupBy -g 1 -c 2 -o sum|sed 's/\./Count/' > results/all_reads.lengthdist

#### Plot reads length distribution using FASTX-Toolkit
##### fastq to fasta
#zcat <sampleX.fastq.gz> | fastq_to_fasta -n -r | gzip > <sampleX.fasta.gz>
##### create the plot
#fasta_clipping_histogram.pl <sampleX.fasta.gz> <sampleX.readlength.png>



##====================
## Genomic mapping
##====================

#### Genomic mapping using STAR
mkdir mapping
cd mapping
STAR --runMode alignReads --runThreadN 30 --outFileNamePrefix Ctrl_ --genomeDir $genomeMappingIndex --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --sjdbGTFfile $GTF --sjdbOverhang $maxReadLength-1 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn ../demultiplex/Ctrl.fastq.gz
STAR --runMode alignReads --runThreadN 30 --outFileNamePrefix KO_rep1_ --genomeDir $genomeMappingIndex --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --sjdbGTFfile $GTF --sjdbOverhang $maxReadLength-1 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn ../demultiplex/KO_rep1.fastq.gz
STAR --runMode alignReads --runThreadN 30 --outFileNamePrefix KO_rep2_ --genomeDir $genomeMappingIndex --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --sjdbGTFfile $GTF --sjdbOverhang $maxReadLength-1 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn ../demultiplex/KO_rep2.fastq.gz
STAR --runMode alignReads --runThreadN 30 --outFileNamePrefix KO_rep3_ --genomeDir $genomeMappingIndex --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --sjdbGTFfile $GTF --sjdbOverhang $maxReadLength-1 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn ../demultiplex/KO_rep3.fastq.gz
STAR --runMode alignReads --runThreadN 30 --outFileNamePrefix WT_rep1_ --genomeDir $genomeMappingIndex --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --sjdbGTFfile $GTF --sjdbOverhang $maxReadLength-1 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn ../demultiplex/WT_rep1.fastq.gz
STAR --runMode alignReads --runThreadN 30 --outFileNamePrefix WT_rep2_ --genomeDir $genomeMappingIndex --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --sjdbGTFfile $GTF --sjdbOverhang $maxReadLength-1 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn ../demultiplex/WT_rep2.fastq.gz
STAR --runMode alignReads --runThreadN 30 --outFileNamePrefix WT_rep3_ --genomeDir $genomeMappingIndex --outFilterMismatchNoverReadLmax 0.04 --outFilterMismatchNmax 999 --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --sjdbGTFfile $GTF --sjdbOverhang $maxReadLength-1 --outReadsUnmapped Fastx --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --readFilesIn ../demultiplex/WT_rep3.fastq.gz
find *.bam|parallel --gnu "samtools index {}"
cd ..


##======================================
## Duplicate removal (deduplication)
##======================================

#### Duplicate removal (deduplication) using UMI-tools
mkdir dedup
find mapping/*.bam|sed 's/mapping\///;s/\.bam//'|xargs -I {} umi_tools dedup -I mapping/{}.bam -L dedup/{}.duprm.log -S dedup/{}.duprm.bam --extract-umi-method read_id --method unique



##==========================================
## Extraction of crosslinked nucleotides
##==========================================

#### Create file of chromosome length using SAMtools (produces a file genome.fasta.fai)
#samtools faidx <genome.fasta>


#### Convert all read locations to intervals in bed file format using BEDTools
source activate umi_tools
mkdir bed
find dedup/*.bam|sed 's/dedup\///;s/\.bam//'|parallel --gnu "bedtools bamtobed -i dedup/{}.bam > bed/{}.bed"
 

#### Shift intervals depending on the strand by 1 bp upstream using BEDTools
find bed/*.bed|sed 's/bed\///;s/\.bed//'|parallel --gnu "bedtools shift -m 1 -p -1 -i bed/{}.bed -g $FAI > bed/{}.shifted.bed"
 

#### Extract the 5' end of the shifted intervals and pile up into coverage track in bedgraph file format (separately for each strand) using BEDTools (in case of RPM-normalised coverage tracks, use additional parameter -scale with 1,000,000/#mappedReads)

find bed/*.shifted.bed|sed 's/bed\///;s/\.bed//'|parallel --gnu "bedtools genomecov -bg -strand + -5 -i bed/{}.bed -g $FAI > bed/{}.plus.bedgraph"
find bed/*.shifted.bed|sed 's/bed\///;s/\.bed//'|parallel --gnu "bedtools genomecov -bg -strand - -5 -i bed/{}.bed -g $FAI > bed/{}.minus.bedgraph"


#### make tdf file
mkdir track
find bed/*.bedgraph|sed 's/bed\///;s/\.bedgraph//'|parallel --gnu "igvtools toTDF bed/{}.bedgraph bed/{}.tdf /home/xfu/igv/genomes/hg38.genome"



#### Optional convertion of bedgraph files to bw file format files using bedGraphToBigWig of the kentUtils suite
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
chmod +x bedGraphToBigWig 

export LC_COLLATE=C
find bed/*.bedgraph|sed 's/\.bedgraph//'|parallel --gnu "sort -k1,1 -k2,2n {}.bedgraph > {}.sorted.bedgraph"
find bed/*.sorted.bedgraph|sed 's/\.bedgraph//'|parallel --gnu "./bedGraphToBigWig {}.bedgraph $FAI {}.bw"

#### Depending on the system and the version of bedGraphToBigWig, it might be necessary to sort the bedgraph files before converting them to bw files:
#export LC_COLLATE=C
#sort -k1,1 -k2,2n <sampleX.strand.bedgraph> > <sampleX.strand.sorted.bedgraph>



##========================================================
## Diagnostic plots and measures of library complexity
##========================================================

###------------------------------
### Duplicate removal summary
###------------------------------
#sampleX=Ctrl
##### Number of crosslink events, i.e. reads after duplicate removal:
#cat bed/${sampleX}_Aligned.sortedByCoord.out.duprm.shifted.*.bedgraph | awk 'BEGIN{ totalcount=0 }{ totalcount += (($3 - $2) * $4) }END{ print "crosslink_events\t"totalcount }' > results/${sampleX}_Xlink_stat.tsv
##### Number of crosslinked nucleotides, i.e. positions harbouring crosslinked nucleotides (if both strands are covered, count as 2):
#cat bed/${sampleX}_Aligned.sortedByCoord.out.duprm.shifted.*.bedgraph | awk 'BEGIN{ totalpos=0 }{ totalpos += ($3 - $2) }END{ print "crosslinked_nucleotides\t"totalpos }' >> results/${sampleX}_Xlink_stat.tsv
##### Number of stacked crosslink events, i.e. crosslink events on positions with >1 crosslink events:
#cat bed/${sampleX}_Aligned.sortedByCoord.out.duprm.shifted.*.bedgraph | awk 'BEGIN{ totalstackedcount=0 }{ if($4 > 1) totalstackedcount += (( $3 - $2) * $4) }END{ print "stacked_crosslink_events\t"totalstackedcount }' >> results/${sampleX}_Xlink_stat.tsv
##### Number of nucleotides with stacked crosslink events, i.e. positions with >1 crosslink events: 
#cat bed/${sampleX}_Aligned.sortedByCoord.out.duprm.shifted.*.bedgraph | awk 'BEGIN{ totalstackedpos=0 }{ if($4 > 1) totalstackedpos += ($3 - $2) }END{ print "nucleotides_with_stacked_crosslink_events\t"totalstackedpos }' >> results/${sampleX}_Xlink_stat.tsv

./Xlink_stat.sh Ctrl
./Xlink_stat.sh KO_rep1
./Xlink_stat.sh KO_rep2
./Xlink_stat.sh KO_rep3
./Xlink_stat.sh WT_rep1
./Xlink_stat.sh WT_rep2
./Xlink_stat.sh WT_rep3

./deduplication_stat.sh Ctrl
./deduplication_stat.sh KO_rep1
./deduplication_stat.sh KO_rep2
./deduplication_stat.sh KO_rep3
./deduplication_stat.sh WT_rep1
./deduplication_stat.sh WT_rep2
./deduplication_stat.sh WT_rep3

# make_stat_table
paste results/dedup_stat*.tsv|cut -f1,2,4,6,8,10,12,14,16 > stat_table.tsv
paste results/Xlink_stat*.tsv|cut -f1,2,4,6,8,10,12,14,16|grep -v 'Ctrl' >> stat_table.tsv

###----------------------------------------
### Reads with insertions and deletions
###----------------------------------------
#### Convert bam to sam file format (using SAMtools)
#samtools view <sampleX.duprm.bam> -o <sampleX.duprm.sam>
#### Number of reads mapped with deletions:
#cut -f6 <sampleX.duprm.sam> | grep D | wc -l
#### Number of reads mapped with insertions:
#cut -f6 <sampleX.duprm.sam> | grep I | wc -l


###---------------------
### iCLIPro analysis
###---------------------
#### iCLIPro analysis and plots
#iCLIPro -r 50 -b 300 -f 50 -g "L15:15,L16:16,L17:17,L18:18,L19:19,L20:20,L21:21,L22:22,L23:23,L24:24,L25:25,L26:26,L27:27,L28:28,L29:29,L30:30,L31:31,L32:32,L33:33,L34:34,L35:35,L36:36,L37:37,L38:38,L39:39,L40:40,R:41" -p "L15-R,L16-R,L17-R,L18-R,L19-R,L20-R,L21-R,L22-R,L23-R,L24-R,L25-R,L26-R,L27-R,L28-R,L29-R,L30-R,L31-R,L32-R,L33-R,L34-R,L35-R,L36-R,L37-R,L38-R,L39-R,L40-R" -o <outdir> <sampleX.duprm.bam>



##===============================
## Peak calling with PureCLIP
##===============================

#### Merge BAM files (sampleX.duprm.bam)
mkdir peak
cd peak

samtools merge -f WT_merged.bam ../dedup/WT_rep1_Aligned.sortedByCoord.out.duprm.bam ../dedup/WT_rep2_Aligned.sortedByCoord.out.duprm.bam ../dedup/WT_rep3_Aligned.sortedByCoord.out.duprm.bam
samtools index WT_merged.bam

samtools merge -f KO_merged.bam ../dedup/KO_rep1_Aligned.sortedByCoord.out.duprm.bam ../dedup/KO_rep2_Aligned.sortedByCoord.out.duprm.bam ../dedup/KO_rep3_Aligned.sortedByCoord.out.duprm.bam
samtools index KO_merged.bam

#### Run PureCLIP on fat node (1T RAM)
GENOME=/home/xfu/Gmatic7/genome/human/GRCh38.fa
pureclip -i KO_merged.bam -bai KO_merged.bam.bai -g $GENOME -ld -nt 10 -o KO_PureCLIP.crosslink_sites.bed -or KO_PureCLIP.crosslink_regions.bed
pureclip -i WT_merged.bam -bai WT_merged.bam.bai -g $GENOME -ld -nt 10 -o WT_PureCLIP.crosslink_sites.bed -or WT_PureCLIP.crosslink_regions.bed


#### Remove 7th column of PureCLIP output file
cat KO_PureCLIP.crosslink_sites.bed | cut -f 1,2,3,4,5,6|grep '^chr' > KO_PureCLIP.crosslink_sites_short.bed
cat WT_PureCLIP.crosslink_sites.bed | cut -f 1,2,3,4,5,6|grep '^chr' > WT_PureCLIP.crosslink_sites_short.bed


bedtools bamtobed -i peak/KO_merged.bam > peak/KO_merged.bed
bedtools bamtobed -i peak/WT_merged.bam > peak/WT_merged.bed

bedtools shift -m 1 -p -1 -i peak/KO_merged.bed -g $FAI > peak/KO_merged.shifted.bed
bedtools shift -m 1 -p -1 -i peak/WT_merged.bed -g $FAI > peak/WT_merged.shifted.bed

bedtools genomecov -bg -strand + -5 -i peak/KO_merged.shifted.bed -g $FAI|grep '^chr' |sortBed > peak/KO_merged.shifted.plus.bedgraph &
bedtools genomecov -bg -strand - -5 -i peak/KO_merged.shifted.bed -g $FAI|grep '^chr' |sortBed > peak/KO_merged.shifted.minus.bedgraph &
bedtools genomecov -bg -strand + -5 -i peak/WT_merged.shifted.bed -g $FAI|grep '^chr' |sortBed > peak/WT_merged.shifted.plus.bedgraph &
bedtools genomecov -bg -strand - -5 -i peak/WT_merged.shifted.bed -g $FAI|grep '^chr' |sortBed > peak/WT_merged.shifted.minus.bedgraph &

./bedGraphToBigWig peak/KO_merged.shifted.plus.bedgraph $FAI peak/KO_merged.shifted.plus.bw &
./bedGraphToBigWig peak/KO_merged.shifted.minus.bedgraph $FAI peak/KO_merged.shifted.minus.bw &
./bedGraphToBigWig peak/WT_merged.shifted.plus.bedgraph $FAI peak/WT_merged.shifted.plus.bw &
./bedGraphToBigWig peak/WT_merged.shifted.minus.bedgraph $FAI peak/WT_merged.shifted.minus.bw &

~/R/3.6.1/bin/Rscript iCLIP2_postprocess_PureCLIP_peak.R WT
~/R/3.6.1/bin/Rscript iCLIP2_postprocess_PureCLIP_peak.R KO

~/R/3.6.1/bin/Rscript iCLIP2_binding_sites_clean.R WT
~/R/3.6.1/bin/Rscript iCLIP2_binding_sites_clean.R KO


~/R/3.6.1/bin/Rscript iCLIP2_binding_sites_clean.R WT

./motif_zscore.sh

~/R/3.6.1/bin/Rscript motif_zscore.R

