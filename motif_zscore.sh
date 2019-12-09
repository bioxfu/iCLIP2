## In order to compare the Z-scores of motifs bwtween two samples, 
## we subsample the same number of crosslink sites, and generate 
## one permutation backgroud for two samples.

GENOME=$HOME/Gmatic7/genome/human/GRCh38.fa
GENOMESIZE=$HOME/Gmatic7/genome/human/GRCh38.fa.fai
INTRON=$HOME/Gmatic7/gene/human/GRCh38_v29_intron.bed
CPU=8
K=5
GTF=/home/xfu/Gmatic7/gene/human/GRCh38_v29.gtf

mkdir zscore
# filter binding sites
./filter_sites_by_RBP_motif.sh $GENOME 'GCAUG|UGCAU' results/WT_binding_site_clean.bed zscore/WT_binding_site_with_RBP_motif.bed
./filter_sites_by_RBP_motif.sh $GENOME 'GCAUG|UGCAU' results/KO_binding_site_clean.bed zscore/KO_binding_site_with_RBP_motif.bed

# define binding regions
cat zscore/WT_binding_site_with_RBP_motif.bed|awk '{print $1"\t"$2-40"\t"$3+40"\t"$4"\t"$5"\t"$6}' > zscore/WT_binding_region.bed
cat zscore/KO_binding_site_with_RBP_motif.bed|awk '{print $1"\t"$2-40"\t"$3+40"\t"$4"\t"$5"\t"$6}' > zscore/KO_binding_region.bed

# sample specific binding sites
bedtools intersect -a zscore/WT_binding_region.bed -b zscore/KO_binding_region.bed -v > zscore/WT_binding_region_specific.bed
bedtools intersect -a zscore/KO_binding_region.bed -b zscore/WT_binding_region.bed -v > zscore/KO_binding_region_specific.bed

# intronic sites
bedtools intersect -a zscore/WT_binding_region_specific.bed -b $INTRON -wa -u -f 1.0 > zscore/WT_binding_region_specific_intron.bed
bedtools intersect -a zscore/KO_binding_region_specific.bed -b $INTRON -wa -u -f 1.0 > zscore/KO_binding_region_specific_intron.bed

# keep same sites
shuf zscore/WT_binding_region_specific_intron.bed|head -7000 > zscore/WT_binding_region_specific_intron_7K.bed
shuf zscore/KO_binding_region_specific_intron.bed|head -7000 > zscore/KO_binding_region_specific_intron_7K.bed

bedtools getfasta -fi $GENOME -bed zscore/WT_binding_region_specific_intron_7K.bed -name -s -fo zscore/WT_binding_region_specific_intron_7K.fa
bedtools getfasta -fi $GENOME -bed zscore/KO_binding_region_specific_intron_7K.bed -name -s -fo zscore/KO_binding_region_specific_intron_7K.fa

jellyfish count -m $K -s 100M -t $CPU zscore/WT_binding_region_specific_intron_7K.fa -o zscore/WT_binding_region_specific_intron_7K_${K}mer.jf
jellyfish count -m $K -s 100M -t $CPU zscore/KO_binding_region_specific_intron_7K.fa -o zscore/KO_binding_region_specific_intron_7K_${K}mer.jf

jellyfish dump zscore/WT_binding_region_specific_intron_7K_${K}mer.jf_0 -c -t|sort > zscore/WT_binding_region_specific_intron_7K_${K}mer.dumps
jellyfish dump zscore/KO_binding_region_specific_intron_7K_${K}mer.jf_0 -c -t|sort > zscore/KO_binding_region_specific_intron_7K_${K}mer.dumps

cat zscore/*${K}mer.dumps|cut -f1|sort|uniq > zscore/crosslink.window.${K}mer.list
rm zscore/*.jf_0

# intron containing iclip regions
bedtools intersect -a $INTRON -b zscore/WT_binding_region.bed -u >  tmp.bed
bedtools intersect -a $INTRON -b zscore/KO_binding_region.bed -u >> tmp.bed
sortBed -i tmp.bed|awk '{print $1"\t"$2"\t"$3"\t.\t.\t"$6}'|uniq|mergeBed -s > intron_with_iclip_regions.bed
rm tmp.bed

#
cat zscore/WT_binding_region.bed zscore/KO_binding_region.bed|sortBed|mergeBed > all_binding_regions.bed

## generate 10 random background
mkdir zscore/${K}mer
for N in {1..10}
do
mkdir zscore/${K}mer/random
for i in {1..100}
do
	echo $i
    bedtools shuffle -i zscore/WT_binding_region_specific_intron_7K.bed -g $GENOMESIZE -incl intron_with_iclip_regions.bed -excl all_binding_regions.bed > zscore/${K}mer/random/random.$i.bed
    bedtools getfasta -fi $GENOME -bed zscore/${K}mer/random/random.$i.bed -name -s -fo zscore/${K}mer/random/random.$i.fa
    jellyfish count -m $K -s 100M zscore/${K}mer/random/random.$i.fa -o zscore/${K}mer/random/random.${K}mer.$i.jf
    jellyfish dump zscore/${K}mer/random/random.${K}mer.$i.jf_0 -c -t |sort > zscore/${K}mer/random/random.${K}mer.$i.dumps
	rm zscore/${K}mer/random/random.$i.bed zscore/${K}mer/random/random.$i.fa zscore/${K}mer/random/random.${K}mer.$i.jf_0
done
mv zscore/${K}mer/random zscore/${K}mer/random${N}
done


