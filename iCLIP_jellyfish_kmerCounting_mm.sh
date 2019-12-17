SITES=$1
K=$2
WINDOW=$3
NAME=$4

GENOME=$HOME/Gmatic7/genome/mouse/GRCm38.fa
GENOMESIZE=$HOME/Gmatic7/genome/mouse/GRCm38.fa.fai
GENES=$HOME/Gmatic7/gene/mouse/GRCm38_vM19_gene.bed

CPU=8

mkdir -p motif/$NAME/$WINDOW
sort -k5 -n -r $SITES|awk '{if($5>0)print}' > motif/$NAME/$WINDOW/crosslink.bed

if [ "$WINDOW" == "21nt.window" ]
then
    cat motif/$NAME/$WINDOW/crosslink.bed|awk '{print $1"\t"$2-10"\t"$3+10"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > motif/$NAME/$WINDOW/crosslink.window.bed
elif [ "$WINDOW" == "11nt.window" ]
then
    cat motif/$NAME/$WINDOW/crosslink.bed|awk '{print $1"\t"$2-5"\t"$3+5"\t"$4"\t.\t"$6}'|sort -k1,1 -k2,2n|uniq > motif/$NAME/$WINDOW/crosslink.window.bed
fi

bedtools getfasta -fi $GENOME -bed motif/$NAME/$WINDOW/crosslink.window.bed -name -s -fo motif/$NAME/$WINDOW/crosslink.window.fa
jellyfish count -m $K -s 100M -t $CPU motif/$NAME/$WINDOW/crosslink.window.fa -o motif/$NAME/$WINDOW/crosslink.window.${K}mer.jf
jellyfish dump motif/$NAME/$WINDOW/crosslink.window.${K}mer.jf_0 -c -t|sort > motif/$NAME/$WINDOW/crosslink.window.${K}mer.dumps

# permutation 100 times
for i in {0..24}; do  echo $i ; done|parallel --gnu "./iCLIP_jellyfish_shuffle_paralle.sh $NAME/$WINDOW $GENOME $GENOMESIZE $GENES $K {}"

cat motif/$NAME/$WINDOW/*.dumps|cut -f1|sort|uniq > motif/$NAME/$WINDOW/crosslink.window.${K}mer.list
./iCLIP_motif_enrich.R motif/$NAME/$WINDOW crosslink.window.${K}mer.list crosslink.window.${K}mer.dumps *random.${K}mer.*.dumps motif/$NAME/$WINDOW/motifs.${K}mer
grep -v 'UUU' motif/$NAME/$WINDOW/motifs.${K}mer.zscore.tsv > motif/$NAME/$WINDOW/motifs.${K}mer.zscore.filtUUU.tsv
./iCLIP_plot_zscore.R motif/$NAME/$WINDOW/motifs.${K}mer.zscore.filtUUU.tsv motif/$NAME/$WINDOW/motifs.${K}mer.zscore.filtUUU.pdf

rm motif/$NAME/$WINDOW/crosslink.window.${K}mer.list motif/$NAME/$WINDOW/crosslink.window.random.${K}mer.*.dumps motif/$NAME/$WINDOW/crosslink.window.${K}mer.dumps motif/$NAME/$WINDOW/crosslink.window.bed
