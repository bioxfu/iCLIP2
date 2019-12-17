WINDOW=$1
GENOME=$2
GENOMESIZE=$3
EXONS=$4
K=$5
P=$6

for i in {0..3}
do
    echo $WINDOW $P.$i
    
    bedtools shuffle -i motif/$WINDOW/crosslink.window.bed -g $GENOMESIZE -incl $EXONS > motif/$WINDOW/crosslink.window.random.$P.$i.bed
    
    bedtools getfasta -fi $GENOME -bed motif/$WINDOW/crosslink.window.random.$P.$i.bed -name -s -fo motif/$WINDOW/crosslink.window.random.$P.$i.fa
    
    jellyfish count -m $K -s 100M motif/$WINDOW/crosslink.window.random.$P.$i.fa -o motif/$WINDOW/crosslink.window.random.${K}mer.$P.$i.jf
    
    jellyfish dump motif/$WINDOW/crosslink.window.random.${K}mer.$P.$i.jf_0 -c -t |sort > motif/$WINDOW/crosslink.window.random.${K}mer.$P.$i.dumps
    
    rm motif/$WINDOW/crosslink.window.random.$P.$i.bed motif/$WINDOW/crosslink.window.random.$P.$i.fa motif/$WINDOW/crosslink.window.random.${K}mer.$P.$i.jf_0
done

