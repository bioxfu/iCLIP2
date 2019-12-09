sampleX=$1

echo -e "\t$sampleX" > results/Xlink_stat_${sampleX}.tsv

#### Number of crosslink events, i.e. reads after duplicate removal:
cat bed/${sampleX}_Aligned.sortedByCoord.out.duprm.shifted.*.bedgraph | awk 'BEGIN{ totalcount=0 }{ totalcount += (($3 - $2) * $4) }END{ print "crosslink_events\t"totalcount }' >> results/Xlink_stat_${sampleX}.tsv

#### Number of crosslinked nucleotides, i.e. positions harbouring crosslinked nucleotides (if both strands are covered, count as 2):
cat bed/${sampleX}_Aligned.sortedByCoord.out.duprm.shifted.*.bedgraph | awk 'BEGIN{ totalpos=0 }{ totalpos += ($3 - $2) }END{ print "crosslinked_nucleotides\t"totalpos }' >> results/Xlink_stat_${sampleX}.tsv

#### Number of stacked crosslink events, i.e. crosslink events on positions with >1 crosslink events:
cat bed/${sampleX}_Aligned.sortedByCoord.out.duprm.shifted.*.bedgraph | awk 'BEGIN{ totalstackedcount=0 }{ if($4 > 1) totalstackedcount += (( $3 - $2) * $4) }END{ print "stacked_crosslink_events\t"totalstackedcount }' >> results/Xlink_stat_${sampleX}.tsv

#### Number of nucleotides with stacked crosslink events, i.e. positions with >1 crosslink events: 
cat bed/${sampleX}_Aligned.sortedByCoord.out.duprm.shifted.*.bedgraph | awk 'BEGIN{ totalstackedpos=0 }{ if($4 > 1) totalstackedpos += ($3 - $2) }END{ print "nucleotides_with_stacked_crosslink_events\t"totalstackedpos }' >> results/Xlink_stat_${sampleX}.tsv
