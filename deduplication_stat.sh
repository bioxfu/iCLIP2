sampleX=$1

echo -e "\t$sampleX" > results/dedup_stat_${sampleX}.tsv

grep 'Number of input reads' mapping/${sampleX}_Log.final.out|sed -r 's/.+Number of input reads/input reads/; s/ \|	/\t/' >> results/dedup_stat_${sampleX}.tsv

grep 'Uniquely mapped reads number' mapping/${sampleX}_Log.final.out|sed -r 's/.+Uniquely mapped reads number/uniquely mapped reads/; s/ \|	/\t/' >> results/dedup_stat_${sampleX}.tsv

grep 'Reads: Input Reads' dedup/${sampleX}_Aligned.sortedByCoord.out.duprm.log|sed -r 's/.+Reads: Input Reads: /deduplication input reads\t/' >> results/dedup_stat_${sampleX}.tsv

grep 'Number of reads out' dedup/${sampleX}_Aligned.sortedByCoord.out.duprm.log|sed -r 's/.+Number of reads out: /deduplication output reads\t/' >> results/dedup_stat_${sampleX}.tsv

