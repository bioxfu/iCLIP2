GENE=/home/xfu/Gmatic7/gene/human/GRCh38_v29_gene.bed
GENE_LEN=/home/xfu/Gmatic7/gene/human/GRCh38_v29_gene.geneLen

######
mkdir occupancy
cat peak/WT_merged.shifted.plus.bedgraph |awk '{print $1"\t"$2"\t"$3"\t.\t.\t+\t"$4}' >  occupancy/WT_binding_site_counts.bed
cat peak/WT_merged.shifted.minus.bedgraph|awk '{print $1"\t"$2"\t"$3"\t.\t.\t-\t"$4}' >> occupancy/WT_binding_site_counts.bed
cat peak/KO_merged.shifted.plus.bedgraph |awk '{print $1"\t"$2"\t"$3"\t.\t.\t+\t"$4}' >  occupancy/KO_binding_site_counts.bed
cat peak/KO_merged.shifted.minus.bedgraph|awk '{print $1"\t"$2"\t"$3"\t.\t.\t-\t"$4}' >> occupancy/KO_binding_site_counts.bed

cat results/WT_binding_site_clean.bed|awk '{print $1"\t"$2-15"\t"$3+15"\t"$4"\t"$5"\t"$6}' > occupancy/WT_binding_site_clean_window.bed
cat results/KO_binding_site_clean.bed|awk '{print $1"\t"$2-15"\t"$3+15"\t"$4"\t"$5"\t"$6}' > occupancy/KO_binding_site_clean_window.bed

bedtools intersect -a occupancy/WT_binding_site_counts.bed -b occupancy/WT_binding_site_clean_window.bed -s -u > occupancy/WT_binding_site_counts_window.bed
bedtools intersect -a occupancy/KO_binding_site_counts.bed -b occupancy/KO_binding_site_clean_window.bed -s -u > occupancy/KO_binding_site_counts_window.bed

bedtools intersect -a occupancy/WT_binding_site_counts_window.bed -b $GENE -s -wa -wb |cut -f1-7,11|sortBed|groupBy -g 1,2,3,6 -c 8 -o collapse -full|grep -v ','|cut -f1-8 >occupancy/WT_binding_site_counts_gene.bed
bedtools intersect -a occupancy/KO_binding_site_counts_window.bed -b $GENE -s -wa -wb |cut -f1-7,11|sortBed|groupBy -g 1,2,3,6 -c 8 -o collapse -full|grep -v ','|cut -f1-8 >occupancy/KO_binding_site_counts_gene.bed

python RBP_occupancy_norm.py $GENE_LEN occupancy/WT_binding_site_counts_gene.bed > occupancy/WT_binding_site_counts_gene_norm.bed
python RBP_occupancy_norm.py $GENE_LEN occupancy/KO_binding_site_counts_gene.bed > occupancy/KO_binding_site_counts_gene_norm.bed
 
## link the normalized crosslink sites with splicing sites
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS.bed   -b occupancy/KO_binding_site_counts_gene_norm.bed -wa -wb -s -F 1.0 > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_KO_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS.bed   -b occupancy/WT_binding_site_counts_gene_norm.bed -wa -wb -s -F 1.0 > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_WT_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS.bed   -b occupancy/KO_binding_site_counts_gene_norm.bed -wa -wb -s -F 1.0 > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_KO_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS.bed   -b occupancy/WT_binding_site_counts_gene_norm.bed -wa -wb -s -F 1.0 > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_WT_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS.bed -b occupancy/KO_binding_site_counts_gene_norm.bed -wa -wb -s -F 1.0 > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_KO_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS.bed -b occupancy/WT_binding_site_counts_gene_norm.bed -wa -wb -s -F 1.0 > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_WT_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS.bed -b occupancy/KO_binding_site_counts_gene_norm.bed -wa -wb -s -F 1.0 > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_KO_occu
bedtools intersect -a MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS.bed -b occupancy/WT_binding_site_counts_gene_norm.bed -wa -wb -s -F 1.0 > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_WT_occu
 
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_KO_occu   |grep -v 'BCS1L-ENSG00000074582.13' |grep -v 'MRPL52-ENSG00000172590.18'|grep -v 'RTF2-ENSG00000022277.12'|awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_KO_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_WT_occu   |grep -v 'BCS1L-ENSG00000074582.13' |grep -v 'MRPL52-ENSG00000172590.18'|grep -v 'RTF2-ENSG00000022277.12'|awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_upstream_3SS_WT_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_KO_occu   |grep -v 'BCS1L-ENSG00000074582.13' |grep -v 'MRPL52-ENSG00000172590.18'|grep -v 'RTF2-ENSG00000022277.12'|awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_KO_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_WT_occu   |grep -v 'BCS1L-ENSG00000074582.13' |grep -v 'MRPL52-ENSG00000172590.18'|grep -v 'RTF2-ENSG00000022277.12'|awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_upstream_3SS_WT_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_KO_occu |grep -v 'BCS1L-ENSG00000074582.13' |grep -v 'MRPL52-ENSG00000172590.18'|grep -v 'RTF2-ENSG00000022277.12'|awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_KO_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_WT_occu |grep -v 'BCS1L-ENSG00000074582.13' |grep -v 'MRPL52-ENSG00000172590.18'|grep -v 'RTF2-ENSG00000022277.12'|awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_enhanced_downstream_5SS_WT_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_KO_occu |grep -v 'BCS1L-ENSG00000074582.13' |grep -v 'MRPL52-ENSG00000172590.18'|grep -v 'RTF2-ENSG00000022277.12'|awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_KO_dist2occu
cat occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_WT_occu |grep -v 'BCS1L-ENSG00000074582.13' |grep -v 'MRPL52-ENSG00000172590.18'|grep -v 'RTF2-ENSG00000022277.12'|awk '{if($6=="+"){print $8-$2"\t"$11} if($6=="-"){print 149-($8-$2)"\t"$11}}'|sort -n|groupBy -g 1 -c 2 -o sum > occupancy/MeCP2_KO_RBFOX2_KI_SE_RPKM_2_silenced_downstream_5SS_WT_dist2occu

Rscript RBP_occupancy_norm_binding_sites.R

