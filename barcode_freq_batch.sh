zcat tmp/x00.fastq.gz | awk -v umi1_len=5 -v exp_bc_len=6 '{ if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len) }' | sort | uniq -c | sort -k1,1rn > tmp/x00.exp_barcodes.detected
zcat tmp/x01.fastq.gz | awk -v umi1_len=5 -v exp_bc_len=6 '{ if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len) }' | sort | uniq -c | sort -k1,1rn > tmp/x01.exp_barcodes.detected
zcat tmp/x02.fastq.gz | awk -v umi1_len=5 -v exp_bc_len=6 '{ if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len) }' | sort | uniq -c | sort -k1,1rn > tmp/x02.exp_barcodes.detected
zcat tmp/x03.fastq.gz | awk -v umi1_len=5 -v exp_bc_len=6 '{ if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len) }' | sort | uniq -c | sort -k1,1rn > tmp/x03.exp_barcodes.detected
zcat tmp/x04.fastq.gz | awk -v umi1_len=5 -v exp_bc_len=6 '{ if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len) }' | sort | uniq -c | sort -k1,1rn > tmp/x04.exp_barcodes.detected
zcat tmp/x05.fastq.gz | awk -v umi1_len=5 -v exp_bc_len=6 '{ if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len) }' | sort | uniq -c | sort -k1,1rn > tmp/x05.exp_barcodes.detected
zcat tmp/x06.fastq.gz | awk -v umi1_len=5 -v exp_bc_len=6 '{ if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len) }' | sort | uniq -c | sort -k1,1rn > tmp/x06.exp_barcodes.detected
zcat tmp/x07.fastq.gz | awk -v umi1_len=5 -v exp_bc_len=6 '{ if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len) }' | sort | uniq -c | sort -k1,1rn > tmp/x07.exp_barcodes.detected
zcat tmp/x08.fastq.gz | awk -v umi1_len=5 -v exp_bc_len=6 '{ if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len) }' | sort | uniq -c | sort -k1,1rn > tmp/x08.exp_barcodes.detected
zcat tmp/x09.fastq.gz | awk -v umi1_len=5 -v exp_bc_len=6 '{ if (FNR%4==2) print substr($1,(umi1_len+1),exp_bc_len) }' | sort | uniq -c | sort -k1,1rn > tmp/x09.exp_barcodes.detected

cat tmp/x0*.exp_barcodes.detected|awk '{print $2"\t"$1}'|sort -k1|groupBy -g 1 -c 2 -o sum|sort -nrk2 > results/exp_barcodes.detected
