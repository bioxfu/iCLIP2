DATA=$1

#split -d -l 44891827 tmp/data.qualFilteredIDs.list 

seqtk subseq $DATA tmp/x00 | sed 's/ /#/g; s/\\//#/g' | gzip > tmp/x00.fastq.gz
seqtk subseq $DATA tmp/x01 | sed 's/ /#/g; s/\\//#/g' | gzip > tmp/x01.fastq.gz
seqtk subseq $DATA tmp/x02 | sed 's/ /#/g; s/\\//#/g' | gzip > tmp/x02.fastq.gz
seqtk subseq $DATA tmp/x03 | sed 's/ /#/g; s/\\//#/g' | gzip > tmp/x03.fastq.gz
seqtk subseq $DATA tmp/x04 | sed 's/ /#/g; s/\\//#/g' | gzip > tmp/x04.fastq.gz
seqtk subseq $DATA tmp/x05 | sed 's/ /#/g; s/\\//#/g' | gzip > tmp/x05.fastq.gz
seqtk subseq $DATA tmp/x06 | sed 's/ /#/g; s/\\//#/g' | gzip > tmp/x06.fastq.gz
seqtk subseq $DATA tmp/x07 | sed 's/ /#/g; s/\\//#/g' | gzip > tmp/x07.fastq.gz
seqtk subseq $DATA tmp/x08 | sed 's/ /#/g; s/\\//#/g' | gzip > tmp/x08.fastq.gz
seqtk subseq $DATA tmp/x09 | sed 's/ /#/g; s/\\//#/g' | gzip > tmp/x09.fastq.gz

