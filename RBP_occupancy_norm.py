import sys

gene_sum = {}
gene_site = {}
gene_len = {}
total_cnt = 0

# gene length
with open(sys.argv[1]) as f:
	for line in f:
		lst = line.strip().split('\t')
		gene = lst[0]
		length = int(lst[1])
		gene_len[gene] = length / 1000

with open(sys.argv[2]) as f:
	for line in f:
		lst = line.strip().split('\t')
		cnt = int(lst[6])
		gene = lst[7]
		total_cnt += cnt

		if gene in gene_sum:
			gene_sum[gene] += cnt
		else:
			gene_sum[gene] = cnt

total_cnt = total_cnt / 1000000

with open(sys.argv[2]) as f:
	for line in f:
		lst = line.strip().split('\t')
		cnt = int(lst[6])
		gene = lst[7]

		#norm_cnt = cnt / gene_sum[gene] / gene_len[gene] / total_cnt * 100
		norm_cnt = cnt / gene_sum[gene] / gene_len[gene] * 100
		print("%s\t%s\t%s\t%s\t%.6f\t%s" % (lst[0], lst[1], lst[2], gene, norm_cnt, lst[5]))
