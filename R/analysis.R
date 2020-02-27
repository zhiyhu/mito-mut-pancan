mut <- read.delim("../data/TCGA-OV.mutect2_snv.tsv", as.is = T)
head(mut)

## retrieve gene list 
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations

listFilters(ensembl)
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = 'GO:0005739', mart = ensembl)

head(gene.data)

mito.gene <- unique(gene.data$hgnc_symbol)
mito.mut <- mut[mut$gene %in% mito.gene,]
head(mito.mut )


length(unique(mut$Sample_ID))

length(unique(mito.mut$Sample_ID))
table(mito.mut$chrom)
# chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 chr21 chr22  chr3  chr4  chr5  chr6 
# 514   212   282   269    94   200   107   209   623    32   209   360   109    23   101   235   146   198   201 
# chr7  chr8  chr9  chrX 
# 213   170   179   146 

length(unique(mut$Sample_ID))
# [1] 436
length(unique(mito.mut$Sample_ID))
# [1] 431

dim(mito.mut)
# [1] 4832   11

dim(mut)
# [1] 75168    11

# write.table(mito.mut, "../data/TCGA-OV.mutect2_snv_mito.txt", sep = "\t", quote = F, row.names = F)
# 