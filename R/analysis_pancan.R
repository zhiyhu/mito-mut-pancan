
library(ggplot2)
### GO - mitochondian

## retrieve gene list 
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations

listFilters(ensembl)
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go', values = 'GO:0005739', mart = ensembl)

head(gene.data)
mito.gene <- unique(gene.data$hgnc_symbol)

filenames <- list.files("~/TCGA/", pattern = "*.tsv")

## Mutation types we will use
mut_type <- c("frameshift_variant","inframe_deletion","inframe_insertion","missense_variant",
              "missense_variant;splice_region_variant","protein_altering_variant","splice_acceptor_variant","splice_donor_variant")

# Calculate the number of mito mutations 
for(itor in 1:length(filenames)){
    data <- read.delim2(file = paste("~/TCGA/", filenames[itor], sep = ""), as.is = T)
    data$Sample_ID <- as.factor(data$Sample_ID)
    
    mut <- data
    mut <- mut[mut$filter == "PASS" & mut$effect %in% mut_type,]
    
    mito.mut <- data[data$gene %in% mito.gene,]
    mito.mut <- mito.mut[mito.mut$filter == "PASS" & mito.mut$effect %in% mut_type,]
    
    freq_all <- table(mut$Sample_ID)
    freq <- table(mito.mut$Sample_ID)
    
    df_tmp <- data.frame(cancer.type = strsplit(filenames[itor], "[.]")[[1]][1],
                         mito_freq = freq,
                         mut_freq = freq_all,
                         sample_ID =names(freq))
    
    if(itor == 1) {
        df <- df_tmp
    } else {
        df <- rbind(df, df_tmp)
    }
}
head(df)

acronym <- read.delim2(header = F, file = "~/TCGA/acronym_list.txt")
acronym$V1 <- gsub(pattern = "GDC TCGA ", replacement = "", x = acronym$V1)

df$cancer.type <- gsub(pattern = "TCGA-", replacement = "", x = df$cancer.type)
df$cancer <- acronym$V1[match(df$cancer.type, acronym$V2)]


ggplot(df, aes(x = cancer, y = freq.Freq))+ geom_jitter(alpha = 0.4, size = 0.75) +coord_flip() + theme_classic() +
    xlab("") + ylab("Frequency of mutations per sample") + 
    theme(axis.text  = element_text(size = 14), axis.title = element_text(size = 16))

# ggsave("plot/jitters_pancan_fullName.png", width = 14, height = 12, dpi = 300)


df_split <- split(df, df$cancer)
sapply(df_split, function(x) return(mean(x$freq.Freq)))
write.csv(as.data.frame(sapply(df_split, function(x) return(mean(x$mito_freq.Freq)))),"average_mito_mut.csv")

sum(df$freq.Freq)

ggplot(df, aes(x = df$mut_freq.Freq,  y = df$mito_freq.Freq)) + 
    geom_point(alpha = 0.4)+facet_wrap(~cancer,ncol = 6, scales = "free") +
    theme_classic() + xlab("Frequency of somatic mutations") + xlim(0,200) + ylim(0,25) +
    ylab("Frequency of somatic mitochondial-located mutations")
# ggsave("plot/scatter_plot_pancan.png", height = 12, width = 12)

df$percentage <- df$mito_freq.Freq/df$mut_freq.Freq


df_filtered <- df[df$mut_freq.Freq <= 200,]

ggplot(df_filtered, aes(x = mut_freq.Freq,  y = mito_freq.Freq)) + 
    geom_point(alpha = 0.4)+facet_wrap(~cancer,ncol = 6, scales = "free") +
    theme_classic() + xlab("Frequency of somatic mutations") + xlim(0,200) + ylim(0,25) +
    ylab("Frequency of somatic mitochondial-located mutations")
# ggsave("plot/scatter_plot_pancan.png", height = 12, width = 12)

median(df_filtered$percentage[df_filtered$cancer == "Kidney Clear Cell Carcinoma"])
median(df_filtered$percentage[df_filtered$cancer != "Kidney Clear Cell Carcinoma"], na.rm = T)
