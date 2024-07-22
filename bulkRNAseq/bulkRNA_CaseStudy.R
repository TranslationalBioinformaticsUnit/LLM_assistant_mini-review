library("edgeR")
library("biomaRt")
library("cqn")
library("sva")

#set environment
count_tables_dir<-"count_tables_human/"
count_tables<-dir(count_tables_dir, pattern="txt")
setwd("bulkRNA/")

metadata <-read.table("bulkRNA_metadata.csv",header = TRUE, sep=",")

###################################
######   1-READ COUNT TABLES ######
###################################

sample_name<-unlist(strsplit(count_tables[1], "_"))
count_table_unified <- read.table(paste(count_tables_dir,count_tables[1], sep=""))
rownames_count_table_unified <- count_table_unified[,1]
count_table_unified <- data.frame(count_table_unified[,-1])
names(count_table_unified) <- sample_name[1]

for (i in 2:length(count_tables)){
  sample_name<-unlist(strsplit(count_tables[i], "_"))
  table<-read.table(paste(count_tables_dir,count_tables[i], sep=""))
  table<-table[,-1]
  count_table_unified<-cbind(count_table_unified, table)
  colnames(count_table_unified)[i] <- sample_name[1]
}
rownames(count_table_unified) <- rownames_count_table_unified

count_table_unified<-count_table_unified[-c(1,2,3,4),]
head(count_table_unified)
rm(count_tables, sample_name, table, i, rownames_count_table_unified, count_tables_dir)

boxplot(log2(count_table_unified), range=0 , main= "Boxplot before CQN", las=2)

##########################################
######    2- NORMALIZATION CQN      ######
##########################################
dataset<-3
count_table_unified<-count_table_unified[,-which(metadata$batch==dataset)]
metadata<-metadata[-which(metadata$batch==dataset),]

count_table_unified<-count_table_unified[,which(metadata$length==76)]
metadata<-metadata[which(metadata$length==76),]
#----FILTERS BEFORE CQN-------------------------------------------
#---2.1-Filter:Delete the genes that have not been mapped to any sample
count_table_unified<-count_table_unified[which(rowSums(as.matrix(count_table_unified))!=0),]
dim(count_table_unified)
#---2.2-Filter:CPM (Keep genes that have 1 read per million in at least three samples)
cpm<-cpm(count_table_unified)
count_table_unified_cpm<-count_table_unified[which((rowSums(cpm>=1))>3), ]
dim(count_table_unified_cpm)

#---2.3-Filter:Genes wiht ensembl information (Ex.ENSG00000236269 ins't in ensembleID list)
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensembleID <- getBM(attributes=c('ensembl_gene_id','percentage_gene_gc_content' ,'start_position', 'end_position', 'hgnc_symbol', 'description'), 
                    filters= c('ensembl_gene_id'), values= rownames(count_table_unified_cpm), mart= mart)

count_table_unified_ensembl<-count_table_unified_cpm[which(rownames(count_table_unified_cpm) %in% ensembleID$ensembl_gene_id),]
ensembleID<-ensembleID[order(ensembleID$ensembl_gene_id),]
count_table_unified_ensembl_order<-count_table_unified_ensembl[order(rownames(count_table_unified_ensembl)),]
all(ensembleID$ensembl_gene_id == rownames(count_table_unified_ensembl_order))
genes_length<-ensembleID$end_position-ensembleID$start_position

dim(count_table_unified_ensembl_order)

count_table_unified.cqn<-cqn(count_table_unified_ensembl_order, lengths = genes_length ,
                             x =ensembleID$percentage_gene_gc_content, 
                             sizeFactors = colSums(as.matrix(count_table_unified_ensembl_order)) , verbose=TRUE)
RPKM.cqn<- count_table_unified.cqn$y + count_table_unified.cqn$offset

#scale matrix
RPKM.zscore<-t(scale(t(RPKM.cqn),center=TRUE, scale=TRUE))
all(colnames(RPKM.zscore) == metadata$IdRun)
