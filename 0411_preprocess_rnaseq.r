### basic shell scripts ###
shuang@rna-garmire:~/integration/ov/tcga/tcga_ov_expr/RNASeq/BCGSC__IlluminaHiSeq_RNASeq/Level_3$ find * -maxdepth 1 -name "*.hg19.gene.quantification.txt" -print -exec cp {} ../../data \; 
shuang@rna-garmire:~/integration/ov/tcga/tcga_ov_expr$ cp FILE_SAMPLE_MAP.txt RNASeq/data   

### R scripts for ovarian cancer gene expression analysis RNASeq ###

row<-list.files()
row<-row[-c(1:2)]
map<-read.table("FILE_SAMPLE_MAP.txt",header=T,sep="\t")
map_subset<-subset(map,map[,"filename"] %in% row)
map_subset<-map_subset[match(row, map_subset$filename),]  ### match the order
colnames(map_subset)[2]<-c("barcode")

tcga_barcode<-substr(as.character(as.matrix(map_subset[,2])),1,16)

###> which(duplicated(tcga_barcode)) ### 
### [1]  141

count_matrix<-matrix(0,nrow=20806,ncol=0)
for(i in 1:length(row))
{
  temp<-read.table(row[i],sep="\t",header=T)
  count_matrix<-cbind(count_matrix,temp[,4])
  colnames(count_matrix)[i]<-tcga_barcode[i]
}

gene_ID<-as.character(as.matrix(temp[,1]))
gene_ID1<-strsplit(gene_ID,"\\|")

gene_symbol<-rep(NA,length(gene_ID))
entrez_symbol<-rep(NA,length(gene_ID))

for( i in 1:length(gene_ID))
{
  gene_symbol[i]<-gene_ID1[[i]][1]
  entrez_symbol[i]<-gene_ID1[[i]][2]
}


row.names(count_matrix)<-gene_symbol

count_matrix<-count_matrix[order(row.names(count_matrix)),]
count_matrix<-count_matrix[-c(1:32),]
#20774 genes
log2_matrix<-log2(count_matrix+1)

write.table(log2_matrix,"0411_OV_rnaseq_log2expr.txt",row.names=T,col.names=T,sep="\t",quote=F)


### batch group information for RNAseq data ###
batch_code<-map_subset[,"barcode"]
batch_code<-substr(as.character(batch_code),22,25)
batch_rnaseq<-cbind(tcga_barcode,batch_code)

write.table(batch_rnaseq,"ov_batch_rnaseq.txt",row.names=F,col.names=T,sep="\t",quote=F)

