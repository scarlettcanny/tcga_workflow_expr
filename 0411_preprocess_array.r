### R scripts for ovarian cancer gene expression analysis ###
row<-list.files()
row<-row[-c(1,2)]
file_manifest<-read.table("file_manifest.txt",header=T,sep="\t")
subfile<-subset(file_manifest,file_manifest[,"File.Name"] %in% row)

subfile<-subfile[match(row, subfile$File.Name),]

tcga_barcode<-as.character(as.matrix(subfile[,"Sample"]))

### duplicate sample #######################
> which(tcga_barcode==c("TCGA-23-1023-01"))
[1] 105 537
############################################

count_matrix<-matrix(0,nrow=17814,ncol=0)

for(i in 1:length(row))
{
  temp<-read.table(row[i],sep="\t",header=T,quote=NULL, comment='')
  temp<-temp[-1,]
  temp[,2]<-as.numeric(as.matrix(temp[,2]))
  count_matrix<-cbind(count_matrix,as.matrix(as.numeric(temp[,2])))
  colnames(count_matrix)[i]<-tcga_barcode[i]
}

gene_ID<-as.character(as.matrix(temp[,1]))
row.names(count_matrix)<-gene_ID

count_matrix1<-count_matrix

### use mean values to replace the NA values
for(i in 1:nrow(count_matrix)){
  count_matrix1[i,is.na(count_matrix1[i,])] <- mean(count_matrix1[i,], na.rm = TRUE)
}

count_matrix1<-count_matrix1[order(row.names(count_matrix1)),]
write.table(count_matrix1,"0411_OV_array_expr.txt",row.names=T,col.names=T,sep="\t",quote=F)

#count_matrix1<-count_matrix[complete.cases(count_matrix),] ### removing NA containing rows ###
#count_matrix1<-count_matrix1[order(row.names(count_matrix1)),]
#write.table(count_matrix1,"0411_OV_array_expr.txt",row.names=T,col.names=T,sep="\t",quote=F)



#gene_ID1<-strsplit(gene_ID,"\\|")
#gene_symbol<-rep(NA,length(gene_ID))
#entrez_symbol<-rep(NA,length(gene_ID))
#for( i in 1:length(gene_ID))
#{
#  gene_symbol[i]<-gene_ID1[[i]][1]
#  entrez_symbol[i]<-gene_ID1[[i]][1]
#}
#count_matrix<-count_matrix[-c(1:29),]  ## delete the anonymous genes ##
#log2_matrix<-log2(count_matrix+1)



### batch group information for array data ###
batch_code<-subfile[,"Barcode"]
batch_code<-substr(as.character(batch_code),22,25)
batch_array<-cbind(tcga_barcode,batch_code)

write.table(batch_array,"ov_batch_array.txt",row.names=F,col.names=T,sep="\t",quote=F)
