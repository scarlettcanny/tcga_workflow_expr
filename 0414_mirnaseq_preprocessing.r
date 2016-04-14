#### 0414 microRNA compiling process ####

#### shell script ####
find * -maxdepth 1 -name "*.hg19.mirbase20.mirna.quantification.txt" -print -exec cp {} ../../data \;  
shuang@rna-garmire:~/integration/ov/tcga/tcga_ov_mirna$ cp FILE_SAMPLE_MAP.txt miRNASeq/data/

  
  ### R script ###
  row<-list.files()
  row<-row[-c(1:2)]
  map<-read.table("FILE_SAMPLE_MAP.txt",header=T,sep="\t")
  map_subset<-subset(map,map[,"filename"] %in% row)
  map_subset<-map_subset[match(row, map_subset$filename),]  ### match the order
  colnames(map_subset)[2]<-c("barcode")
  
  tcga_barcode<-substr(as.character(as.matrix(map_subset[,2])),1,16)
  
  count_matrix<-matrix(0,nrow=1870,ncol=0)
  for(i in 1:length(row))
  {
    temp<-read.table(row[i],sep="\t",header=T)
    count_matrix<-cbind(count_matrix,as.numeric(as.matrix(temp[,3])))
    colnames(count_matrix)[i]<-tcga_barcode[i]
  }
  
  mirna_ID<-as.character(as.matrix(temp[,1]))
  
  row.names(count_matrix)<-mirna_ID
  
  count_matrix<-count_matrix[order(row.names(count_matrix)),]
  
  log2_miRNA<-log2(count_matrix+1)
  write.table(log2_miRNA,"0415_OV_miRNASeq_log2_matrix.txt",row.names=T,col.names=T,sep="\t",quote=F)
  
  
  ### batch group information for methylation data ###
  batch_code<-map_subset[,"barcode"]
  batch_code<-substr(as.character(batch_code),22,25)
  batch_miRNA<-cbind(tcga_barcode,batch_code)
  
  write.table(batch_miRNA,"ov_batch_miRNASeq.txt",row.names=F,col.names=T,sep="\t",quote=F)