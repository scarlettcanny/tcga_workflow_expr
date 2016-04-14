#### 0414 microRNA array compiling process ####

#### shell script ####
find * -maxdepth 1 -name "*.gene.tcga_level3.data.txt" -print -exec cp {} ../../data \;  
shuang@rna-garmire:~/integration/ov/tcga/tcga_ov_mirna$ cp FILE_SAMPLE_MAP.txt Expression-miRNA/data/
  
  
  ### R script ###
  row<-list.files()
  row<-row[-c(1:2)]
  map<-read.table("FILE_SAMPLE_MAP.txt",header=T,sep="\t")
  map_subset<-subset(map,map[,"filename"] %in% row)
  map_subset<-map_subset[match(row, map_subset$filename),]  ### match the order
  colnames(map_subset)[2]<-c("barcode")
  
  tcga_barcode<-substr(as.character(as.matrix(map_subset[,2])),1,16)
  
  
  count_matrix<-matrix(0,nrow=799,ncol=0)
  for(i in 1:length(row))
  {
    temp<-read.table(row[i],sep="\t",header=T)
    count_matrix<-cbind(count_matrix,as.numeric(as.matrix(temp[-1,2])))
    colnames(count_matrix)[i]<-tcga_barcode[i]
  }
  
  temp<-read.table(row[1],sep="\t",header=T)
  colnames(temp)<-as.character(as.matrix(temp[1,]))
  temp<-temp[-1,]
  mirna_ID<-as.character(as.matrix(temp[,1]))
  
  row.names(count_matrix)<-mirna_ID
  
  count_matrix<-count_matrix[order(row.names(count_matrix)),]
  
  write.table(count_matrix,"0415_OV_miRNAarray_matrix.txt",row.names=T,col.names=T,sep="\t",quote=F)
  
  
  ### batch group information for methylation data ###
  batch_code<-map_subset[,"barcode"]
  batch_code<-substr(as.character(batch_code),22,25)
  batch_miRNA<-cbind(tcga_barcode,batch_code)
  
  write.table(batch_miRNA,"ov_batch_miRNAarray.txt",row.names=F,col.names=T,sep="\t",quote=F)