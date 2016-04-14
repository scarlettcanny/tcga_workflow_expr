### 0414 methylation data preprocessing ###
# shell script #

shuang@rna-garmire:~/integration/ov/tcga/tcga_ov_methylation/DNA_Methylation/JHU_USC__HumanMethylation27/Level_3$ find * -maxdepth 1 -name "jhu-usc.edu_OV.HumanMethylation27.*" -print -exec cp {} ../../data \;
shuang@rna-garmire:~/integration/ov/tcga/tcga_ov_methylation/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3$ find * -maxdepth 1 -name "jhu-usc.edu_OV.HumanMethylation450.*" -print -exec cp {} ../../data \;

### R script ###
row<-list.files()
row<-row[-c(1:2)]
map<-read.table("FILE_SAMPLE_MAP.txt",header=T,sep="\t")
map_subset<-subset(map,map[,"filename"] %in% row)
map_subset<-map_subset[match(row, map_subset$filename),]  ### match the order
colnames(map_subset)[2]<-c("barcode")

tcga_barcode<-substr(as.character(as.matrix(map_subset[,2])),1,16)

count_matrix<-matrix(0,nrow=27578,ncol=0)
for(i in 1:length(row))
{
  temp<-read.table(row[i],sep="\t",header=T)
  colnames(temp)<-as.character(as.matrix(temp[1,]))
  temp<-temp[-1,]
  count_matrix<-cbind(count_matrix,as.numeric(as.matrix(temp[,2])))
  colnames(count_matrix)[i]<-tcga_barcode[i]
}

gene_ID<-as.character(as.matrix(temp[,3]))
methyl_ID<-as.character(as.matrix(temp[,1]))

row.names(count_matrix)<-methyl_ID

count_matrix<-count_matrix[order(row.names(count_matrix)),]

count_matrix1<-count_matrix[rowSums(is.na(count_matrix))!=615, ]  ### delete rows that are NAs in all platforms
count_matrix1<-count_matrix1[rowSums(is.na(count_matrix1))!=605, ]  ### delete rows that are NAs in 27k platforms

write.table(count_matrix1,"0414_OV_methyl_matrix.txt",row.names=T,col.names=T,sep="\t",quote=F)


### batch group information for methylation data ###
batch_code<-map_subset[,"barcode"]
batch_code<-substr(as.character(batch_code),22,25)
batch_methyl<-cbind(tcga_barcode,batch_code)

write.table(batch_methyl,"ov_batch_methyl.txt",row.names=F,col.names=T,sep="\t",quote=F)