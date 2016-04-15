#### 0415 miRNA batch effect removal ####
### install packages of Combat function (for microarray) ###
source("https://bioconductor.org/biocLite.R")
biocLite("sva")
library(ggfortify)
### sample PCA plot ###
df<-iris[c(1,2,3,4)]
autoplot(prcomp(df), data = iris, colour = 'Species')

### merge duplicate samples ###
#count_matrix_dupli<-rowMeans(count_matrix[,210:211])
#count_matrix2<-cbind(count_matrix[,1:209],count_matrix_dupli)
#colnames(count_matrix2)[210]<-c("TCGA-23-1023-01A")
#count_matrix2<-cbind(count_matrix2,count_matrix[,212:587])

### miRNA array PCA plot ###
miRNA_batch<-cbind(t(count_matrix2),batch_code)
autoplot(prcomp(t(count_matrix2)), data = miRNA_batch, colour = 'batch_code')

### removing batch effect with ComBat ###
library(sva)

combat_mirna_array<-ComBat(dat = count_matrix2,batch = batch_code, mod=NULL, par.prior = TRUE, prior.plots = TRUE) ### removing batch effect

miRNA_afterbatch<-cbind(t(combat_mirna_array),batch_code)
autoplot(prcomp(t(combat_mirna_array)), data = miRNA_afterbatch, colour = 'batch_code',frame=T)
