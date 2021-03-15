install.packages("stringr")                        # Install stringr package
library("stringr")                                 # Load stringr
cancer_proteomes_CPTAC_itraq=read.csv('C:\\Users\\Noura\\Desktop\\4th_year_1st_Term\\Done\\R\\done\\ass3\\77_cancer_proteomes_CPTAC_itraq.csv')
clinical_data_breast_cancer=read.csv('C:\\Users\\Noura\\Desktop\\4th_year_1st_Term\\Done\\R\\done\\ass3\\clinical_data_breast_cancer.csv')
#############################################################
#rem TCGA-
for(r in 1:nrow(clinical_data_breast_cancer)){
  clinical_data_breast_cancer[r,'Complete.TCGA.ID']<-paste(substr(clinical_data_breast_cancer$Complete.TCGA.ID[r],6,7),substr(clinical_data_breast_cancer$Complete.TCGA.ID[r],9,12))
  
}
#rem .35TCGA
for(r in 4:ncol(cancer_proteomes_CPTAC_itraq)){
  names(cancer_proteomes_CPTAC_itraq)[r]<-paste(substr(names(cancer_proteomes_CPTAC_itraq)[r],1,2),substr(names(cancer_proteomes_CPTAC_itraq)[r],4,7))
  
}
#########################################################################################################
#matches 
matches=list()
for(r in names(cancer_proteomes_CPTAC_itraq)[4:ncol(cancer_proteomes_CPTAC_itraq)]){
  t=table(str_detect(clinical_data_breast_cancer$Complete.TCGA.ID, r))
  dff=as.data.frame(t)
  if(dff[2,2]>0){
    matches = append(matches,r)
  }}
clincal_mathces=clinical_data_breast_cancer[clinical_data_breast_cancer$Complete.TCGA.ID %in% matches,]
prot_matches=cancer_proteomes_CPTAC_itraq[,clincal_mathces$Complete.TCGA.ID]
prot_matches[colnames(cancer_proteomes_CPTAC_itraq)[1]]=cancer_proteomes_CPTAC_itraq[,1]
prot_matches[colnames(cancer_proteomes_CPTAC_itraq)[2]]=cancer_proteomes_CPTAC_itraq[,2]
prot_matches[colnames(cancer_proteomes_CPTAC_itraq)[3]]=cancer_proteomes_CPTAC_itraq[,3]
dim(prot_matches)
dim(clincal_mathces)

###########################################################################################################
#Subset the rows of proteomes table by removing rows with missing values
prot_matches=na.omit(prot_matches)
dim(prot_matches)
##############################################################################################################
u=1:nrow(clincal_mathces)
for( q in u){
  if(clincal_mathces$HER2.Final.Statu[q]=="Negative"){
    clincal_mathces$HER2.Final.Statu[q]="0"
  }else if(clincal_mathces$HER2.Final.Statu[q]=="Equivocal"){
    clincal_mathces$HER2.Final.Statu[q]="2"}else{
      clincal_mathces$HER2.Final.Statu[q]="1" }}
#############################################################################
#1 (1) Apply Pearson correlation coefficient between the each protein (each protein is considered as feature)
#from the protein table and the HER2 Final Status column in the clinical data table
prot_matchesT=as.data.frame(t(prot_matches[,1:77]))
dfff<-prot_matchesT[order(rownames(prot_matchesT)), ]
#dfff<-sapply(dfff, as.numeric)
dfff2<-clincal_mathces[order(clincal_mathces$Complete.TCGA.ID), ]
dim(dfff2)
corr1=cor(as.numeric(dfff2$HER2.Final.Statu),dfff,method = c("pearson"))
length(corr1)
table(dfff2$HER2.Final.Statu)
######################################################################
#2(1) Order the result from the highest correlation (regardless positive or negative) to the lowest correlation
yt=as.data.frame(corr1)
yt=as.data.frame(t(yt))
yt=data.frame(prot_matches$RefSeq_accession_number,yt)
colnames(yt)[2]="corelation"
yt <- yt[ sort(abs(as.numeric(yt$corelation)),decreasing=TRUE,index.return=TRUE)[[2]], ]

###################################################
#threshold
thewsold=0.09
thewsold=as.double(thewsold)
filtprotein=yt[yt$corelation>thewsold,]
###################################################
#split 
dfff5=cbind(dfff2$HER2.Final.Status,dfff)
dim(dfff)
dfff5$`dfff2$HER2.Final.Status`
group1=dfff5[dfff5$`dfff2$HER2.Final.Status`=="Negative",]
group2=dfff5[dfff5$`dfff2$HER2.Final.Status`=="Positive",]

group1<-group1[,-c(1)]
group2<-group2[,-c(1)]
t.test(group1,group2)
list_test=list()

for (p in 1:6934){
  temp<-t.test(as.numeric(group1[,p]),as.numeric((group2[,p]))) 
  yu=as.list(temp)
  list_test=append(list_test,yu[1])
}

list_test=as.data.frame(list_test)
list_test2=as.data.frame(t(list_test))
list_test2=cbind(prot_matches$RefSeq_accession_number,list_test2)
list_test2 <- list_test2[ sort(abs(as.numeric(list_test2$t)),decreasing=TRUE,index.return=TRUE)[[2]], ]
hist(list_test2$t)
thewsold2=-1
thewsold2=as.double(thewsold2)
mean(abs(yf$test))
median(abs(yf$test))
range(abs(yf$test))

demo_test<-list_test2[list_test2$t>thewsold2,]

length(intersect(demo_test$`prot_matches$RefSeq_accession_number`,filtprotein$prot_matches.RefSeq_accession_number))
##################################################
