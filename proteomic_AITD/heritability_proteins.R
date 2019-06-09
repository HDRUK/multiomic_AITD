## Tiphaine Martin - 6/6/2019
## Script developped for the paper doi: https://doi.org/10.1101/662957

## if you use this script as template, please cite this paper
## doi: https://doi.org/10.1101/662957

rm(list=ls())
library(mets) # This loads the mets package contributed to R

load("proteomics.rdata")
meta_proteomics <- read.table(file='TwinsUK_SOMAscan_Protein_info.csv',sep=",",header=T)

## remove samples having comments from SomaLogic ("Hemolyzed","Lipid")
proteomics_QC <- proteomics[-which(proteomics$SampleNotes %in% c("Hemolyzed","Lipid")),]

dat <- proteomics_QC[,c(13:1141)]
dat_continous <- dat

sex<-proteomics_QC$Gender
FamID<-proteomics_QC$FamilyID
age<-proteomics_QC$Age
zyg<-factor(proteomics_QC$ACTUALTZYG)

summary(age)
table(zyg)

res<-{}
res_noD<-{}

## for continuous values
for (i in 1:ncol(dat_continous)) {
  methy<-dat_continous[,i]
  data_tmp<-data.frame(methy,zyg,age,sex,FamID)
  data_tmp <- data_tmp[complete.cases(data_tmp),]
  data <- data_tmp
  na<- which(data_tmp$methy %in% ("NA"))
  if(length(na) > 0){
    data <- data_tmp[-na,]
  }
  data$methy <- log10(as.numeric(data$methy))
  data$age <- as.numeric(data$age)
  data$FamID <- as.numeric(data$FamID)
  data$zyg <- as.factor(data$zyg)
  num_twins <- nrow(data)
  num_fam <- table(data$zyg)
  adce <- summary(twinlm(methy~age, id="FamID", DZ="DZ", zyg="zyg",data=data,type="adce"))
  h_adce <- (adce$acde[1,]+adce$acde[2,])/(adce$acde[1,]+adce$acde[2,]+adce$acde[3,]+adce$acde[4,])
  ace <- summary(twinlm(methy~age, id="FamID", DZ="DZ", zyg="zyg",data=data,type="ace"))
  h_ace <- (ace$acde[1,])/(ace$acde[1,]+ace$acde[2,]+ace$acde[3,])
  ade <- summary(twinlm(methy~age, id="FamID", DZ="DZ", zyg="zyg",data=data,type="ade"))
  h_ade <- (ade$acde[1,]+ade$acde[2,])/(ade$acde[1,]+ade$acde[2,]+ade$acde[3,])
  dce <- summary(twinlm(methy~age, id="FamID", DZ="DZ", zyg="zyg",data=data,type="dce"))
  h_dce <- (dce$acde[2,])/(dce$acde[1,]+dce$acde[2,]+dce$acde[3,])
  ae <- summary(twinlm(methy~age, id="FamID", DZ="DZ", zyg="zyg",data=data,type="ae"))
  h_ae <- (ae$acde[1,])/(ae$acde[1,]+ae$acde[2,])
  de <- summary(twinlm(methy~age, id="FamID", DZ="DZ", zyg="zyg",data=data,type="de"))
  h_de <- (de$acde[1,])/(de$acde[1,]+de$acde[2,])
  ce <- summary(twinlm(methy~age, id="FamID", DZ="DZ", zyg="zyg",data=data,type="ce"))
  h_ce <- c(0,0,0)
  e <- summary(twinlm(methy~age, id="FamID", DZ="DZ", zyg="zyg",data=data,type="e"))
  h_e <- c(0,0,0)
  res1<-
    rbind(c("adce",adce$acde[1,],adce$estimate[23],adce$acde[2,],adce$estimate[24],adce$acde[3,],adce$estimate[25],adce$acde[4,],adce$estimate[26],h_adce,adce$AIC),
          c("ace",ace$acde[1,],ace$estimate[20],ace$acde[2,],ace$estimate[21], 0,0,0,1,ace$acde[3,],ace$estimate[22],h_ace,ace$AIC),
          c("ade",ade$acde[1,],ade$estimate[20], 0,0,0,1,ade$acde[2,],ade$estimate[21],ade$acde[3,],ade$estimate[22],h_ade,ade$AIC),
          c("dce",0,0,0,1,dce$acde[1,],dce$estimate[20], dce$acde[2,],dce$estimate[21],dce$acde[3,],dce$estimate[22],h_dce,dce$AIC),
          c("ae", ae$acde[1,], ae$estimate[17],0,0,0,1,0,0,0,1,ae$acde[2,],ae$estimate[18],h_ae,ae$AIC),
          c("ce", 0,0,0,1,ce$acde[1,], ce$estimate[17], 0,0,0,1,ce$acde[2,],ce$estimate[18],h_ce,ce$AIC),
          c("de", 0,0,0,1,0,0,0,1,de$acde[1,], de$estimate[17],de$acde[2,],de$estimate[18],h_de,de$AIC),
          c("e", 0,0,0,1,0,0,0,1, 0,0,0,1,e$acde[1,],e$estimate[14],h_e,e$AIC))
  res11<-res1[as.numeric(res1[,ncol(res1)])==min(as.numeric(res1[,ncol(res1)])),]
  if(is.matrix(res11)){
    res11 <- res11[1,]
  }
  res<-rbind(res,c(colnames(dat_continous)[i],"continous",num_twins,num_fam,res11))
  
  res1_noD<-
    rbind(c("adce",adce$acde[1,],adce$estimate[23],adce$acde[2,],adce$estimate[24],adce$acde[3,],adce$estimate[25],adce$acde[4,],adce$estimate[26],h_adce,adce$AIC),
          c("ace",ace$acde[1,],ace$estimate[20],ace$acde[2,],ace$estimate[21], 0,0,0,1,ace$acde[3,],ace$estimate[22],h_ace,ace$AIC),
          c("ade",ade$acde[1,],ade$estimate[20], 0,0,0,1,ade$acde[2,],ade$estimate[21],ade$acde[3,],ade$estimate[22],h_ade,ade$AIC),
          c("ae", ae$acde[1,], ae$estimate[17],0,0,0,1,0,0,0,1,ae$acde[2,],ae$estimate[18],h_ae,ae$AIC),
          c("ce", 0,0,0,1,ce$acde[1,], ce$estimate[17], 0,0,0,1,ce$acde[2,],ce$estimate[18],h_ce,ce$AIC),
          c("e", 0,0,0,1,0,0,0,1, 0,0,0,1,e$acde[1,],e$estimate[14],h_e,e$AIC))
  res11_noD<-res1_noD[as.numeric(res1_noD[,ncol(res1_noD)])==min(as.numeric(res1_noD[,ncol(res1_noD)])),]
  if(is.matrix(res11_noD)){
    res11_noD <- res11_noD[1,]
  }
  res_noD<-rbind(res_noD,c(colnames(dat_continous)[i],"continous",num_twins,num_fam,res11_noD))
  
  
}
colnames(res) <- c("Name_feature","type_feature","Number_twins","Num_DZ","Num_MZ","best_model",
                   "A_estimation","A_CI_2.50","A_CI_97.50","A_Pvalue",
                   "C_estimation","C_CI_2.50","C_CI_97.50","C_Pvalue",
                   "D_estimation","D_CI_2.50","D_CI_97.50","D_Pvalue",
                   "E_estimation","E_CI_2.50","E_CI_97.50","E_Pvalue",
                   "h2_estimation","h2_CI_2.50","h2_CI_97.50","AIC")

file_path <- "proteomics_heritability_log10.csv"
write.table(res,file=file_path,col.names=T,row.names=F,quote=F,sep=",")

colnames(res_noD) <- c("Name_feature","type_feature","Number_twins","Num_DZ","Num_MZ","best_model",
                   "A_estimation","A_CI_2.50","A_CI_97.50","A_Pvalue",
                   "C_estimation","C_CI_2.50","C_CI_97.50","C_Pvalue",
                   "D_estimation","D_CI_2.50","D_CI_97.50","D_Pvalue",
                   "E_estimation","E_CI_2.50","E_CI_97.50","E_Pvalue",
                   "h2_estimation","h2_CI_2.50","h2_CI_97.50","AIC")

file_path_noD <- "proteomics_heritability_noDalone_log10.csv"
write.table(res_noD,file=file_path_noD,col.names=T,row.names=F,quote=F,sep=",")
