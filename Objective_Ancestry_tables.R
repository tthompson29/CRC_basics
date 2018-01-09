############################################Objective ancestry Tables################################################

######################Reading in pertinant data######################
### These tables have all of the genes not just the ones specific to the mutation I'm interested in I need the short talbes
SNVH_fn<- readRDS("/Volumes/home/greally-lab/T_Trial/SNVHIGH_fn.Rds")
dim(SNVH_fn)
SNVL_fn<- readRDS("/Volumes/home/greally-lab/T_Trial/SNVLOW_fn.Rds")
SNVM_fn<- readRDS("/Volumes/home/greally-lab/T_Trial/SNVMOD_fn.Rds")
INDELHIGH_tfn<- readRDS("/Volumes/home/greally-lab/T_Trial/INDELHIGH_fn.Rds")
INDELLow_tfn<- readRDS("/Volumes/home/greally-lab/T_Trial/INDELLOW_fn.Rds")
INDELMOD_tfn<- readRDS("/Volumes/home/greally-lab/T_Trial/INDELMOD_fn.Rds")
#####################Reading in mutation specific tables##############################
SNVH_n<-readRDS("/Volumes/home/greally-lab/T_Trial/SNVHIGH_n.Rds")
INDELH_n<-readRDS("/Volumes/home/greally-lab/T_Trial/INDELHIGH_n.Rds") # *_n tables are named and mutation specific

META_data = read.table("/Volumes/home/greally-lab/T_Trial/TCGA_data_CRC/CRC_clinical.txt")
OB_AFR = read.table("/Volumes/home/greally-lab/T_Trial/OB_AFR_300")
sum(nrow(OB_AFR))# only 39 patients(these are the observed from the 300 patients)
OB_EUR = read.table("/Volumes/home/greally-lab/T_Trial/OB_EUR_300")
sum(nrow(OB_EUR))#199 patients
OB_EAS = read.table("/Volumes/home/greally-lab/T_Trial/OB_EAS_300")
OB_SAS_AMR = read.table("/Volumes/home/greally-lab/T_Trial/OB_SAS_AMR_300")
OB_ADMIX= read.table("/Volumes/home/greally-lab/T_Trial/OB_ADMIX_300")
######################################Making shorter tables#############################################

sum(OB_EUR$pcolor =="UNK") # 46 of the OB_EUR were unkown
sum(OB_EUR$pcolor =="EUR")#153 already self identified as CEU/EUR
sum(OB_EUR$pcolor == "HIS")# 0 none
META_data$bcr_patient_barcode = as.factor(META_data$bcr_patient_barcode)
OB_EUR$bcr_patient_barcode = as.factor(OB_EUR$bcr_patient_barcode)
OB_EUR_short= cbind(OB_EUR[1:4],OB_EUR[15],OB_EUR[16],OB_EUR[20])
OB_AFR_short= cbind(OB_AFR[1:4],OB_AFR[15],OB_AFR[16],OB_AFR[20])
OB_EAS_short = cbind(OB_EAS[1:4],OB_EAS[15],OB_EAS[16],OB_EAS[20])
OB_SAS_AMR_short= cbind(OB_SAS_AMR[1:4],OB_SAS_AMR[15],OB_SAS_AMR[16],OB_SAS_AMR[20])
OB_ADMIX_short= cbind(OB_ADMIX[1:4],OB_ADMIX[15],OB_ADMIX[16],OB_ADMIX[20])
View(OB_AFR_short)

#########Unecessary tables#######
META_OB_CEU = merge(x = OB_EUR_short, y= META_data, by= "bcr_patient_barcode" ,all.x = T, all.y = F) 
Meta_OB_AFR = merge(x = OB_AFR_short, y= META_data, by = "bcr_patient_barcode" ,all.x = T, all.y = F) 
Meta_OB_EAS =merge(x = OB_EAS_short, y= META_data, by = "bcr_patient_barcode" ,all.x = T, all.y = F)
Meta_OB_SAS_AMR =merge(x = OB_SAS_AMR_short, y= META_data, by = "bcr_patient_barcode" ,all.x = T, all.y = F)
Meta_OB_ADMIX =merge(x = OB_ADMIX_short, y= META_data, by = "bcr_patient_barcode" ,all.x = T, all.y = F)

OB_AFR_bcr_SUB = OB_AFR_short[,6:7]

########################Making table for SNV table organized by patient ancestry#########################################
SNVH_n_transposed = t(SNVH_n)
dim(SNVH_n_transposed) #658 6645
SNVH_n_transposed = as.data.frame(SNVH_n_transposed)
SNVH_n_transposed[1:5,1:10]
row.names(SNVH_n_transposed)
SNVH_n_transposed$patients = NA
SNVH_n_transposed$patients = paste0("TCGA-",row.names(SNVH_n_transposed))
SNVH_n_transposed[1:5,6640:6646]
dim(SNVH_n_transposed)
SNVH_n_OB_AFR_300 = NULL
SNVH_n_OB_AFR_300 = merge(x= OB_AFR_short[,6:7], y= SNVH_n_transposed, by.x = "bcr_patient_barcode", by.y = "patients", all.x = F, all.y = F)
dim(SNVH_n_OB_AFR_300) #48 6647( one of the patients is missing from the gene tables but A6-6650 is dupicated several times)
colnames(SNVH_n_OB_AFR_300[3:6647]) = SNVH_n_transposed[1,1:6645]
SNVH_n_OB_AFR_300 = SNVH_n_OB_AFR_300[!duplicated(SNVH_n_OB_AFR_300),]
dim(SNVH_n_OB_AFR_300) #41 6647
SNVH_n_OB_AFR_300[1:41,1:20]

SNVH_n_OB_AFR_300[9,6643:6647]
head(SNVH_n)
SNVH_n[1:5,1:5]
which(SNVH_n_transposed$patients =="TCGA-A6-6650")


SNVH_n_OB_AFR_300 <-as.data.frame(SNVH_n[1])
SNVH_n_OB_CEU_300 <- as.data.frame(SNVH_n[1])
SNVH_n_OB_EAS_300 <- as.data.frame(SNVH_n[1])
SNVH_n_OB_SAS_AMR_300 <- as.data.frame(SNVH_n[1])
SNVH_n_OB_ADMIX_300 <- as.data.frame(SNVH_n[1])

options(stringsAsFactors = FALSE)

class(SNVH_n[,5])

x<-16

for(x in 2:ncol(SNVH_n)){
  if(match(paste0("TCGA-",colnames(SNVH_n)[x]),OB_AFR_short$bcr_patient_barcode, nomatch = 0) > 0){
    tempcol <- as.data.frame(SNVH_n[,x],stingsAsFactors = FALSE)
    cname<- colnames(SNVH_n)[x]
    colnames(tempcol) <- cname
    SNVH_n_OB_AFR_300 <- cbind(SNVH_n_OB_AFR_300,tempcol)
  }
  if(match(paste0("TCGA-",colnames(SNVH_n)[x]),OB_EUR_short$bcr_patient_barcode,nomatch = 0) >0 ){
    tempcol <- as.data.frame(SNVH_n[,x])
    cname <- colnames(SNVH_n)[x]
    colnames(tempcol) <- cname
    SNVH_n_OB_CEU_300 <- cbind(SNVH_n_OB_CEU_300,tempcol)
  }
  if(match(paste0("TCGA-",colnames(SNVH_n)[x]),OB_EAS_short$bcr_patient_barcode,nomatch = 0) >0){
    tempcol <- as.data.frame(SNVH_n[,x])
    cname<- colnames(SNVH_n)[x]
    colnames(tempcol) <- cname
    SNVH_n_OB_EAS_300 <- cbind(SNVH_n_OB_EAS_300,tempcol)
  }
  if(match(paste0("TCGA-",colnames(SNVH_n)[x]),OB_ADMIX_short$bcr_patient_barcode,nomatch = 0) >0 ){
    tempcol <- as.data.frame(SNVH_n[,x])
    cname <- colnames(SNVH_n)[x]
    colnames(tempcol) <- cname
    SNVH_n_OB_ADMIX_300 <- cbind(SNVH_n_OB_ADMIX_300,tempcol)
  }
  if(match(paste0("TCGA-",colnames(SNVH_n)[x]),OB_SAS_AMR_short$bcr_patient_barcode,nomatch = 0) >0 ){
    tempcol <- as.data.frame(SNVH_n[,x])
    cname <- colnames(SNVH_n)[x]
    colnames(tempcol) <- cname
    SNVH_n_OB_SAS_AMR_300 <- cbind(SNVH_n_OB_SAS_AMR_300,tempcol)
  }
} 

dim(SNVH_n_OB_ADMIX_300) #6645 7
dim(SNVH_n_OB_AFR_300) #6645 43
dim(SNVH_n_OB_CEU_300) #6645 231
dim(SNVH_n_OB_EAS_300) # 6645 8
dim(SNVH_n_OB_SAS_AMR_300) #6645 43
head(SNVH_n_OB_AFR_300) # one additional is for the gene names in each of theese samples
unique(row.names(t(SNVH_n_OB_SAS_AMR_300)))
unique(row.names(t(SNVH_n_OB_ADMIX_300)))
unique(OB_SAS_AMR_short$bcr_patient_barcode)
unique(row.names(t(SNVH_n_OB_CEU_300)))
unique(row.names(t(SNVH_n_OB_EAS_300)))
unique(row.names(t(SNVH_n_OB_AFR_300)))
################Removing duplicated rows#####################
SNVH_n_OB_CEU_300_nd = SNVH_n_OB_CEU_300[,which(!duplicated(colnames(SNVH_n_OB_CEU_300)))] #6645 195(still 30 additional patients I know some have multiple runs)
SNVH_n_OB_AFR_300_nd = SNVH_n_OB_AFR_300[,which((!duplicated(colnames(SNVH_n_OB_AFR_300))))] #6645 38

SNVH_n_OB_ADMIX_300_nd = SNVH_n_OB_ADMIX_300[,which((!duplicated(colnames(SNVH_n_OB_ADMIX_300))))] #6645 7 removed any multiple
SNVH_n_OB_EAS_300_nd = SNVH_n_OB_EAS_300[,which((!duplicated(colnames(SNVH_n_OB_EAS_300))))]# 6645 8 no duplicates
SNVH_n_OB_SAS_AMR_300_nd = SNVH_n_OB_SAS_AMR_300[,which(!duplicated(colnames(SNVH_n_OB_SAS_AMR_300)))] #6645 43
sum(ncol(SNVH_n_OB_ADMIX_300_nd),ncol(SNVH_n_OB_AFR_300_nd),ncol(SNVH_n_OB_CEU_300_nd),ncol(SNVH_n_OB_EAS_300_nd), ncol(SNVH_n_OB_SAS_AMR_300_nd)- 5)
#286 columns when you remove all the duplicates. Some of the samples that are in the whole data set do not have  
##############################Merge races in order######################################
SNVH_OB_Ordered <- cbind(SNVH_n_OB_AFR_300_nd,SNVH_n_OB_CEU_300_nd[,2:195],SNVH_n_OB_EAS_300_nd[2:8],SNVH_n_OB_SAS_AMR_300_nd[,2:43], SNVH_n_OB_ADMIX_300_nd[,2:7])
rownames(SNVH_OB_Ordered) <- SNVH_OB_Ordered$GENES
SNVH_OB_Ordered[1:5,1:5]
SNVH_OB_Ordered_short <-SNVH_OB_Ordered[-1]
write.table(SNVH_OB_Ordered_short,"/Volumes/home/greally-lab/T_Trial/SVNH_300_ordered.txt",sep = "\t")
dim(SNVH_OB_Ordered_short)
#6645 286 - missing 3 samples. Still not sure what happened = SNVH_n(16,359,624)
colnames(SNVH_n_OB_ADMIX_300_nd[2]) #A6-6652
colnames(SNVH_n_OB_AFR_300_nd[2])#A6-A565
colnames(SNVH_n_OB_CEU_300_nd[2])#G4-6627
colnames(SNVH_n_OB_EAS_300_nd[2])# CA-6717
colnames(SNVH_n_OB_SAS_AMR_300_nd[2])#AA-3531

saveRDS(SNVH_OB_Ordered_short, file ="/Volumes/home/greally-lab/T_Trial/SVNH_ordered_300.Rds")

#################################################Permuation############################################################

AFR_mean_OB_300_SNVH = NULL
CEU_mean_OB_300_SNVH = NULL
SAS_AMR_mean_OB_300_SNVH = NULL
EAS_mean_OB_300_SNVH = NULL
ADMIX_mean_OB_300_SNVH = NULL


#read in the data each time
#SNVH_OB_Ordered_short <- read.table("/Volumes/home/greally-lab/T_Trial/SVNH_ordered.txt")
dim(SNVH_OB_Ordered_short) #6645 286
SNVH_OB_Ordered_short[1:5,1:5]
#row.names(SNVH_OB_Ordered_short) <- SNVH_OB_Ordered_short$GENES
#SNVH_OB_Ordered_short <-SNVH_OB_Ordered_short[-1]

class(SNVH_OB_Ordered_short[1,1])
which(class(SNVH_OB_Ordered_short)!= "numeric")
sum(SNVH_OB_Ordered_short[,2])
#apply(SNVH_OB_Ordered_short,2,as.numeric(SNVH_OB_Ordered_short))
#SNVH_OB_Ordered_short[,2] <- as.numeric(levels(SNVH_OB_Ordered_short[,2]))[SNVH_OB_Ordered_short[,2]]

options(stringsAsFactors = F)
set.seed(1)
for (i in 1:1000){
  if(i%%100==0) print(i)
  #remember to skip the first column it is genes
  tmp_data = t(apply(SNVH_OB_Ordered_short,1,function(x){as.numeric(sample(x))}))
  tmp_data = as.data.frame(tmp_data)
  #class(tmp_data[,1])
  
  # here instead of collect the permuated data everytime,
  # we will only collect statistical values to save memory
  # we will collect median/mean of number of SNVs in that gene
  # in each population (or other values you are interested in) 
  
  sub_AFR_mean_OB_300_SNVH = apply(tmp_data[,1:37],1,mean)
  sub_CEU_mean_OB_300_SNVH = apply(tmp_data[,38:231],1,mean)
  sub_EAS_mean_OB_300_SNVH = apply(tmp_data[,232:238],1,mean)
  sub_SAS_AMR_mean_OB_300_SNVH = apply(tmp_data[,239:280],1,mean)
  sub_ADMIX_mean_OB_300_SNVH = apply(tmp_data[,281:286],1,mean)
  ###Testing to see some sample sums
  #sum(sub_ASIAN_mean) 23 .02066
  #sum(sub_ASIAN_mean) = 25
  #sum(sub_CEU_mean) = 23.21429
  
  AFR_mean_OB_300_SNVH = cbind(AFR_mean_OB_300_SNVH, sub_AFR_mean_OB_300_SNVH)
  CEU_mean_OB_300_SNVH = cbind(CEU_mean_OB_300_SNVH, sub_CEU_mean_OB_300_SNVH)
  EAS_mean_OB_300_SNVH = cbind(EAS_mean_OB_300_SNVH,sub_EAS_mean_OB_300_SNVH)
  ADMIX_mean_OB_300_SNVH = cbind(ADMIX_mean_OB_300_SNVH, sub_ADMIX_mean_OB_300_SNVH)
  SAS_AMR_mean_OB_300_SNVH = cbind(SAS_AMR_mean_OB_300_SNVH, sub_SAS_AMR_mean_OB_300_SNVH)
}


## the result say AFR_mean
## each column is the statistics from one permutation, each row is the sampled distribution of randomness
# Commented out the row names becase we care about the gene names
AFR_mean_OB_300_SNVH= as.data.frame(AFR_mean_OB_300_SNVH)
CEU_mean_OB_300_SNVH= as.data.frame(CEU_mean_OB_300_SNVH)
EAS_mean_OB_300_SNVH= as.data.frame(EAS_mean_OB_300_SNVH)
SAS_AMR_mean_OB_300_SNVH= as.data.frame(SAS_AMR_mean_OB_300_SNVH)
ADMIX_mean_OB_300_SNVH= as.data.frame(ADMIX_mean_OB_300_SNVH)

colnames(AFR_mean_OB_300_SNVH) = paste0("permute",seq(1:1000))
AFR_mean_OB_300_SNVH[1:10,1:10]
rownames(AFR_mean_OB_300_SNVH) = row.names(SNVH_OB_Ordered_short)
colnames(CEU_mean_OB_300_SNVH) = paste0("permute",seq(1:1000))
CEU_mean_OB_300_SNVH[1:10,1:10]
rownames(CEU_mean_OB_300_SNVH) = row.names(SNVH_OB_Ordered_short)
colnames(EAS_mean_OB_300_SNVH) = paste0("permute",seq(1:1000))
rownames(EAS_mean_OB_300_SNVH) = row.names(SNVH_OB_Ordered_short)
colnames(SAS_AMR_mean_OB_300_SNVH) = paste0("permute",seq(1:1000))
rownames(SAS_AMR_mean_OB_300_SNVH) = row.names(SNVH_OB_Ordered_short)
colnames(ADMIX_mean_OB_300_SNVH) = paste0("permute",seq(1:1000))
rownames(ADMIX_mean_OB_300_SNVH) = row.names(SNVH_OB_Ordered_short)
head(ADMIX_mean_OB_300_SNVH)
class(AFR_mean_OB_300_SNVH[1,])

SNVH_OB_Ordered_short= as.data.frame(SNVH_OB_Ordered_short)
class(SNVH_OB_Ordered_short[3,])
### next we gonna compute the observation statistics we have
ob_AFR_mean_OB_300_SNVH = apply(SNVH_OB_Ordered_short[,1:37],1,function(x){mean(as.numeric(x))})
ob_CEU_mean_OB_300_SNVH = apply(SNVH_OB_Ordered_short[,38:231],1,function(x){mean(as.numeric(x))})
ob_EAS_mean_OB_300_SNVH = apply(SNVH_OB_Ordered_short[,232:238],1,function(x){mean(as.numeric(x))})
ob_SAS_AMR_mean_OB_300_SNVH = apply(SNVH_OB_Ordered_short[,239:280],1,function(x){mean(as.numeric(x))})
ob_ADMIX_mean_OB_300_SNVH = apply(SNVH_OB_Ordered_short[,281:286],1,function(x){mean(as.numeric(x))})

### so we can compute the p value now
p_AFR_300_enrichment = sapply(1:6645,function(x)sum(ob_AFR_mean_OB_300_SNVH[x]<AFR_mean_OB_300_SNVH[x,])/500)
p_CEU_300_enrichment = sapply(1:6645,function(x)sum(ob_CEU_mean_OB_300_SNVH[x]<CEU_mean_OB_300_SNVH[x,])/500)
p_EAS_300_enrichment = sapply(1:6645,function(x)sum(ob_EAS_mean_OB_300_SNVH[x]<EAS_mean_OB_300_SNVH[x,])/500)
p_SAS_AMR_300_enrichment = sapply(1:6645,function(x)sum(ob_SAS_AMR_mean_OB_300_SNVH[x]<SAS_AMR_mean_OB_300_SNVH[x,])/500)
p_ADMIX_300_enrichment = sapply(1:6645,function(x)sum(ob_ADMIX_mean_OB_300_SNVH[x]<ADMIX_mean_OB_300_SNVH[x,])/500)
p_AFR_300_enrichment[1:20]
# we will see if where the enrichment are, you can set your sinificant level, 
# here I will use 1 out of 500


#seperate actual enrichment from simply zero, moved the p value up to p of 0.02 for most of the samples except the CEU which for some reason hadno signigifant values. I think this is because it is the majority of the table.
SNVH_AFR_enrich_300 <- which(p_AFR_300_enrichment< 10/500 & p_AFR_300_enrichment > 0/1)
length(SNVH_AFR_enrich_300)
#291 gnes of interest
SNVH_CEU_enrich_300 <-which(p_CEU_300_enrichment<100/500 & p_CEU_300_enrichment > 0/1) 
SNVH_EAS_enrich_300 <-which(p_EAS_300_enrichment<10/500 & p_EAS_300_enrichment > 0/1)
SNVH_SAS_AMR_enrich_300 <-which(p_SAS_AMR_300_enrichment<50/500 & p_SAS_AMR_300_enrichment >0/1)
SNVH_ADMIX_enrich_300 <-which(p_ADMIX_300_enrichment<10/500 & p_ADMIX_300_enrichment >0/1)
#SNVH_AFR_enrich_short_genes <- rownames(SNVH_OB_Ordered_short)[which(p_AFR_enrichment<1/500)]
#SNVH_CEU_enrich_short_genes <- rownames(SNVH_OB_Ordered_short)[which(p_CEU_enrichment<1/500)]
#SNVH_ASIAN_enrich_short_genes <- rownames(SNVH_OB_Ordered_short)[which(p_ASIAN_enrichment<1/500)]
#SNVH_OTHER_enrich_short_genes <- rownames(SNVH_OB_Ordered_short)[which(p_OTHER_enrichment<1/500)]

##############  Save enrichment table####################
saveRDS(SNVH_AFR_enrich_300, file = "/Volumes/Macintosh HD/Users/TVictoria/Documents/Git/CRC_basics/Tables/SNVH_AFR_300_enriched.Rds")
saveRDS(SNVH_CEU_enrich_300, file = "/Volumes/Macintosh HD/Users/TVictoria/Documents/Git/CRC_basics/Tables/SNVH_CEU_300_enriched.Rds")
saveRDS(SNVH_EAS_enrich_300, file = "/Volumes/Macintosh HD/Users/TVictoria/Documents/Git/CRC_basics/Tables/SNVH_EAS_300_enriched.Rds")
saveRDS(SNVH_SAS_AMR_enrich_300, file = "/Volumes/Macintosh HD/Users/TVictoria/Documents/Git/CRC_basics/Tables/SNVH_SAS_AMR_300_enriched.Rds")
saveRDS(SNVH_ADMIX_enrich_300, file = "/Volumes/Macintosh HD/Users/TVictoria/Documents/Git/CRC_basics/Tables/SNVH_ADMIX_300_enriched.Rds")
rownames(SNVH_OB_Ordered_short)[129]

######       SNVH_ethnicity_enrich_short is the row number of significance for that particlar ethnicity##########
####################### Making significant tables for each group ######################
SNVH_300_sig_AFR <- as.matrix.data.frame(SNVH_OB_Ordered_short[SNVH_AFR_enrich_300,])
SNVH_300_sig_CEU <- as.matrix.data.frame(SNVH_OB_Ordered_short[SNVH_CEU_enrich_300,])
SNVH_300_sig_EAS <- as.matrix.data.frame(SNVH_OB_Ordered_short[SNVH_EAS_enrich_300,])
SNVH_300_sig_SAS_AMR <- as.matrix.data.frame(SNVH_OB_Ordered_short[SNVH_SAS_AMR_enrich_300,])
SNVH_300_sig_ADMIX <- as.matrix.data.frame(SNVH_OB_Ordered_short[SNVH_ADMIX_enrich_300,])
head(SNVH_300_sig_AFR)

SNVH_300_sig_AFR[1:9,1:5]
heatmap(SNVH_300_sig_AFR)


class(SNVH_300_sig_AFR[,1])

apply(SNVH_300_sig_AFR,1,sum)
which(sapply(SNVH_300_sig_AFR,is.numeric)==FALSE)

#1805 (number of genes of interest for CEU)
write.table(SNVH_300_sig_AFR, file = "/Volumes/Macintosh HD/Users/TVictoria/Documents/Git/CRC_basics/Tables/SNVH_300_sig_table_AFR.txt", sep = "\t")
write.table(SNVH_300_sig_CEU, file = "/Volumes/Macintosh HD/Users/TVictoria/Documents/Git/CRC_basics/Tables/SNVH_300_sig_table_CEU.txt", sep = "\t")
write.table(SNVH_300_sig_EAS, file = "/Volumes/Macintosh HD/Users/TVictoria/Documents/Git/CRC_basics/Tables/SNVH_300_sig_table_EAS.txt", sep = "\t")
write.table(SNVH_300_sig_SAS_AMR, file = "/Volumes/Macintosh HD/Users/TVictoria/Documents/Git/CRC_basics/Tables/SNVH_300_sig_table_SAS_AMR.txt", sep = "\t")
write.table(SNVH_300_sig_ADMIX, file = "/Volumes/Macintosh HD/Users/TVictoria/Documents/Git/CRC_basics/Tables/SNVH_300_sig_table_ADMIX.txt", sep = "\t")

apply(SNVH_sig_CEU,1,sum)
which(sapply(SNVH_sig_CEU,is.numeric)==FALSE)

# let's take a look at the index I set to be enriched in CEU again
which(p_AFR_enrichment<1/500)
which(p_OTHER_enrichment<1/500)


##############################Double checking the values of the tables####################################
############### Variables for table
SNVH_AFR = read.table("/Volumes/home/greally-lab/T_Trial/SNVH_p002_AFR.txt")
SNVH_CEU = read.table("/Volumes/home/greally-lab/T_Trial/SNVH_p002_CEU.txt")
SNVH_ASIAN = read.table("/Volumes/home/greally-lab/T_Trial/SNVH_p002_ASIAN.txt")
SNVH_OTHER = read.table("/Volumes/home/greally-lab/T_Trial/SNVH_p002_OTHER.txt")

rownames(SNVH_OTHER) <- SNVH_OTHER$GENES

SNVH_OTHER <- SNVH_OTHER[-1]
SNVH_OTHER[1:10,1:10]


rownames(SNVH_CEU) = SNVH_CEU$GENES
SNVH_CEU = SNVH_CEU[-1]
SNVH_CEU[1:10,1:10]

dim(SNVH_n_AFR)
#6645 71(1st column gene names)
dim(SNVH_sig_AFR)
#291 71
dim(SNVH_CEU)
#1805 329

## Check that at least 2 people have the gene
SNVH_CEU2plus = which(apply(SNVH_CEU,1,function(x)sum(x)>1))
#########################Making plots###################################################
##black 1- 37 white 38 to 231
library(ComplexHeatmap)
W_B_300_AFR_HM = as.data.frame(rbind(SNVH_300_sig_AFR[,1:231],SNVH_300_sig_CEU[,1:231]))
W_B_300_AFR = as.data.frame(rbind(SNVH_300_sig_AFR[,1:231],SNVH_300_sig_CEU[,1:231]))
W_300_AFR = as.data.frame(W_B_300_AFR[,38:231])
B_300_AFR = as.data.frame(W_B_300_AFR[,1:37])
SNVH_Mat_list =list(SNVH_300_sig_AFR[,1:231], SNVH_300_sig_CEU[,1:231])
color = c( 0 = "blue", "1" = "red", "2" = "#008000")
alter_fun = list(
  "0" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  #SNVLOW = function(x, y, w, h) {
  #  grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "blue", col = NA))
  #  },
  "1" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "red", col = NA))
  },
  "2" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
  }
)

####### Changing the 0 an 1 to zero and one#######
x = 1

for(x in 1:nrow(W_B_300_AFR)){
  y = W_B_300_AFR[x,]
       W_B_300_AFR[x,] = replace(W_B_300_AFR[x,], list = which(W_B_300_AFR[x,] == 0), values = "")
       W_B_300_AFR[x,] = replace(W_B_300_AFR[x,], list = which(W_B_300_AFR[x,] == 1), values = "SNV_HIGH;")
}


#col_2 = c("NO_MUT" = "azure3", "SNV_HIGH" = "orangered3")
col_2 = c( "SNV_HIGH" = "orangered3")
alter_fun_2 = list(
  background = function(x,y,w,h){ 
    grid.rect(x, y, w-unit(0.5,"mm"), h-unit(0.5, "mm"), gp = gpar(fill ="#CCCCCC"), col = NA)
  },
  SNV_HIGH = function(x,y,w,h){
    grid.rect(x, y, w-unit(0.5,"mm"), h-unit(0.5, "mm"), gp = gpar(fill ="orangered3", col = NA))
  }
)
 
B_300_AFR = as.matrix(B_300_AFR)
W_300_AFR = as.matrix(W_300_AFR)
mat_list = list(B_300_AFR,W_300_AFR) 
rownames(W_B_300_AFR)=cat(rownames(SNVH_300_sig_AFR), rownames(SNVH_300_sig_CEU))
rownames(W_B_300_AFR) = rbind(rownames(SNVH_300_sig_AFR), rownames(SNVH_300_sig_CEU))
AFR_300 = oncoPrint(W_B_300_AFR, get_type = function(x) strsplit(x,";")[[1]],
          alter_fun = alter_fun_2, 
          col= col_2,
          #remove_empty_columns = T,
          column_title = "CRC SNV High Impact Mutations",
          heatmap_legend_param = list(title = "Mutations", at = c("SNV_HIGH"), 
                                      labels = c("SNV")))

draw(AFR_300, heatmap_legend_side = "bottom")
head(W_B_300_AFR)
Heatmap(W_B_300_AFR_HM)
colnames(W_B_300_AFR) = colnames(W_B_300_AFR_HM)
SNVH_AFR_300_
