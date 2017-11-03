options(stringsAsFactors = FALSE)

#Individualized specil tables that don't have all the genes, only the significant genes 
SNVHIGH_table <- readRDS("/Volumes/home/greally-lab/T_Trial/Tables/snvh_table.Rds")
SNVLOW_table <- readRDS("/Volumes/home/greally-lab/T_Trial/Tables/snvl_table.Rds")
SNVMOD_table <- readRDS("/Volumes/home/greally-lab/T_Trial/Tables/snvm_table.Rds")

INDELHIGH_table<- readRDS("/Volumes/home/greally-lab/T_Trial/Tables/IDH_table.Rds")
INDELLOW_table<- readRDS("/Volumes/home/greally-lab/T_Trial/Tables/IDL_table.Rds")
INDELMOD_table<- readRDS("/Volumes/home/greally-lab/T_Trial/Tables/IDM_table.Rds")

name_table <- read.delim("/Volumes/home/greally-lab/T_Trial/File_name.txt", header = F)
View(name_table)
name_table_2 <- as.data.frame(cbind(name_table,0))


x <-2
for(x in 1:nrow(name_table)){
  if(x%%2 == 0){
    y <- x-1
    name_table_2[y,2] <- as.character(name_table[x,])
  }
}
View(name_table_2)
name_table_3 <- NULL
for(x in 1:nrow(name_table_2)){
  if(name_table_2[x,2] != 0){
    temp_row <-name_table_2[x,]
    name_table_3 <- rbind(name_table_3,temp_row)
  }
}
dim(name_table_3)
#657 2
write.csv(name_table_3,file = "/Volumes/home/greally-lab/T_Trial/Name_convert.csv")
name_table_3 <- read.csv("/Volumes/home/greally-lab/T_Trial/Name_convert.csv")
name_table_3 <- name_table_3[-1]

colnames(INDELHIGH_table_f)
INDELLow_tfn <- as.data.frame(INDELLow_table_f)
INDELHIGH_tfn <-as.data.frame(INDELHIGH_table_f)
INDELMOD_tfn <-as.data.frame(INDELMOD_table_f)
#match(colnames(INDELHIGH_table_f)[-1],name_table_3[,1])
#match(colnames(INDELMOD_table_f)[-1],name_table_3[,1])
#Test to make sure the column names are in the same order as the current columns
#These should be in the same order because I always read in the files same way each time
colnames(INDELHIGH_tfn)[-1]<-name_table_3[,2]
colnames(INDELLow_tfn)[-1]<-name_table_3[,2]
colnames(INDELMOD_tfn)[-1]<-name_table_3[,2]
head(INDELHIGH_tfn)

saveRDS(INDELHIGH_tfn, file = "/Volumes/home/greally-lab/T_Trial/INDELHIGH_fn.Rds")
saveRDS(INDELLow_tfn, file = "/Volumes/home/greally-lab/T_Trial/INDELLOW_fn.Rds")
saveRDS(INDELMOD_tfn, file = "/Volumes/home/greally-lab/T_Trial/INDELMOD_fn.Rds")

SNVH_fn <- as.data.frame(snvh_table_f)
SNVM_fn <- as.data.frame(snvl_table_f)
SNVL_fn <- as.data.frame(snvm_table_f)
#match(colnames(snvm_table_f)[-1],name_table_3[,1])
colnames(SNVH_fn)[-1]<-name_table_3[,2]
colnames(SNVL_fn)[-1]<-name_table_3[,2]
colnames(SNVM_fn)[-1]<-name_table_3[,2]
saveRDS(SNVH_fn, file = "/Volumes/home/greally-lab/T_Trial/SNVHIGH_fn.Rds")
saveRDS(SNVL_fn, file = "/Volumes/home/greally-lab/T_Trial/SNVLOW_fn.Rds")
saveRDS(SNVM_fn, file = "/Volumes/home/greally-lab/T_Trial/SNVMOD_fn.Rds")



#############REPEAT WITH INDIVIDUAL TABLES(NOT FULL) #############################
name_table_3 <- read.csv("/Volumes/home/greally-lab/T_Trial/Name_convert.csv")
name_table_3 <- name_table_3[-1]

colnames(INDELHIGH_table)
INDELLOW_tn <- as.data.frame(INDELLOW_table)
INDELHIGH_tn <-as.data.frame(INDELHIGH_table)
INDELMOD_tn <-as.data.frame(INDELMOD_table)
#match(colnames(INDELHIGH_table)[-1],paste0(name_table_3[,1],"_pINDELHIGH"))
#match(colnames(INDELMOD_table_f)[-1],name_table_3[,1])
#Test to make sure the column names are in the same order as the current columns
#These should be in the same order because I always read in the files same way each time
colnames(INDELHIGH_tn)[-1]<-name_table_3[,2]
colnames(INDELLOW_tn)[-1]<-name_table_3[,2]
colnames(INDELMOD_tn)[-1]<-name_table_3[,2]
head(INDELHIGH_tfn)

saveRDS(INDELHIGH_tn, file = "/Volumes/home/greally-lab/T_Trial/INDELHIGH_n.Rds")
saveRDS(INDELLOW_tn, file = "/Volumes/home/greally-lab/T_Trial/INDELLOW_n.Rds")
saveRDS(INDELMOD_tn, file = "/Volumes/home/greally-lab/T_Trial/INDELMOD_n.Rds")

SNVH_n <- as.data.frame(SNVHIGH_table)
SNVL_n <- as.data.frame(SNVLOW_table)
SNVM_n <- as.data.frame(SNVMOD_table)
#Make sure they are inthe same order
#match(colnames(snvh_table)[-1],paste0(name_table_3[,1],"_pSNVHIGH"))
#shows that 1=1,2=2 etc So they have the same file names which are matched to sample names in name_table_3[,2]
colnames(SNVH_n)[-1]<-name_table_3[,2]
colnames(SNVL_n)[-1]<-name_table_3[,2]
colnames(SNVM_n)[-1]<-name_table_3[,2]
saveRDS(SNVH_n, file = "/Volumes/home/greally-lab/T_Trial/SNVHIGH_n.Rds")
saveRDS(SNVL_n, file = "/Volumes/home/greally-lab/T_Trial/SNVLOW_n.Rds")
saveRDS(SNVM_n, file = "/Volumes/home/greally-lab/T_Trial/SNVMOD_n.Rds")
SNVH_n[1:5,1:5]
name_table_3[1:5,]

class(SNVHIGH_table[,3])
as.numeric(SNVHIGH_table[,2])


apply(SNVHIGH_table,2,as.numeric)

sum(SNVHIGH_table[,3])
