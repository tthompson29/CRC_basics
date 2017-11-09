META_data = read.table("/Volumes/home/greally-lab/T_Trial/TCGA_data_CRC/CRC_clinical.txt")
familial_tumorPT = as.data.frame(matrix(nrow = 0,ncol = 76))
colnames(familial_tumorPT)= colnames(META_data)

for (n in 1:nrow(META_data)){
  y = META_data[n,]
  if (y$microsatellite_instability == "YES"){
    familial_tumorPT = rbind(familial_tumorPT,y)
  }
  if (y$loss_expression_of_mismatch_repair_proteins_by_ihc == "YES"){
    familial_tumorPT = rbind(familial_tumorPT,y)
  }
  if (y$synchronous_colon_cancer_present =="YES"){
    familial_tumorPT = rbind(familial_tumorPT,y)
  }
}
dim(familial_tumorPT) # 178 76
##### Check for duplicates
exclusion_pts = unique(sort(familial_tumorPT$bcr_patient_barcode))
length(exclusion_pts) # 89 unique pts that have familial tumors

########exclusion_pts = gsub("TCGA-","", familal_tumorPT$bcr_patient_barcode)

View(exclusion_pts)
write.table(exclusion_pts,file = "/Volumes/home/greally-lab/T_Trial/Familail_tumor_pts.txt", sep = "\t")
