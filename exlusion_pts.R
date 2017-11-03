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
    familal_tumorPT = rbind(familial_tumorPT,y)
  }
}
dim(familal_tumorPT)
exclusion_pts = familial_tumorPT$bcr_patient_barcode
exclusion_pts = gsub("TCGA-","", exclusion_pts)
