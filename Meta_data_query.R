source("http://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")
library("TCGAbiolinks")

## get clinical info
query <- GDCquery(project = "TCGA-READ", 
                  data.category = "Clinical")
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")

write.csv(clinical,file="/Volumes/home/greally-lab/Taylor_Yu/READ_clinical_info.csv")


