#libraries
library(dplyr)
library(survival)
library(survminer)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(RTCGA)
library(RTCGA.clinical)
library(cgdsr)
library("annotate")
library("hgu133a.db") 

#
# Survival analysis with TCGA clinical data
#

# Making a GDC Query for TCGA-BRCA Clinical data
query_BRCA_clinical <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Clinical")
GDCdownload(query_BRCA_clinical)

#Assigning the clinical data to a data frame
BRCA_clinical <- GDCprepare_clinic(query_BRCA_clinical, clinical.info = "patient")


#Keeping only clinical entries that have "Positive" or "Negative" in the variables corrsponding to PR, ER and HER2 

trial<-BRCA_clinical[ BRCA_clinical$breast_carcinoma_progesterone_receptor_status == "Positive" 
                     |BRCA_clinical$breast_carcinoma_progesterone_receptor_status == "Negative" , ]

trialv2<-trial[trial$breast_carcinoma_estrogen_receptor_status == "Positive" |trial$breast_carcinoma_estrogen_receptor_status == "Negative",]
trialv3<-trialv2[trialv2$lab_proc_her2_neu_immunohistochemistry_receptor_status == "Positive" |trialv2$lab_proc_her2_neu_immunohistochemistry_receptor_status == "Negative",]

#Adding columns with regards to breast cancer type as well as a survival status column that is compatible with survivalTCGA
BRCA_clinical_final = within(trialv3, {
  luminal_A = ifelse((breast_carcinoma_estrogen_receptor_status == "Positive" | breast_carcinoma_progesterone_receptor_status== "Positive") 
                     & lab_proc_her2_neu_immunohistochemistry_receptor_status == "Negative"
                     , "True", "False")
  luminal_B = ifelse((breast_carcinoma_estrogen_receptor_status == "Positive" | breast_carcinoma_progesterone_receptor_status== "Positive") 
                       & lab_proc_her2_neu_immunohistochemistry_receptor_status == "Positive"
                       , "True", "False")
  basal = ifelse(breast_carcinoma_estrogen_receptor_status == "Negative" & breast_carcinoma_progesterone_receptor_status== "Negative"
                 & lab_proc_her2_neu_immunohistochemistry_receptor_status == "Negative"
                 , "True", "False")
  her2_plus = ifelse(breast_carcinoma_estrogen_receptor_status == "Negative" & breast_carcinoma_progesterone_receptor_status== "Negative"
                 & lab_proc_her2_neu_immunohistochemistry_receptor_status == "Positive"
                 , "True", "False")
  surv_status = ifelse(vital_status == "Alive"
                     , "alive", "dead")
})

#Generating a survival object including ER information to do survival analysis

clin_surv_ER <- survivalTCGA(BRCA_clinical_final, extract.cols = "breast_carcinoma_estrogen_receptor_status", extract.names = FALSE,
                  barcode.name = "bcr_patient_barcode",
                  event.name = "surv_status",
                  days.to.followup.name = "days_to_last_followup",
                  days.to.death.name = "days_to_death")

#Model the surivial object using ER and plot the Kaplan-Meier graph with log rank test p-value

sfit_ER <- survfit(Surv(times, surv_status)~breast_carcinoma_estrogen_receptor_status, data=clin_surv)
ggsurvplot(sfit_ER, pval=TRUE , legend.labs=c("Negative", "Positive" ), legend.title="Estrogen Receptor Status")

#Generating a survival object including PR information to do survival analysis

clin_surv_PR <- survivalTCGA(BRCA_clinical_final, extract.cols = "breast_carcinoma_progesterone_receptor_status", extract.names = FALSE,
                             barcode.name = "bcr_patient_barcode",
                             event.name = "surv_status",
                             days.to.followup.name = "days_to_last_followup",
                             days.to.death.name = "days_to_death")


#Model the surivial object using PR and plot the Kaplan-Meier graph with log rank test p-value

sfit_PR <- survfit(Surv(times, surv_status)~breast_carcinoma_progesterone_receptor_status, data=clin_surv_PR)
ggsurvplot(sfit_PR, pval=TRUE,  legend.labs=c("Negative", "Positive" ), legend.title="Progesteron Receptor Status")

#Generating a survival object including HER2 information to do survival analysis

clin_surv_HER <- survivalTCGA(BRCA_clinical_final, extract.cols = "lab_proc_her2_neu_immunohistochemistry_receptor_status", extract.names = FALSE,
                             barcode.name = "bcr_patient_barcode",
                             event.name = "surv_status",
                             days.to.followup.name = "days_to_last_followup",
                             days.to.death.name = "days_to_death")

#Model the surivial object using HER2 and plot the Kaplan-Meier graph with log rank test p-value

sfit_HER <- survfit(Surv(times, surv_status)~lab_proc_her2_neu_immunohistochemistry_receptor_status, data=clin_surv_HER)
ggsurvplot(sfit_HER, pval=TRUE, legend.labs=c("Negative", "Positive" ), legend.title="HER2 Receptor Status")

#Generating a survival object including basal information to do survival analysis (Not used in the project report)


clin_surv_3N <- survivalTCGA(BRCA_clinical_final, extract.cols = "basal", extract.names = FALSE,
                              barcode.name = "bcr_patient_barcode",
                              event.name = "surv_status",
                              days.to.followup.name = "days_to_last_followup",
                              days.to.death.name = "days_to_death")

#Model the surivial object using HER2 and plot the Kaplan-Meier graph with log rank test p-value   (Not used in the project report)

sfit_3N <- survfit(Surv(times, surv_status)~basal, data=clin_surv_3N)
ggsurvplot(sfit_3N, pval=TRUE, legend.labs=c("Negative", "Positive" ), legend.title="Basal")

#Generating a survival object including age information to do survival analysis (Used in the presentation but not used in the project report)
clin_surv_age <- survivalTCGA(BRCA_clinical_final, extract.cols = "age_at_initial_pathologic_diagnosis", extract.names = FALSE,
                              barcode.name = "bcr_patient_barcode",
                              event.name = "surv_status",
                              days.to.followup.name = "days_to_last_followup",
                              days.to.death.name = "days_to_death")

#Modeling using cox porportional hazard model and reporting the results Used in the presentation but not used in the project report)

surv_age.cox <- coxph(Surv(times, surv_status)~age_at_initial_pathologic_diagnosis, data=clin_surv_age)
summary(surv_age.cox)

#
#Replication of TCGA results using HT-Seq
#

#Making a GDC Query for TCGA-BRCA HT-seq expression data

query.exp.hg38 <- GDCquery(project = "TCGA-BRCA", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "HTSeq - FPKM-UQ")

#Downloading the data using the query and saving it as the file is very huge

GDCdownload(query.exp.hg38)
expdat <- GDCprepare(query = query.exp.hg38,
                     save = TRUE, 
                     save.filename = "exp.rda")


#Assigning the expression data 
BRCAMatrix <- assay(expdat,"HTSeq - FPKM-UQ")
#df_BRCA_Matrix <- as.data.frame(BRCAMatrix)


#names(df_BRCA_Matrix) <- substring(names(df_BRCA_Matrix),1,12)


#ENSG codes of IL6 CSF2 CCL5 VEGFA VEGFC, TCGA uses ENSG codes for genes so create a vector for interested genes
include_list <- c("ENSG00000136244","ENSG00000164400","ENSG00000271503","ENSG00000112715", "ENSG00000150630",
                  "ENSG00000134352", "ENSG00000198223", "ENSG00000100368", "ENSG00000160791","ENSG00000128052",
                  "ENSG00000099250","ENSG00000037280","ENSG00000118257")

#we search for interested genes by using the following vector
sub_BRCA_matrix<-BRCAMatrix[include_list,]

#removing the big matrix to save space
rm(BRCAMatrix)

#converting to a dataframe
df_BRCA_genes<-as.data.frame(sub_BRCA_matrix)

#renaming the row names
interested_genes <- c("IL6", "CSF2","CCl5","VEGFA","VEGFC", "GP130", "GMRA", "GMRB", "CCR5", "VEGFR2", "NRP1", "VEGFR3", "NRP2")
rownames(df_BRCA_genes) <- interested_genes


#Extracting patient id belonging to a type of breast cancer
basal_patient_list <-BRCA_clinical_final[BRCA_clinical_final$basal == "True",]$bcr_patient_barcode 
her2_plus_patient_list <-BRCA_clinical_final[BRCA_clinical_final$her2_plus == "True",]$bcr_patient_barcode 
luminal_A_patient_list <-BRCA_clinical_final[BRCA_clinical_final$luminal_A == "True",]$bcr_patient_barcode 
luminal_B_patient_list <-BRCA_clinical_final[BRCA_clinical_final$luminal_B == "True",]$bcr_patient_barcode 


patient.list<-list(basal_patient_list,her2_plus_patient_list,luminal_A_patient_list,luminal_B_patient_list)
patient.categories <- c("basal_patient_list","her2_plus_patient_list","luminal_A_patient_list","luminal_B_patient_list")

#Removing the tail portions of tail portions of names so that they can be matched
names(df_BRCA_genes) <- substring(names(df_BRCA_genes),1,12)

#Converting 0 to a small number such that taking the log does not result in negative infinity
df_BRCA_genes[df_BRCA_genes == 0] <- 0.000001
df_BRCA_genes.log= log2(df_BRCA_genes)

#t-test between expression profiles of subtypes
tcga_htseq_mat<-matrix(,ncol=12)
for (k in 1:13){
  allp<-c()
  for (i in 1:4){
    p_val<-c()
    log2_fold_change<-c()
    sub_type<-c()
    for (j in 1:4){
      if (i == j){
        next
      }
      x <- as.numeric(df_BRCA_genes.log[k,unlist(patient.list[i])[unlist(patient.list[i]) %in% names(df_BRCA_genes.log)]])
      y <- as.numeric(df_BRCA_genes.log[k,unlist(patient.list[j])[unlist(patient.list[j]) %in% names(df_BRCA_genes.log)]])
      p_val<-append(p_val,t.test(x, y, alternative="g")$p.value)
    
    }
    p_adjust <- p.adjust(p_val, "BH")
    allp <-append(allp, p_adjust)
  }
  tcga_htseq_mat<-rbind(tcga_htseq_mat,allp)
}
#The row names IL6  CSF2  CCL5  VEGFA VEGFC GP130(IL6ST)  GMRA(CSF2RA) GMRB(CSF2RB) CCR5  VEGFR2  NRP1  VEGFR3 NRP2
#The column names Basal vs HER2+	Basal vs Luminal A	Basal vs Luminal B	HER2+ vs Basal	HER2+ vs Luminal A	Her2+ vs Luminal B	Luminal A vs Basal	Luminal A vs HER2+	Luminal A vs Luminal B	Luminal B vs Basal	Luminal B vs HER2+	Luminal B vs Luminal A
#output the p-values to current working directory

write.csv(tcga_htseq_mat,file="TCGA_HTSeq_p_values.txt",row.names=FALSE)

#Making a CGDA object and obtaining expression profiles of interested genes from all patients

mycgds = CGDS("http://www.cbioportal.org/")
getGeneticProfiles(mycgds,'brca_tcga')
brca_expression_rna_seq<-t(getProfileData(mycgds,c("IL6", "CSF2","CCl5","VEGFA","VEGFC", "IL6ST", "CSF2RA", "CSF2RB", "CCR5", "VEGFR2", "NRP1", "VEGFR3", "NRP2"), 
                                          "brca_tcga_rna_seq_v2_mrna", "brca_tcga_all"))
df_BRCA_expression_rna_seq<-as.data.frame(brca_expression_rna_seq)

#modifying columns to allow matching
names(df_BRCA_expression_rna_seq) <- gsub(x = names(df_BRCA_expression_rna_seq),
                                          pattern = "\\.",
                                          replacement = "-")

#reorder rows and rename them
test_include<-c("IL6", "CSF2","CCL5","VEGFA","VEGFC", "IL6ST", "CSF2RA", "CSF2RB", "CCR5", "KDR", "NRP1", "FLT4", "NRP2")
rename_include<-c("IL6", "CSF2","CCL5","VEGFA","VEGFC", "IL6ST", "CSF2RA", "CSF2RB", "CCR5", "VEGFR2", "NRP1", "VEGFR3", "NRP2")
df_BRCA_expression_rna_seqv2<-df_BRCA_expression_rna_seq[test_include,]
rownames(df_BRCA_expression_rna_seqv2) <- rename_include

#trimming the suffix of column names to allow matching
names(df_BRCA_expression_rna_seqv2) <- substring(names(df_BRCA_expression_rna_seqv2),1,12)

#Converting 0 to a small number such that taking the log does not result in negative infinity
df_BRCA_expression_rna_seqv2[df_BRCA_expression_rna_seqv2 == 0] <- 0.000001
df_BRCA_expression_rna_seqv2.log= log2(df_BRCA_expression_rna_seqv2)

#t-test between expression profiles of subtypes
tcga_RNASeq_mat<-matrix(,ncol=12)
for (k in 1:13){
  allp<-c()
  for (i in 1:4){
    p_val<-c()
    log2_fold_change<-c()
    
    sub_type<-c()
    for (j in 1:4){
      if (i == j){
        next
      }
      x <- as.numeric(df_BRCA_expression_rna_seqv2.log[k,unlist(patient.list[i])[unlist(patient.list[i]) %in% names(df_BRCA_expression_rna_seqv2.log)]])
      
      y <- as.numeric(df_BRCA_expression_rna_seqv2.log[k,unlist(patient.list[j])[unlist(patient.list[j]) %in% names(df_BRCA_expression_rna_seqv2.log)]])
      p_val<-append(p_val,t.test(x, y, alternative="g")$p.value)
      
    }
    p_adjust <- p.adjust(p_val, "BH")
    allp <-append(allp, p_adjust)
  }
  tcga_RNASeq_mat<-rbind(tcga_RNASeq_mat,allp)
}

#The row names IL6  CSF2  CCL5  VEGFA VEGFC GP130(IL6ST)  GMRA(CSF2RA) GMRB(CSF2RB) CCR5  VEGFR2  NRP1  VEGFR3 NRP2
#The column names Basal vs HER2+	Basal vs Luminal A	Basal vs Luminal B	HER2+ vs Basal	HER2+ vs Luminal A	Her2+ vs Luminal B	Luminal A vs Basal	Luminal A vs HER2+	Luminal A vs Luminal B	Luminal B vs Basal	Luminal B vs HER2+	Luminal B vs Luminal A
#output the p-values to current working directory


write.csv(tcga_RNASeq_mat,file="TCGA_RNASeq_p_values.txt",row.names=FALSE)



#
#
#METABRIC REPLICATION
#
#


#I manually extracted the expression matrix of METABRIC with only the interested gene into data_expression_interested.txt 
#Read the file into R
metabric_expression_interested_genes<- read.table("data_expression_interested.txt", sep ="\t",  header =  TRUE, row.names = 1)
df_metabric_breast_expression_interested<-as.data.frame(metabric_expression_interested_genes)

#Reordering and renaming the rows
metabric_include_list<-c("IL6", "CSF2","CCL5","VEGFA","DQ896666", "IL6ST", "CSF2RA", "CSF2RB", "CCR5", "KDR", "NRP1", "FLT4", "NRP2")
df_metabric_breast_expression_interested<-df_metabric_breast_expression_interested[metabric_include_list,]
rename_include<-c("IL6", "CSF2","CCL5","VEGFA","VEGFC", "GP130", "CSF2RA", "CSF2RB", "CCR5", "VEGFR2", "NRP1", "VEGFR3", "NRP2")
rownames(df_metabric_breast_expression_interested) <- rename_include

#Converting column names to allow matching
names(df_metabric_breast_expression_interested) <- gsub(x = names(df_metabric_breast_expression_interested),
                                                        pattern = "\\.",
                                                        replacement = "-")

#Reading the case names of breast cancer patients studied in METABRIC study
cases<- as.data.frame(t(read.table("cases_nature_2012.txt", nrows = 1)))
cases$V1 <- as.character(cases$V1)

#Removing patient entries taht are not part of the METABRIC study
df_metabric_only_breast_expression_interested<-df_metabric_breast_expression_interested[,unlist(cases$V1)[unlist(cases$V1) 
                                                                                                          %in% names(df_metabric_breast_expression_interested)]]
#reading the clinical data of METABRIC
metabric_breast_clinical<- read.table("data_clinical_supp_patient.txt", sep ="\t",  header =  TRUE)
df_metabric_breast_clinical<-as.data.frame(metabric_breast_clinical)

#Adding columns with regards to breast cancer type
metabric_clinical_final = within(df_metabric_breast_clinical, {
  luminal_A = ifelse(CLAUDIN_SUBTYPE == "LumA"
                     , "Positive", "Negative")
  luminal_B = ifelse(CLAUDIN_SUBTYPE == "LumB"
                     , "Positive", "Negative")
  basal = ifelse(CLAUDIN_SUBTYPE == "Basal"
                 , "Positive", "Negative")
  her2_plus = ifelse(CLAUDIN_SUBTYPE == "Her2"
                     , "Positive", "Negative")
  
})

metabric_clinical_final$PATIENT_ID <- as.character(metabric_clinical_final$PATIENT_ID)

#Extracting patient id belonging to a type of breast cancer
metabric_basal_patient_list <-metabric_clinical_final[metabric_clinical_final$basal == "Positive",]$PATIENT_ID
metabric_her2_plus_patient_list <-metabric_clinical_final[metabric_clinical_final$her2_plus == "Positive",]$PATIENT_ID 
metabric_luminal_A_patient_list <-metabric_clinical_final[metabric_clinical_final$luminal_A == "Positive",]$PATIENT_ID
metabric_luminal_B_patient_list <-metabric_clinical_final[metabric_clinical_final$luminal_B == "Positive",]$PATIENT_ID



metabric_patient.list<-list(metabric_basal_patient_list,metabric_her2_plus_patient_list,metabric_luminal_A_patient_list,metabric_luminal_B_patient_list)
metabric_patient.categories <- c("metabric_basal_patient_list","metabric_her2_plus_patient_list","metabric_luminal_A_patient_list","metabric_luminal_B_patient_list")

#t-test between expression profiles of subtypes

metabric_mat<-matrix(,ncol=12)
for (k in 1:13){
  allp<-c()

  for (i in 1:4){
    
    p_val<-c()
    sub_type<-c()
    for (j in 1:4){
      if (i == j){
        next
      }
      x <- as.numeric(df_metabric_only_breast_expression_interested[k,unlist(metabric_patient.list[i])[unlist(metabric_patient.list[i]) 
                                                                                                       %in% names(df_metabric_only_breast_expression_interested)]])

      y <- as.numeric(df_metabric_only_breast_expression_interested[k,unlist(metabric_patient.list[j])[unlist(metabric_patient.list[j]) 
                                                                                                       %in% names(df_metabric_only_breast_expression_interested)]])

      p_val<-append(p_val,t.test(x, y, alternative="g")$p.value)

      
    }
    p_adjust <- p.adjust(p_val, "BH")
    allp <-append(allp, p_adjust)
  }
  metabric_mat<-rbind(metabric_mat,allp)
}

#The row names IL6  CSF2  CCL5  VEGFA VEGFC GP130(IL6ST)  GMRA(CSF2RA) GMRB(CSF2RB) CCR5  VEGFR2  NRP1  VEGFR3 NRP2
#The column names Basal vs HER2+	Basal vs Luminal A	Basal vs Luminal B	HER2+ vs Basal	HER2+ vs Luminal A	Her2+ vs Luminal B	Luminal A vs Basal	Luminal A vs HER2+	Luminal A vs Luminal B	Luminal B vs Basal	Luminal B vs HER2+	Luminal B vs Luminal A
#output the p-values to current working directory

write.csv(metabric_mat,file="METABRIC_p_values.txt",row.names=FALSE)

#
#Extending to GSE5847
#


#Manually extracted the expression matrix from GSE5847 data set and reading it

gse5847_expression_genes<- read.table("GSE5847.txt", sep ="\t",  header =  TRUE, row.names = 1)

#obtain a character vector of probe names

PROBES<- as.character(row.names(gse5847_expression_genes))

#map probe name to symbols

OUT <- select(hgu133a.db, PROBES, c("SYMBOL", "ENTREZID", "GENENAME"))
gse_rename_include<-OUT$SYMBOL


#manually extracted the expression matrix of interested gene using the mapping from above. The probe with the alphabetically largest value is chosen for genes that have
#multiple probes and reading it

gse5847_expression_interested_genes<- read.table("GSE5847_interested_genes.txt", sep ="\t",  header =  TRUE, row.names = 1)
df_gse5847_expression_interested_genes<-as.data.frame(gse5847_expression_interested_genes)

#renaming rows
rename_include<-c("IL6", "CSF2","CCL5","VEGFA","VEGFC", "GP130", "CSF2RA", "CSF2RB", "CCR5", "VEGFR2", "NRP1", "VEGFR3", "NRP2")
rownames(df_gse5847_expression_interested_genes) <- rename_include

#manually extracted the diagnosis information associated with GSE samples and reading it. Also replaced "Diagnosis: IBC" to "IBC" and "Diagnosis: non-IBC" to "non-IBC"
df_gse_sample_diagnosis<- as.data.frame(t(read.table("gse_sample_disagnosis.txt", sep ="\t")))
df_gse_sample_diagnosis$V1 <- as.character(df_gse_sample_diagnosis$V1)
df_gse_sample_diagnosis$V2 <- as.character(df_gse_sample_diagnosis$V2)

#Extracting patient id belonging to a type of breast cancer
gse_5847_IBC_patient_list <-df_gse_sample_diagnosis[df_gse_sample_diagnosis$V2 == "IBC",]$V1
gse_5847_non_IBC_patient_list <-df_gse_sample_diagnosis[df_gse_sample_diagnosis$V2 == "non-IBC",]$V1

gse_patient.list<-list(gse_5847_IBC_patient_list,gse_5847_non_IBC_patient_list)
gse_patient.categories <- c("gse_5847_IBC_patient_list","gse_5847_non_IBC_patient_list")

#t-test between expression profiles of subtypes

gse_mat<-matrix(,ncol=2)
for (k in 1:13){
  allp<-c()
  
  for (i in 1:2){
    
    p_val<-c()

    for (j in 1:2){
      if (i == j){
        next
      }
      x <- as.numeric(df_gse5847_expression_interested_genes[k,unlist(gse_patient.list[i])[unlist(gse_patient.list[i]) 
                                                                                           %in% names(df_gse5847_expression_interested_genes)]])

      y <- as.numeric(df_gse5847_expression_interested_genes[k,unlist(gse_patient.list[j])[unlist(gse_patient.list[j]) 
                                                                                           %in% names(df_gse5847_expression_interested_genes)]])

      p_val<-append(p_val,t.test(x, y, alternative="g")$p.value)

    }
    p_adjust <- p.adjust(p_val, "BH")

    #print (log2_fold_change)
    allp <-append(allp, p_adjust)
  }
  gse_mat<-rbind(gse_mat,allp)
}

#The row names IL6  CSF2  CCL5  VEGFA VEGFC GP130(IL6ST)  GMRA(CSF2RA) GMRB(CSF2RB) CCR5  VEGFR2  NRP1  VEGFR3 NRP2
#The column names IBC vs non-IBC non-IBC vs IBC
#output the p-values to current working directory

write.csv(gse_mat,file="GSE_p_values.txt",row.names=FALSE)
