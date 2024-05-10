##### TCGA data processing
options(stringsAsFactors = F)
if (!requireNamespace("readr", quietly = TRUE))
  install.packages("readr")
library(readr)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(limma)
library(biomaRt)
library(tibble)
library(igraph)
library(readr)
#####
library(SummarizedExperiment)
library(dplyr)
library(lubridate)
library(tidyverse)
library(EDASeq)
############
library(Rtsne)
library(ggplot2)
install.packages("plotly")
library(plotly)
#####################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("EDASeq", quietly = TRUE))
  BiocManager::install("EDASeq")



#####################################################workround###############################################################
#https://github.com/BioinformaticsFMRP/TCGAbiolinks/issues/627
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinksGUI.data")
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")

#############################################################################################################################
#### Collect gene expression data#################################################################
setwd("C:/Users/gklizh/Documents/Workspace/code_and_data19/")
dataset_TCGA <- c("BLCA", "BRCA", "COAD", "ESCA", "KICH", "KIRC", 
                  "KIRP", "LIHC", "LUAD", "LUSC", "PAAD", "PRAD", 
                  "READ", "SKCM", "STAD", "THCA", "THYM", "UCEC")

tried_groups<-c("BLCA",  "COAD",  "KICH", "KIRC", "KIRP","SKCM","LIHC","LUAD","LUSC","UCEC","READ","PAAD","BRCA")
still_groups<-c( )
not_installed<-c("ESCA","THYM","THCA","STAD")

#print(cancer_type)
in_pool<-"THYM"
doing_groups<-c("THYM")
##############################################################################
cancer_type<-"SKCM"
query_exp <- GDCquery(project = paste0("TCGA-",cancer_type), 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")#
GDCdownload(query_exp) 
data_exp <- GDCprepare(query_exp) 
print(cancer_type)
##########################################################################################
#https://www.biostars.org/p/9544856/
fpkm_uq_unstrand<-assays(data_exp)[["fpkm_uq_unstrand"]]
exp_assay<-as.data.frame(fpkm_uq_unstrand)
exp_rowRanges <- as.data.frame(rowRanges(data_exp)) # Gene annotation
saveRDS(exp_assay,paste0("./data/tcga_data/",cancer_type,"_exp_assay.RData"))
exp_assay<-readRDS(paste0("./data/tcga_data/",cancer_type,"_exp_assay.RData"))
saveRDS(exp_rowRanges, paste0("./data/tcga_data/",cancer_type,"_exp_rowRangess.RData"))
#exp_rowRanges <- as.data.frame(rowRanges(exp_assay))
############################################################################################################
#### Collect DNA methylation data
query_mty <- GDCquery(project = paste0("TCGA-",cancer_type), 
                      data.category = "DNA Methylation",
                      platform = "Illumina Human Methylation 450",
                      data.type = "Methylation Beta Value")
GDCdownload(query_mty)
data_mty <- GDCprepare(query_mty)
mty_assay <- as.data.frame(assay(data_mty)) # DNA methylation matrix
mty_colData <- as.data.frame(colData(data_mty)) # Patient annotation(475 patients)
#View(mty_colData)
mty_rowRanges <- as.data.frame(rowRanges(data_mty)) # cg probe annotation
#######################Save files##########################################
saveRDS(mty_colData, paste0("./data/tcga_data/",cancer_type,"_mty_colData.RData"))
saveRDS(mty_assay, paste0("./data/tcga_data/",cancer_type,"_mty_assay.RData"))
saveRDS(mty_rowRanges, paste0("./data/tcga_data/",cancer_type,"_mty_rowRanges.RData"))

mty_colData<-readRDS(paste0("./data/tcga_data/",cancer_type,"_mty_colData.RData"))
mty_assay<-readRDS(paste0("./data/tcga_data/",cancer_type,"_mty_assay.RData"))
mty_rowRanges<-readRDS(paste0("./data/tcga_data/",cancer_type,"_mty_rowRanges.RData"))
#########################################Masked Somatic Mutation ################################################
query <- GDCquery(
  project = paste0("TCGA-", cancer_type),
  data.category = "Simple Nucleotide Variation", 
  access = "open",
  data.type = "Masked Somatic Mutation", 
  #workflow.type = "MuSE Variant Aggregation and Masking"
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query)
print(cancer_type)
######################################################################################
# file_paths <- getResults(query)$file_id
# # #file_path<-"0b0de7da-a900-4595-bf65-805513874ecf"
# # library(readr)
# # 
#  data <- read_tsv("./GDCdata/TCGA-ESCA/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/1cdbcfc2-2cfc-435b-b1f3-ea1824983b1b/0695bd65-8748-4396-b47b-f3fa238895fb.wxs.aliquot_ensemble_masked.maf.gz", comment = "#")
# # 
#  View(data)
# # data$Tumor_Seq_Allele2 <- as.character(data$Tumor_Seq_Allele2)
# # #View(data)
# # #View(data)
# # 
# gz_files <- c()
# for (file_path in file_paths) {
# #   # 使用 list.files 函数找到所有 .gz 文件
#   gz_files <- c(gz_files, list.files(path = paste0("./GDCdata/TCGA-ESCA/Simple_Nucleotide_Variation/Masked_Somatic_Mutation/",file_path), pattern = "\\.gz$", full.names = TRUE))
#  }
# # 
# print(gz_files)
# # 
# # 
# # for (gz_file in gz_files) {
# #   df <- read_tsv(gz_file,comment = "#",show_col_types = FALSE)
# #   #df <- read.table(file_path, header = TRUE, stringsAsFactors = FALSE,comment = "#"))
# #   df$Tumor_Seq_Allele2 <- as.character(df$Tumor_Seq_Allele2)
# #   #write.table(df, file_path, row.names = FALSE)
# #   temp_file <- tempfile(fileext = ".gz")
# #   write.table(df, gzfile(temp_file), row.names = FALSE, sep = "\t", quote = FALSE)
# #   #gzip(temp_file, destname = "dbbd21e0-3f6a-43cd-868e-5d3ed3e7a024.wxs.aliquot_ensemble_masked.maf.gz", overwrite = TRUE)
# #   file.copy(temp_file, gz_file, overwrite = TRUE)
# # }
###################################################################################
#View(as.data.frame(file_paths))
maf <- GDCprepare(query)
#View(maf)
#View(clinicalInfo)
filename1 <- paste0("./data/tcga_data/", cancer_type, "_maf_test06.RData")
#print(filename1)
saveRDS(maf, filename1)
maf<-readRDS(filename1)
#############################################################"Copy Number Variation"###################################################################################
print(cancer_type)
query_cnv <- GDCquery(project = paste0("TCGA-", cancer_type),
                      data.category = "Copy Number Variation",
                      data.type = "Gene Level Copy Number")

GDCdownload(query_cnv)
potential_df <- query_cnv$results[[1]]
potential_df$created_datetime <- ymd_hms(potential_df$created_datetime)
potential_df$timestamp <- as.numeric(potential_df$created_datetime)

#######################################################################
unique_cases_df <- potential_df %>%
  distinct(cases, .keep_all = TRUE) %>%
  group_by(cases) %>%
  mutate(datetime = ymd_hms(created_datetime)) %>%  # 将created_datetime转换为日期时间对象
  filter(datetime == max(datetime)) %>%            # 过滤出最新的日期时间
  ungroup() %>%
  dplyr::select(-datetime)%>%
  dplyr::select(-timestamp)
#####################################################
query_cnv$results[[1]] <- unique_cases_df
##########################################################################
saveRDS(query_cnv,paste0("./data/tcga_data/",cancer_type,"_query_cnv_test06.RData"))
query_cnv<-readRDS(paste0("./data/tcga_data/",cancer_type,"_query_cnv_test06.RData"))
write.csv(query_cnv$results[[1]], file = paste0("./data/tcga_data/",cancer_type,"_query_results_test06.csv", row.names = TRUE))
############################################################
query_cnv_copy<-query_cnv
query_cnv_copy$results[[1]] <- query_cnv_copy$results[[1]] %>%
  distinct(sample.submitter_id, .keep_all = TRUE)

GDCdownload(query_cnv_copy)
data_cnv <- GDCprepare(query_cnv_copy)
data_cnv_savage<-data_cnv
saveRDS(data_cnv_savage,"./data/tcga_data/data_cnv_savage_GDCprepare_test.RData")
data_temp<-data_cnv_savage#
View(data_temp)
data_temp <- as.data.frame(assay(data_temp))
#View(data_temp)
dim(data_temp)
saveRDS(data_temp,paste0("./data/tcga_data/",cancer_type,"_data_temp_savage_test.csv"))
data_temp<-readRDS(paste0("./data/tcga_data/",cancer_type,"_data_temp_savage_test.csv"))
#################################################### Collect clinical data################################################################################################
clinical <- GDCquery_clinic(project = paste0("TCGA-",cancer_type), type = "clinical")
#View(clinical)
####################################################Collect clinical radiation and drug therapy data###############################################
#print(cancer_type)
query <- GDCquery(project = paste0("TCGA-",cancer_type), 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab"
)
GDCdownload(query)
sckm.tab.all <- GDCprepare(query)
column_name0 <- paste0("clinical_drug_", tolower(cancer_type))
therapy <- sckm.tab.all[[column_name0]] # clinical drug therapy ->NULL
#View(therapy)
saveRDS(therapy, paste0("./data/tcga_data/",cancer_type,"_therapy_test.RData"))
therapy<-readRDS(paste0("./data/tcga_data/",cancer_type,"_therapy_test.RData"))
therapy$pharmaceutical_therapy_type # Therapy types
column_name1<-paste0("clinical_radiation_", tolower(cancer_type))
radiation <- sckm.tab.all[[column_name1]] # Clinical radiation therapy
#View(radiation)
saveRDS(radiation, paste0("./data/tcga_data/",cancer_type,"_radiation_test.RData"))
radiation<-readRDS(paste0("./data/tcga_data/",cancer_type,"_radiation_test.RData"))
##########################################################################################################################

########################################################### Process gene expression data#################################

exp_assay_frame<-exp_assay
dim(exp_assay_frame)
colnames(exp_assay_frame) <- substr(colnames(exp_assay_frame), 1, 16)
true_index<-!(is.na(exp_rowRanges[row.names(exp_assay_frame), "gene_name"]))
exp_assay_frame$SYMBOL <- exp_rowRanges[row.names(exp_assay_frame), "gene_name"]
exp_assay_frame<-exp_assay_frame[true_index,]
#which(is.na(exp_assay_frame$SYMBOL))
saveRDS(exp_assay_frame,paste0("./data/tcga_data/",cancer_type,"_exp_assay_frame01.RData"))
####################################################################################################################################################
num_rows <- nrow(exp_assay_frame)
num_cols <- ncol(exp_assay_frame)

exp_assay_test1 <- exp_assay_frame[,c(num_cols,1:num_cols-1)]

exp_assay_test2 <- as.matrix(exp_assay_test1)

row.names(exp_assay_test2) <- exp_assay_test2[,1]

exp_assay_test2 <- exp_assay_test2[,2:num_cols] # Convert Ensemble ID to corresponding gene symbols

row_names_tcga<-rownames(as.data.frame(exp_assay_test2)) 

exp_assay_test2 <- matrix(as.numeric(exp_assay_test2), nrow = nrow(exp_assay_test2), dimnames = list(row.names(exp_assay_test2), colnames(exp_assay_test2)))
any(is.na(row.names(exp_assay_test2)))
exp_assay_test2 <- avereps(exp_assay_test2) # If the gene corresponds to multiple gene expression values, take the average value,remove duplicated gene name

exp_assy_trial<-exp_assay_test2
#############################################################################################################################

exp_intgr <- t(exp_assy_trial) # Gene expression data matrix
exp_intgr <- log10(exp_intgr + 1) # 1og10 conversion
saveRDS(exp_intgr, paste0("./data/tcga_data_processed/",cancer_type,"_exp_intgr_all01.RData"))
write.csv(exp_intgr,paste0("./data/tcga_data_processed/",cancer_type,"_exp_intgr_all01.csv"))
##############################################################################################################################
#### Process DNA methylation data

data_mty <- subset(data_mty, subset = (rowSums(is.na(assay(data_mty)))==0)) 

data_mty<- subset(data_mty, subset =!is.na(as.data.frame(rowRanges(data_mty))$gene_HGNC))
dim(data_mty)

#############################################################

mty_assay <- as.data.frame(assay(data_mty))
mty_rowRanges <- as.data.frame(rowRanges(data_mty))
mty_symbol <- strsplit(mty_rowRanges$gene_HGNC, split = ";")
names(mty_symbol) <- row.names(mty_rowRanges)
mty_symbol <- lapply(mty_symbol, FUN = function(x){x<-unique(x)})
mty_symbol <- as.matrix(unlist(mty_symbol))
row.names(mty_symbol) <- substr(row.names(mty_symbol), 1, 10)
mty_symbol <- data.frame("probe" = row.names(mty_symbol),"SYMBOL" = mty_symbol[,1])#cg14528386.1 ->H19_1
mty_assay$probe <- row.names(mty_assay)
mty_assay <- merge(mty_assay, mty_symbol, by.x = "probe", by.y = "probe") 
num_rows <- nrow(mty_assay)

# 获取列数
num_cols <- ncol(mty_assay)
mty_mat <- as.matrix(mty_assay[,2:num_cols])
colnames(mty_mat) <- substr(colnames(mty_mat), 1, 16)
row.names(mty_mat) <- mty_assay$SYMBOL # Convert probe ID to corresponding gene symbols
mty_mat <- matrix(as.numeric(mty_mat), nrow = nrow(mty_mat), dimnames = list(row.names(mty_mat), colnames(mty_mat)))
mty_mat <- avereps(mty_mat) # If the gene corresponds to multiple methylation values, take the average value
dim(mty_mat)
write.csv(mty_mat,paste0("./data/",cancer_type,"_mty_mat_all.csv"))
saveRDS(mty_mat,paste0("./data/",cancer_type,"_mty_mat_all.RData"))
mty_mat_trial<-readRDS(paste0("./data/",cancer_type,"_mty_mat_all.RData"))
mty_mat_trial <- t(mty_mat_trial) # Gene expression data matrix
saveRDS(mty_mat_trial, paste0("./data/tcga_data_processed/",cancer_type,"_mty_mat_all01.RData"))
write.csv(mty_mat_trial,paste0("./data/tcga_data_processed/",cancer_type,"_mty_mat_all01.csv"))
samples <- intersect(colnames(exp_intgr), colnames(mty_mat_trial)) 
print(samples)
#dim(samples)
saveRDS(samples, paste0("./data/tcga_data_processed/",cancer_type,"_samples.RData"))

#### Integrate the genomic data into the network
exp_intgr <- exp_intgr[,samples]
mty_mat_trial <- mty_mat_trial[,samples]
rownamesSamples<-intersect(rownames(exp_intgr),rownames(mty_mat_trial))
length(rownamesSamples)
#class(exp_intgr)# "matrix" "array
exp_intgr_trial<-exp_intgr[rownamesSamples,]
mty_mat_trial<-mty_mat_trial[rownamesSamples,]
dim(exp_intgr_trial)
#View(exp_intgr_trial)
################################################################################################################################################################

clinicalInfo <- mty_colData
clinicalInfo<-clinicalInfo[!duplicated(clinicalInfo$sample),]
indices <- which(clinicalInfo$shortLetterCode == "NT",clinicalInfo$sample)
#print(indices)
#View(clinicalInfo)

NT_remove_name<-(rownames(clinicalInfo[indices,]))
print(NT_remove_name)
# View(NT_remove_name)
# new_clinicalInfo<-clinicalInfo[samples,]
# length(samples)
# #dim(new_clinicalInfo)
# new_clinicalInfo<-clinicalInfo[-indices,]
# View(new_clinicalInfo)
# saveRDS(new_clinicalInfo, paste0("./data/tcga_data_processed/",cancer_type,"_clinical_info_test06.RData"))
############################################################
new_exp_intgr<-exp_intgr_trial[!(rownames(exp_intgr_trial) %in% NT_remove_name),]

new_mty_intgr<-mty_mat_trial[!(rownames(mty_mat_trial) %in% NT_remove_name),]

saveRDS(new_exp_intgr, paste0("./data/tcga_data_processed/",cancer_type,"_exp_intgr_all.RData"))
write.csv(new_exp_intgr,paste0("./data/tcga_data_processed/",cancer_type,"_exp_intgr_all.csv"))

write.csv(new_mty_intgr,paste0("./data/tcga_data_processed/",cancer_type,"_mty_mat_all_test.csv"))
saveRDS(new_mty_intgr,paste0("./data/tcga_data_processed/",cancer_type,"_mty_mat_all_test.RData"))

#dim(new_exp_intgr)
#dim(new_mty_intgr)

########################################################MAF start for data_cnv######################################################################

maf <- maf[,c("Tumor_Sample_Barcode","Hugo_Symbol","Gene","Variant_Classification")]

rnames <- unique(maf$Hugo_Symbol)
cnames <- unique(maf$Tumor_Sample_Barcode)
snv_count <- matrix(data = 0, nrow = length(rnames), ncol = length(cnames), dimnames = list(rnames,cnames)) 

# Calculate the frequency of genes' variants
for(i in 1:nrow(maf)){
  rname <- maf[i,]$Hugo_Symbol
  cname <- maf[i,]$Tumor_Sample_Barcode
  snv_count[rname, cname] <- snv_count[rname,cname] + 1
}
dim(snv_count)
#https://www.notion.so/57782b3efe144e91a93ae602f85a1cbd

#View(snv_count)
colnames(snv_count) <- substr(colnames(snv_count), 1, 16)

###############################################################################################!!!!!!!!!!!!!!!

data_temp$ensembl<-rownames(data_temp)
#View(as.data.frame(rownames(data_temp)))

rownames(data_temp)<-seq_along(rownames(data_temp))

data_temp$ensembl<-substr((data_temp$ensembl),1, 15)

na_count <- function(x) sum(is.na(x))

data_temp$na_count <- apply(data_temp, 1, na_count)

data_temp <- data_temp[order(data_temp$ensembl, data_temp$na_count), ]

data_temp <- data_temp[!duplicated(data_temp$ensembl), ]
#remove na_count
data_temp$na_count <- NULL

row.names(data_temp)<-data_temp[,ncol(data_temp)]

#
#dim(data_temp)
data_cnv_tmp<-data_temp
#dim(data_cnv)
colnames(data_cnv_tmp) <- substr(colnames(data_cnv_tmp), 1, 16)

#View(data_cnv_tmp)

###############################################################################################!!!!!!!!!!!!!!!!!!!!
col_data_cnv<-ncol(data_cnv_tmp)
#View(as.data.frame(col_data_cnv))#146
##########################################################################################################
# Convert Ensemble ID to corresponding gene symbols
# Using biomaRt for gene ID conversion will cause some loss
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")


#class(ensembl)

#View(ensembl)
#View(listFilters(ensembl))

#	hgnc_symbol  HGNC symbol(s) [e.g. A1BG]
#ensembl_gene_id  Gene stable ID(s) [e.g. ENSG00000000003]

#print(row.names(data_cnv))
#View(data_cnv)

#View(as.data.frame(row.names(data_cnv_tmp)))


cnv_df <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                filters = c("ensembl_gene_id"),
                values = row.names(data_cnv_tmp),
                mart = ensembl)

#View(cnv_df)
##############################################################
cnv_df <- cnv_df[which(cnv_df$hgnc_symbol !=''),]
#View(cnv_df)

cnv_df<-unique(cnv_df)

colnames(data_cnv_tmp)[colnames(data_cnv_tmp)=="ensembl"]<-"emsembl"

unique_cols<-!duplicated(t(data_cnv_tmp))
data_cnv_tmp<-data_cnv_tmp[,unique_cols]
#############################################################################################################

data_cnv_tmp <- merge(x = data_cnv_tmp, y = cnv_df, by.x = "emsembl", by.y = "ensembl_gene_id")
data_cnv_tmp <- as.matrix(data_cnv_tmp)
dim(data_cnv_tmp)
#############################################################################################################
##########################################################################################################

row.names(data_cnv_tmp) <- data_cnv_tmp[,ncol(data_cnv_tmp)]
#View(data_cnv_tmp)
data_cnv_tmp <- data_cnv_tmp[,2:ncol(data_cnv_tmp)-2]  # Convert Ensemble ID to corresponding gene symbols
data_cnv_tmp <- matrix(as.numeric(data_cnv_tmp), nrow = nrow(data_cnv_tmp), dimnames = list(row.names(data_cnv_tmp), colnames(data_cnv_tmp)))
data_cnv_tmp <- data_cnv_tmp[!duplicated(row.names(data_cnv_tmp)),] # Only one gene PRAMEF7 has duplicate copy number variation value, and the duplicate value is the same
data_cnv_tmp<-data_cnv_tmp[,-1]

print(rownamesSamples)
samples_f1<-intersect(rownamesSamples, colnames(snv_count))
print(samples_f1)

print(colnames(snv_count))
samples_f2 <- intersect(samples_f1, colnames(data_cnv_tmp))
print(samples_f2)
print(cancer_type)
saveRDS(samples_f2,paste0("./data/tcga_data_processed/",cancer_type,"_samples.RData"))

sample_f2_check<-readRDS(paste0("./data/tcga_data_processed/",cancer_type,"_samples.RData"))
print(sample_f2_check)
#length(samples_f2)

# print(samples_f2)
# 
# class(new_exp_intgr)
# class(samples_f2)
# new_exp_intgr_tmp1 <- new_exp_intgr[,samples_f2]
# new_mty_intgr_tmp1 <- new_mty_intgr[,samples_f2]
# snv_count_tmp1 <- snv_count[,samples_f2]
# data_cnv_tmp1 <- data_cnv_tmp[,samples_f2]
# 
# print(samples_f2[2])

# dim(snv_count)
# length(samples)
# print(samples[ZZZ3])
# View(samples)
# dim(data_cnv_tmp)
# View(as.data.frame(samples))
# View(data_cnv_tmp)
# max(data_cnv_tmp)
#class(samples)
#rows_to_keep <- row.names(data_cnv_tmp) %in% samples
#class(rows_to_keep)
#data_cnv_tmp <- data_cnv_tmp[rows_to_keep, ]
#cols_to_keep<-colnames(data_cnv_tmp) %in% rownamesSamples


#class(cols_to_keep)
#data_cnv_samfilter<-data_cnv_tmp[,cols_to_keep]
#data_cnv_samfilter<-data_cnv_tmp[samples,]
#dim(data_cnv_tmp)
#max(samples)
#length(samples)
#dim(data_cnv_samfilter)



#View(as.data.frame(samples))
#View(data_cnv_tmp)
#snv_count_samfilter<-snv_count[samples,]
#snv_count_samfilter<-snv_count[,cols_to_keep] 
#dim(snv_count)
#data_cnv_samfilter<-t(data_cnv_samfilter)
#snv_count_samfilter<-t(snv_count_samfilter)

#### Keep patient samples with four genomic profiles

#View(data_cnv)

#data_cnv_samfilter<- data_cnv[,samples]

#dim(data_cnv_samfilter)


print(samples_f2)

print(colnames(new_exp_intgr))

expSample_to_keep<-rownames(new_exp_intgr) %in% samples_f2
mytSample_to_keep<-rownames(new_mty_intgr) %in% samples_f2
snvSample_to_keep<-colnames(snv_count) %in% samples_f2 
dncSample_to_keep<-colnames(data_cnv_tmp) %in% samples_f2
print(samples_f2)
#View(snv_count)
#View(data_cnv_tmp)

#View(new_exp_intgr)
# length(expSample_to_keep)
# print(expSample_to_keep)
# print(dncSample_to_keep)
#print(expSample_to_keep)
###########################################
new_exp_intgr_tmp01<-new_exp_intgr[expSample_to_keep,]

new_mty_intgr_tmp01<-new_mty_intgr[mytSample_to_keep,]
# dim(new_mty_intgr_tmp01)
###########################################

snv_count_tmp01<-snv_count[,snvSample_to_keep]
data_cnv_tmp01<-data_cnv_tmp[,dncSample_to_keep]

# dim(snv_count_tmp01)
# dim(data_cnv_tmp01)

snv_count_tmp01<-t(snv_count_tmp01)
data_cnv_tmp01<-t(data_cnv_tmp01)


#View(NT_remove_name)
#class(clinicalInfo) data.frame
new_clinicalInfo<-clinicalInfo[samples_f2,]
View(as.data.frame(samples_f2))
#View(new_clinicalInfo)
#length(samples_f2)
#dim(new_clinicalInfo)
#print(indices)
#length(indices)

if (length(indices) > 0) {
  new_clinicalInfo <- clinicalInfo[-indices,]
} else {
  print("The length of NT indices is 0")  # 或者进行其他操作
}

#View(indices)
#View(new_clinicalInfo)
saveRDS(new_clinicalInfo, paste0("./data/tcga_data_processed/",cancer_type,"_clinical_info_test01.RData"))

View(new_clinicalInfo)

#View(snv_count_tmp01)
###########################################


#saveRDS(snv_count, paste0("./data/tcga_data_processed/",cancer_type,"_snv_intgr_test06.RData"))
#saveRDS(data_cnv_samfilter, paste0("./data/tcga_data_processed/",cancer_type,"_cnv_intgr_test06.RData"))
####################################################################added extra information####################
# graph_comp<-readRDS("./data/network_processed/graph_comp.RData")
# graph_comp_updated <- upgrade_graph(graph_comp)
# 
# genes <- as.character(unlist(vertex.attributes(graph_comp_updated))) # Nodes in the maximum connected subgraph
# #########################################################################################################
# y1 <- which(row.names(exp_assay) %in% genes)
# #https://www.notion.so/R-in-16bd93aeb8f940fe9a132f5d9d7d3ec8
# #View(as.data.frame(exp_assay))
# #class(y1)
# #head(y1)
# #length(y1)
# exp_intgr <- exp_assay[y1,] 
# exp_intgr <- t(exp_intgr) # Gene expression data matrix
# exp_intgr <- log10(exp_intgr + 1) # 1og10 conversion
# saveRDS(exp_intgr, paste0("./data/tcga_data_processed/",cancer_type,"_exp_intgr.RData"))
# write.csv(exp_intgr,paste0("./data/tcga_data_processed/",cancer_type,"_exp_intgr.csv"))
# 
# 
# 
# 
# 
# y2 <- which(row.names(mty_mat) %in% genes)
# mty_intgr <- mty_mat[y2,]
# mty_intgr <- t(mty_intgr) # DNA methylation data matrix
# saveRDS(mty_intgr, paste0("./data/tcga_data_processed/",cancer_type,"_mty_intgr.RData"))
# write.csv(mty_intgr,paste0("./data/tcga_data_processed/",cancer_type,"_mty_intgr.csv"))
# 
# y3 <- which(row.names(snv_count) %in% genes)
# snv_intgr <- snv_count[y3,]
# snv_intgr <- t(snv_intgr) # Gene mutation data matrix
# saveRDS(snv_intgr, paste0("./data/tcga_data_processed/",cancer_type,"_snv_intgr.RData"))
# write.csv(snv_intgr,paste0("./data/tcga_data_processed/",cancer_type,"_snv_intgr.CSV"))
# 
# saveRDS(snv_intgr, paste0("./data/tcga_data_processed/",cancer_type,"_snv_intgr.RData"))
# 
# 
# #write.csv(snv_count,paste0("./data/tcga_data_processed/",cancer_type,"_snv_intgr_test06.CSV"))
# 
# 
# 
# y4 <- which(row.names(data_cnv) %in% genes)
# cnv_intgr <- data_cnv[y4,]
# cnv_intgr <- t(cnv_intgr) # Copy number variation data matrix
# saveRDS(cnv_intgr, paste0("./data/tcga_data_processed/",cancer_type,"_cnv_intgr.RData"))
# #write.csv(clinicalInfo,paste0("./data/tcga_data_processed/",cancer_type,"_cnv_intgr.CSV"))
# saveRDS(data_cnv, paste0("./data/tcga_data_processed/",cancer_type,"_cnv_intgr_test06.RData"))
# 
# 
# #write.csv(clinicalInfo,paste0("./data/tcga_data_processed/",cancer_type,"_cnv_intgr.CSV")
# 
# 
# #clinicalInfo <- mty_colData
# #View(mty_colData)
# #row.names(clinicalInfo) <- clinicalInfo$sample
# #clinicalInfo <- clinicalInfo[samples,]
# #saveRDS(clinicalInfo, paste0("./data/tcga_data_processed/",cancer_type,"_clinical_info_test06.RData"))
# 
# 
# #write.csv(clinicalInfo,paste0("./data/tcga_data_processed/",cancer_type,"_clinical_info.csv"))
# # clinicalInfo <- mty_colData
# # row.names(clinicalInfo) <- clinicalInfo$sample
# # clinicalInfo <- clinicalInfo[load_samples,]
# # saveRDS(clinicalInfo, paste0("./data/tcga_data_processed/",cancer_type,"_clinical_info_test06.RData"))
# # print(paste0("./data/tcga_data_processed/",cancer_type,"_clinical_info_test06.RData"))
# # #write.csv(clinicalInfo,paste0("./data/tcga_data_processed/",cancer_type,"_clinical_info.csv"))
# setwd("C:/Users/gklizh/Documents/Workspace/code_and_data19/")
# save.image(file = 'my_work02_20240419.RData')
# #load('my_work02_20240419.RData')

######################################section2 singlass for communities##########################################################

##### The melanoma network and network community detection
options(stringsAsFactors = F)
library(igraph)
library(visNetwork)

library(OmnipathR)#Get interactions data
library(igraph)#network analysis
library(ggraph)#network visualization
library(RColorBrewer)#Network visualization, heat map drawing

library(factoextra)#principal component analysis
library(FactoMineR)#principal component analysis

library(DESeq2)#differential expression analysis
library(SANTA)#Node Score Calculation

library(ComplexHeatmap)#heatmap drawing
library(stringr)#regular expression

library(randomForest)#Random Forest
library(caret)#Cross-validation
library(pROC)#AUC

library(clusterProfiler)#enrichment analysis
library(org.Hs.eg.db)
library(enrichplot)#enrichment analysis
library(msigdbr)#enrichment analysis
library(dplyr)

library(survival)#survival analysis
library(survminer)#survival analysis
library(glmnet)# cox regression
library(timeROC)#AUC

library(reshape2)#ggplot2 drawing

library(coin)#permutation test

library(networkD3)#sankey diagram

library(janitor)# clean data
library(here)# home directory setting
library(ggplot2) # draw figures
library(tidyr) # data processing
if (!require("visNetwork")) install.packages("visNetwork")
if (!require("igraph")) install.packages("igraph")
library(visNetwork)
library(igraph)
### The melanoma network
#cancer_type<-"BLCA"
print(cancer_type)
setwd("C:/Users/gklizh/Documents/Workspace/fan_project/code_data")
edgetest<-read.csv(paste0("Raw data/network_merge/edge/", cancer_type, ".csv"))
vertextest<-read.csv(paste0("Raw data/network_merge/vertex/", cancer_type, ".csv"))
interactions_E <- edgetest
interactions_V <- vertextest
graph_comp <- graph_from_data_frame(interactions_E, directed = TRUE, vertices = interactions_V)
graph_comp <- decompose(graph_comp, min.vertices = 15)[[1]]
graph_comp <- set_vertex_attr(graph_comp, "degree", index = V(graph_comp)$name, degree(graph_comp))

print(class(degree(graph_comp)))
print(length(degree(graph_comp)))
print(class(V(graph_comp)$log2FoldChange))
print(length(V(graph_comp)$log2FoldChange))
#print(graph_comp)
#################################################################

setwd("C:/Users/gklizh/Documents/Workspace/code_and_data19/")
spg01 <- list() # List of spinglass simulations
spg_mod01 <- numeric() # List of modularity simulations
count01<-1
for (k in 1:50){
  print(Sys.time())
  print(k)
  spg01[[k]] <- cluster_spinglass(graph_comp,
                                  spins = 25,
                                  weights = E(graph_comp)$weight,
                                  implementation = "neg")
  
  spg_mod01[k] <- spg01[[k]]$modularity
  #count = count + 1
  print(Sys.time())
}

max_index01 <- which.max(spg_mod01)
length(spg_mod01)
best_spg01 <- spg01[[max_index01]]
best_spg01_copy<-best_spg01
saveRDS(spg01, paste0("./data/spinglass",cancer_type,"_spg01.RData"))
saveRDS(spg_mod01, paste0("./data/spinglass",cancer_type,"_spg_mod01.RData"))
best_spg01<-as.list(communities(best_spg01))
saveRDS(best_spg01, paste0("./data/spinglass/",cancer_type,"_melanet_cmt_test.RData"))

print(best_spg01)
#class(best_spg01) array
genes <- as.character(unlist(vertex.attributes(graph_comp)))
#######################################################################################################
new_exp_intgr_tmp01<-new_exp_intgr[expSample_to_keep,]
View(as.data.frame(expSample_to_keep))
new_mty_intgr_tmp01<-new_mty_intgr[mytSample_to_keep,]
snv_count_tmp01<-snv_count[,snvSample_to_keep]
data_cnv_tmp01<-data_cnv_tmp[,dncSample_to_keep]

#View(new_exp_intgr_tmp01)

y1<-which(colnames(new_exp_intgr_tmp01) %in% genes)
y2<-which(colnames(new_mty_intgr_tmp01) %in% genes)

#print(y2)


y3<-which(colnames(snv_count_tmp01) %in% genes)
y4<-which(colnames(data_cnv_tmp01) %in% genes)


new_exp_intgr_tmp02<-new_exp_intgr_tmp01[,y1]
new_mty_intgr_tmp02<-new_mty_intgr_tmp01[,y2]
snv_count_tmp02<-snv_count_tmp01[,y3]
data_cnv_tmp02<-data_cnv_tmp01[,y4]
# dim(new_exp_intgr_tmp02)
# dim(new_mty_intgr_tmp02)
# dim(snv_count_tmp02)
# dim(data_cnv_tmp02)

#print(y4)

#write.csv(clinicalInfo,paste0("./data/tcga_data_processed/",cancer_type,"_clinical_info.csv"))
####################################################################################################




####################################################################################################
##################################section3 community mapping to TCGA genes###############################################

options(stringsAsFactors = F)
setwd("C:/Users/gklizh/Documents/Workspace/code_and_data19/")
#exp_intgr <- readRDS("./data/tcga_data_processed/exp_intgr.RData")
#mty_intgr <- readRDS("./data/tcga_data_processed/mty_intgr.RData")
#error solved:mty_intgr <- readRDS("./data/tcga_data_processed/exp_intgr.RData")

hfeat <- colnames(new_exp_intgr_tmp02) # The features of gene expression data
mfeat <- colnames(new_mty_intgr_tmp02) # The features of DNA methylation data
#melanet_cmt<-readRDS(paste0("./data/spinglass/",cancer_type,"_melanet_cmt_test.RData"))
melanet_cmt<-best_spg01
selected_features <- matrix(ncol = 3, nrow=30000) # Matrix of 3 columns; column1: community, column2: genomic type, column3: mapped component
len <- 0
j <- 1
indx <- NULL
for(i in 1:length(melanet_cmt)){
  cmt = melanet_cmt[[i]]
  # Start mapping
  # Gene expression
  j = j + length(indx)
  indx = NULL
  indx = which(hfeat %in% cmt)
  if (length(indx) != 0){
    len = len +length(indx)
    selected_features[j:len,1] = i
    selected_features[j:len,2] = 1
    selected_features[j:len,3] = indx
  }
  
  # DNA methylation
  j = j+length(indx)
  indx = NULL
  indx = which(mfeat %in% cmt)
  if (length(indx) != 0){
    len = len + length(indx)
    selected_features[j:len,1] = i
    selected_features[j:len,2] = 2
    selected_features[j:len,3] = indx
  }
}
print(cancer_type)
selected_features <- na.omit(selected_features)
write.csv(selected_features, paste0("./data/python_related/data/",cancer_type,"_selected_features01.csv"), row.names = F)
write.csv(colnames(new_exp_intgr_tmp02), paste0("./data/python_related/data/",cancer_type,"_exp_feature_names01.csv"), row.names = F)
write.csv(colnames(new_mty_intgr_tmp02), paste0("./data/python_related/data/",cancer_type,"_mty_feature_names01.csv"), row.names = F) 
write.csv(new_exp_intgr_tmp02, paste0("./data/python_related/data/",cancer_type,"_exp_intgr01.csv"))
write.csv(new_mty_intgr_tmp02, paste0("./data/python_related/data/",cancer_type,"_mty_intgr01.csv"))
check_new_mty_intgr_tmp02<-read.csv(paste0("./data/python_related/data/",cancer_type,"_mty_intgr01.csv"))
print(cancer_type)
nrow(check_new_mty_intgr_tmp02)
#print(best_spg01)
#print(cancer_type)

######################################### section 4 You need to get score profile to proceed ##############################################################################
##### Patient stratification
options(stringsAsFactors = F)

library(NbClust)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(maftools)
print(cancer_type)
#cancer_type<-"BRCA"
setwd("C:/Users/gklizh/Documents/Workspace/code_and_data19/")
print(cancer_type)
#samples <- readRDS(paste0("./data/tcga_data_processed/",cancer_type,"_samples.RData"))


#length(samples_f2)
cmtScores <- read.csv(paste0("./data/python_related/result/community/combine_",cancer_type,"_score_profile_test01.csv"), check.names = F, header = F)
#dim(cmtScores)
#View(cmtScores)
#exp_intgr <- readRDS("./data/tcga_data_processed/SKCM_exp_intgr_all01.RData")
#mty_intgr <- readRDS("./data/tcga_data_processed/SKCM_mty_intgr_all01.RData")
#snv_intgr <- readRDS("./data/tcga_data_processed/SKCM_snv_intgr_test06.RData")
#cnv_intgr <- readRDS("./data/tcga_data_processed/SKCM_cnv_intgr.RData")

#clinicalInfo <- readRDS("./data/tcga_data_processed/SKCM_clinical_info_test06.RData")
#write.csv(as.data.frame(clinicalInfo),"./data/tcga_data_processed/SKCM_clinical_info_test06.csv")
# 
# exp_intgr<-
# mty_intgr<-
# snv_intgr<-
# cnv_intgr<-

#therapy <- readRDS(paste0("./data/tcga_data/",cancer_type,"_therapy_test06.RData"))
#radiation <- readRDS(paste0("./data/tcga_data/",cancer_type,"_radiation_test06.RData"))
#melanet_cmt <- readRDS(paste0("./data/spinglass/",cancer_type,"_melanet_cmt_test.RData"))

#View(therapy)

#print(melanet_cmt)
#dim(samples_f2)
#print(samples_f2)
#View(clinicalInfo)
clinicalInfo_tmp<-new_clinicalInfo

View(clinicalInfo_tmp)
#View(new_clinicalInfo)
row.names(cmtScores) <- samples_f2
colnames(cmtScores) <- paste0("cmt", 1:ncol(cmtScores))
saveRDS(cmtScores, paste0("./data/",cancer_type,"_community_scores.RData"))

#print(cmtScores)
### Determine the best number of clusters
nc <- NbClust(scale(cmtScores), distance = "euclidean", min.nc = 2, max.nc = 10, method = "complete", index = "all")

####################################################################################

if (!require(kernlab)) install.packages("kernlab")
library(kernlab)  # 加载kernlab包
# 可视化（如果数据是二维或三维的话）
library(ggplot2)

if (!require("Rtsne")) install.packages("Rtsne")
library(Rtsne)


n_rows<-nrow(cmtScores)
normalized_data<-scale(cmtScores)


set.seed(42)  # 为了可重复性设置随机种子

tsne_results <- Rtsne(normalized_data, dims = 2, perplexity = ceiling(sqrt(n_rows)), check_duplicates = FALSE, pca = FALSE, verbose = TRUE)


tsne_data <- as.data.frame(tsne_results$Y)

# 使用 ggplot2 进行可视化

ggplot(tsne_data, aes(x = V1, y = V2)) +
  geom_point(alpha = 0.6) +  # 点的透明度
  theme_minimal() +
  labs(title = "t-SNE Visualization of cmtScores")


#-------------------------------------------------------------------------------------------
tsne_results01 <- Rtsne(normalized_data, dims = 3, perplexity =45 , check_duplicates = FALSE, pca = FALSE, verbose = TRUE)
tsne_data01 <- as.data.frame(tsne_results01$Y)
tsne_data_knn<-as.data.frame(tsne_results01$Y)
tsne_dbscan<-as.data.frame(tsne_results01$Y)
# ggplot(tsne_data01, aes(x = V1, y = V2,z=v3)) +
#   geom_point(alpha = 0.6) +  # 点的透明度
#   theme_minimal() +
#   labs(title = "t-SNE Visualization of cmtScores")



row_names <- rownames(cmtScores)
rownames(tsne_data01) <- row_names

plot_ly(data = tsne_data01, x = ~V1, y = ~V2, z = ~V3, type = 'scatter3d', mode = 'markers',
        marker = list(size = 5, opacity = 0.6)) %>%
  layout(title = "3D t-SNE Visualization of cmtScores",
         scene = list(xaxis = list(title = 'Component 1'),
                      yaxis = list(title = 'Component 2'),
                      zaxis = list(title = 'Component 3')))

print(tsne_data01)
#sc_result <- specc(scale(cmtScores), centers = 3)


# 计算欧式距离矩阵
distance_matrix <- as.matrix(dist(tsne_data01))

# 转换为相似性矩阵，sigma 是高斯核的宽度参数
sigma <- mean(distance_matrix)  # 示例，实际可能需要调整
similarity_matrix <- exp(-distance_matrix^2 / (2 * sigma^2))
k <- 4
spectral_clust <- specc(similarity_matrix, centers = k)
clusters <- spectral_clust@.Data

#rownames(final_data) <- row_names
tsne_data01$SampleID <- row_names
tsne_data01$cluster <- as.factor(clusters)  # 添加聚类结果到数据框

# plot_ly(data = tsne_data01, x = ~V1, y = ~V2, z = ~V3, color = ~cluster, type = 'scatter3d', mode = 'markers',
#         marker = list(size = 5, opacity = 0.6)) %>%
#   layout(title = "3D t-SNE with Spectral Clustering",
#          scene = list(xaxis = list(title = 'Component 1'),
#                       yaxis = list(title = 'Component 2'),
#                       zaxis = list(title = 'Component 3')))
# 
# print(tsne_data01)
#---------------------------------------------

plot_3d <- plot_ly(data = tsne_data01, x = ~V1, y = ~V2, z = ~V3,
                   text = ~SampleID,  # 显示行名称
                   type = 'scatter3d', mode = 'markers',
                   marker = list(size = 5, opacity = 0.6),
                   color = ~cluster,  # 根据聚类结果上色
                   hoverinfo = 'text+name')  # 显示文本和簇名称

plot_3d <- plot_3d %>%
  layout(title = "3D Visualization of Spectral Clustering",
         scene = list(xaxis = list(title = 'Component 1'),
                      yaxis = list(title = 'Component 2'),
                      zaxis = list(title = 'Component 3')))

# 显示图形
plot_3d


View(tsne_data01)
##############################k-means#####################################

#--------------------------------------------------------
k1 <- 4  # 设定聚类数量
kmeans_result <- kmeans(tsne_data_knn, centers = k1)


if (!require("scatterplot3d")) install.packages("scatterplot3d")
library(scatterplot3d)

# 查看聚类结果
print(kmeans_result$cluster)

# 可视化
plot(tsne_data_knn, col = kmeans_result$cluster, pch = 19, cex = 1.5)
points(kmeans_result$centers, col = 1:k, pch = 8, cex = 2)
#---------------------------------------------------------------------------------
s3d <- scatterplot3d(tsne_data$V1, tsne_data$V2, tsne_data$V3, pch=19, cex.symbols=1.5, color=kmeans_result$cluster)

# 标记聚类中心
s3d$points3d(kmeans_result$centers[,1], kmeans_result$centers[,2], kmeans_result$centers[,3], col=1:k, pch=8, cex=2)



if (!require("rgl")) install.packages("rgl")
library(rgl)

plot3d(tsne_data$V1, tsne_data$V2, tsne_data$V3, col = rainbow(3)[kmeans_result$cluster], size = 5)

# 添加聚类中心点
centers <- kmeans_result$centers
spheres3d(centers[,1], centers[,2], centers[,3], radius = 0.1, color = "black", alpha = 0.8)
# 将三维视图转换为交互式网页组件
if (!require("htmlwidgets")) install.packages("htmlwidgets")
library(htmlwidgets)
rglwidget()


#-------------------------------------------------------------
# 创建一个简单的三维散点图
plot3d(x = tsne_data_knn$V1, y = tsne_data_knn$V2, z = tsne_data_knn$V3,
       col = kmeans_result$cluster, size = 3)

# 添加 Kmeans 聚类中心
points3d(kmeans_result$centers[,1], kmeans_result$centers[,2], kmeans_result$centers[,3],
         col = "red", size = 10)

# 添加坐标轴标签
#xlabel3d("Component 1")
#ylabel3d("Component 2")
#zlabel3d("Component 3")

# 添加标题
#title3d("3D Visualization of Kmeans Clustering")

# 添加图例
legend3d("topright", legend = unique(kmeans_result$cluster), col = 1:length(unique(kmeans_result$cluster)), pch = 16, cex = 1.5)

# 旋转和交互
rglwidget()
#-----------------------------------------------
install.packages("dbscan")
library(dbscan)
dbscan_result <- dbscan(tsne_dbscan, eps = 0.5, MinPts = 5)

# 获取聚类结果
clusters <- dbscan_result$cluster

# 如果您想要查看噪声点的数量，可以查看 cluster 值为 0 的数据
noise_points <- sum(clusters == 0)
print(paste("Number of noise points:", noise_points))

# 可视化
library(rgl)
plot3d(tsne_dbscan, col = clusters, size = 5)

###############################################################


best_clusters <- nc$Best.nc[1, ]  # 获取所有指数的推荐集群数量

# 计算每个集群数量的出现频率
cluster_counts <- table(best_clusters)

View(cluster_counts)
# 找到频率最高的集群数量
most_frequent_cluster <- as.numeric(names(which.max(cluster_counts)))

# 打印最频繁推荐的集群数量
print(most_frequent_cluster)

pdf(paste0("./figure/",cancer_type,"_best_number_of_clusters01.pdf", width = 7, height = 7))
ggplot(data.frame(cluster = factor(nc$Best.nc[1,])), aes(x = cluster)) + 
  geom_bar(stat = "count", fill = "#C1BFBF") + 
  labs(x = "Number of clusters", y = "Number of criteria", title = "Number of clusters chosen by 26 criteria") + 
  theme(text = element_text(size = 18), plot.title = element_text(hjust = 0.5, face = "bold"), panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black')) + 
  scale_y_continuous(breaks = seq(0,14,2), limits = c(0,14))
dev.off()

#View(clinicalInfo)
### Community scores of clustered patients
tumorType <- clinicalInfo_tmp[,"shortLetterCode"]



#tumorStage <- clinicalInfo_tmp[,"tumor_stage"]
clinicalInfo_col<-as.data.frame(colnames(clinicalInfo))
write.csv(clinicalInfo_col,'./data/clinicalInfo_col.csv')


# tumorStage[-grep("^stage", tumorStage)] <- NA
# tumorStage <- gsub("^stage ", "", tumorStage)
# tumorStage <- gsub("[a-c]$", "", tumorStage)
# 
# tumorStage <- clinicalInfo_tmp[,"tumor_stage"]
# tumorStage[-grep("^stage", tumorStage)] <- NA
# tumorStage <- gsub("^stage ", "", tumorStage)
# tumorStage <- gsub("[a-c]$", "", tumorStage)

therapy <- therapy[3:nrow(therapy),]
ifTherapy <- substr(samples, 1, 12) %in% therapy$bcr_patient_barcode
ifTherapy <- ifelse(ifTherapy, "Yes", "No")

radiation <- radiation[3:nrow(radiation),]
ifRadiation <- substr(samples, 1, 12) %in% radiation$bcr_patient_barcode
ifRadiation <- ifelse(ifRadiation, "Yes", "No")

therapy_type <- sapply(substr(samples, 1, 12), function(x){paste(sort(unique(therapy$pharmaceutical_therapy_type[therapy$bcr_patient_barcode == x])),collapse = ";")})

tr_df <- data.frame(sample = samples, theray = ifTherapy, radiation = ifRadiation, therapy_type = therapy_type)
saveRDS(tr_df, paste0("./data/",cancer_type,"_therapy_radiation_df.RData"))

tumorType_col_fun <- c("TM" = "#CC79A7", "TP" = "#0072B2")
tumorStage_col_fun <- c("0" = "#FAFCC2", "i" = "#FFEFA0", "ii" = "#FFD57E", "iii" = "#FCA652", "iv" = "#AC4B1C")
ifTherapy_col_fun <- c("Yes" = "red", "No" = "gray")
ifRadiation_col_fun <- c("Yes" = "red", "No" = "gray")
#cancer_type<-"COAD"
#print(most_frequent_cluster)
print(cancer_type)
#topAnno <- HeatmapAnnotation(Therapy = ifTherapy, Radiation = ifRadiation, `Tumor type` = tumorType,  col = list(Therapy = ifTherapy_col_fun, Radiation = ifRadiation_col_fun, `Tumor type` = tumorType_col_fun, `Tumor stage` = tumorStage_col_fun), border = T, show_annotation_name = T)
ht = Heatmap(t(scale(cmtScores)), 
             name = "Community score", 
             show_column_names = F,
             # top_annotation = topAnno,
             clustering_distance_columns = "euclidean",
             #clustering_distance_columns = "spearman",
             clustering_method_columns = "complete",
             column_split = most_frequent_cluster,
             #column_split = 2,
             column_title = "%s",
)#`Tumor stage` = tumorStage, removed

pdf(paste0("./figure/",cancer_type,"_heatmap_cmtScores01_05.pdf", width = 7, height = 7))
draw(ht, merge_legends = TRUE)
dev.off()


ht = draw(ht)
rowOrder <- row_order(ht)
colOrder <- column_order(ht)
print(rowOrder)
print(colOrder)
View(as.data.frame(rowOrder))
View(as.data.frame(colOrder))

#View(ht)
#View(colOrder)
#View(rowOrder)

#View(ht)
##############################################################################################
samplePartition <- data.frame(cluster = rep(1:length(colOrder), lengths(colOrder)), sampleID = unlist(colOrder))
samplePartition <- samplePartition[order(samplePartition$sampleID),]
saveRDS(samplePartition, paste0("./data/",cancer_type,"_sample_partition01.RData"))
#View(data.frame(cluster = rep(1:length(colOrder), lengths(colOrder)), sampleID = unlist(colOrder)))
#View(samplePartition)
#nrow(samplePartition)
#View(samplePartition)
################################################################################
##
##
survivalInfo_tmp_test<-survivalInfo
survivalInfo_tmp_test$PatientID <- rownames(survivalInfo)
#View(as.data.frame(substr((survivalInfo_tmp_test$PatientID),1,15)))
survivalInfo_tmp_test$PatientID<-substr((survivalInfo_tmp_test$PatientID),1,16)
#View(tsne_data01)
combined_data <- merge(survivalInfo_tmp_test, tsne_data01, by.x = "PatientID", by.y = "SampleID")
#colnames(combined_data)[colnames(data) == "cluster"] <- "hc"
#View(combined_data)
#View(survivalInfo_tmp_test)
###################################Section 5################################################################


##### Survival analysis
options(stringsAsFactors = F)
library(TCGAbiolinks)
#cancer_type<-"BLCA"
setwd("C:/Users/gklizh/Documents/Workspace/code_and_data19/")
print(cancer_type)
#setwd("C:/Users/gklizh/Documents/Workspace/code_and_data19")
#clinicalInfo <- readRDS(paste0("./data/tcga_data_processed/",cancer_type,"_new_clinicalInfo_test06.RData"))
#View(new_clinicalInfo)
#clinicalInfo_tmp<-new_clinicalInfo
#samplePartition <- readRDS(paste0("./data/",cancer_type,"_sample_partition01.RData"))
#mty_colData<-readRDS(paste0("./data/tcga_data/",cancer_type,"_mty_colData.RData"))
#mty_colData<-paste0("./data/tcga_data/",cancer_type,"_mty_colData.RData")
#clinicalInfo <- mty_colData
#View(as.data.frame(clinicalInfo))
#row.names(clinicalInfo) <- clinicalInfo$sample
#clinicalInfo <- clinicalInfo[samples,]


#saveRDS(clinicalInfo, paste0("./data/tcga_data_processed/",cancer_type,"_clinical_info_test06.RData"))
#write.csv(clinicalInfo,paste0("./data/tcga_data_processed/",cancer_type,"_clinical_info.csv"))
#dim(samplePartition)
#View(clinicalInfo)
#dim(clinicalInfo)
#View(samplePartition)
#View(samplePartition)
#dim(samplePartition)
#View(clinicalInfo_tmp)
#nrow(clinicalInfo_tmp)
survivalInfo <- clinicalInfo_tmp[,c("shortLetterCode", "vital_status","days_to_death","days_to_last_follow_up")]
# nrow(samplePartition)
# [1] 294
# > nrow(survivalInfo)
# [1] 308
# > nrow(samplePartition)
# [1] 294
# > nrow(survivalInfo)
# [1] 308

View(survivalInfo)

survivalInfo_tmp <- survivalInfo[1:nrow(samplePartition),]

View(survivalInfo_tmp)

survivalInfo_tmp$hc <- samplePartition$cluster # ???˲??ξ???????


combined_data$V1<-NULL
combined_data$V2<-NULL
combined_data$V3<-NULL


View(combined_data)


#View(samplePartition)
#View(survivalInfo_tmp)

# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage 0")] <- 0
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage i")] <- 1
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage ia")] <- 1
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage ib")] <- 1
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage ic")] <- 1
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage ii")] <- 2
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage iia")] <- 2
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage iib")] <- 2
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage iic")] <- 2
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage iii")] <- 3
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage iiia")] <- 3
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage iiib")] <- 3
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage iiic")] <- 3
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "stage iv")] <- 4
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "not reported")] <- "NA"
# survivalInfo$tumor_stage[which(survivalInfo$tumor_stage == "i/ii nos")] <- "NA"
# survivalInfo$tumor_stage[which(is.na(survivalInfo$tumor_stage))] <- "NA"
#print(cancer_type)
TCGAanalyze_survival(survivalInfo_tmp, clusterCol = "hc", color = c("#33A02C","#1F78B4","#E31A1C","#AC4B1C","#FFD57E","#EE3388","#345565"), filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis01_update.pdf"), conf.int = F, width = 7, height = 7)

TCGAanalyze_survival(survivalInfo_tmp, clusterCol = "shortLetterCode", filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis_tumorType01_update.pdf"), conf.int = F, width = 7, height = 7)

#TCGAanalyze_survival(survivalInfo, clusterCol = "tumor_stage", filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis_tumorStage.pdf"), conf.int = F, width = 7, height = 7)
#print(cancer_type)
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
TCGAanalyze_survival(combined_data, clusterCol = "cluster", color = c("#33A02C","#1F78B4","#E31A1C","#AC4B1C","#FFD57E","#EE3388","#345565"), filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis01_update.pdf"), conf.int = F, width = 7, height = 7)

TCGAanalyze_survival(combined_data, clusterCol = "shortLetterCode", filename = paste0("./figure/surv_analysis/",cancer_type,"_survival_analysis_tumorType01_update.pdf"), conf.int = F, width = 7, height = 7)


######################################
# 
# 
# data_scaled <- scale(cmtScores)
# data_vis <- as.data.frame(scale(cmtScores))
# colnames(data_vis) <- colnames(cmtScores)
# data_vis$cluster <- as.factor(clusters)
# #View(data_vis)
# 
# #View(cmtScores)
# data_vis$cluster <- as.factor(clusters)
# # ggplot(data_vis, aes(x = V1, y = V2, color = cluster)) +
# #   geom_point(alpha = 0.6) +
# #   theme_minimal() +
# #   labs(title = "Spectral Clustering Results")
# 
# normalized_data<-scale(cmtScores)
# 
# n_rows<-nrow(cmtScores)
# tsne_results <- Rtsne(normalized_data, dims = 3, perplexity = ceiling(sqrt(n_rows)), verbose = TRUE)
# 
# tsne_data <- data.frame(tsne_results$Y, Cluster = as.factor(clusters))
# 
# #print(tsne_data$Cluster)
# sc_result <- specc(scale(cmtScores), centers = 3)
# 
# # 获取聚类结果
# clusters <- sc_result@.Data
# 
# print(clusters)
# 
# #绘制三维 t-SNE 结果
# plot_ly(data = tsne_data, x = ~X1, y = ~X2, z = ~X3, color = ~Cluster, colors = RColorBrewer::brewer.pal(length(unique(tsne_data$Cluster)), "Set1"),
#         type = "scatter3d", mode = "markers") %>%
#   layout(title = "3D t-SNE Visualization of Clusters",
#          scene = list(xaxis = list(title = 't-SNE 1'),
#                       yaxis = list(title = 't-SNE 2'),
#                       zaxis = list(title = 't-SNE 3')))
# 
# 



#print(tsne_data)
#--------------------------------------------------------------------------------------------#
# sc_result <- specc(tsne_results$Y, centers = 3)  # 假设将数据聚为 4 类
# 
# 
# print(sc_result)
# # 获取聚类结果
# clusters <- sc_result@.Data
# print(clusters)
# 
# # 将聚类结果添加到数据框中
# #tsne_data <- cbind(sc_result, Cluster_spectral = as.factor(clusters))
# 
# #print(tsne_data)
# #(print(tsne_data$Cluster_spectral))
# 
# 
# tsne_data <- data.frame(tsne_results$Y, Cluster = as.factor(clusters))
# 
# # 将谱聚类的结果添加到数据框中
# tsne_data <- cbind(tsne_data, Cluster_spectral = as.factor(clusters))
# 
# print(tsne_data$Cluster_spectral)
# 
# # 可视化谱聚类结果
# plot_ly(data = tsne_data, x = ~X1, y = ~X2, z = ~X3, color = ~Cluster_spectral, 
#         colors = RColorBrewer::brewer.pal(length(unique(tsne_data$Cluster_spectral)), "Set1"), 
#         type = "scatter3d", mode = "markers") %>%
#   layout(title = "3D t-SNE Visualization with Spectral Clustering",
#          scene = list(xaxis = list(title = 't-SNE 1'),
#                       yaxis = list(title = 't-SNE 2'),
#                       zaxis = list(title = 't-SNE 3')))
# plot <- plot_ly(data = tsne_data01, x = ~V1, y = ~V2, z = ~V3, 
#                 text = ~SampleID,  # 指定悬停文本显示行名称
#                 type = 'scatter3d', 
#                 mode = 'markers',
#                 marker = list(size = 5, opacity = 0.6),
#                 hoverinfo = 'text')  # 只显示 text 信息，不显示 x, y, z 值
# 
# plot <- plot %>%
#   layout(title = "3D t-SNE Visualization",
#          scene = list(xaxis = list(title = 'Component 1'),
#                       yaxis = list(title = 'Component 2'),
#                       zaxis = list(title = 'Component 3')))
# 
# # 显示图形
# plot

# combined_data <- combined_data %>%
#   select(-V1, -V2, -V3)
