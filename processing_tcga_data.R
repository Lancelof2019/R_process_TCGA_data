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
#####
library(SummarizedExperiment)
library(dplyr)
library(lubridate)
#### Collect gene expression data

cancer_type<-"SKCM"
setwd("/Workspace/code_and_data18")
query_exp <- GDCquery(project = "TCGA-SKCM", 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")#
GDCdownload(query_exp) 
data_exp <- GDCprepare(query_exp) 
exp_assay <- as.data.frame(assay(data_exp)) # Gene expression matrix

exp_assay_temp<-as.data.frame(assay(data_exp)) # Gene expression matrix

dim(exp_assay_temp)#60660 x 473

exp_colData <- as.data.frame(colData(data_exp)) # Patient annotation(472 patients)
exp_rowRanges <- as.data.frame(rowRanges(data_exp)) # Gene annotation

#saveRDS(exp_rowRanges,"./data/tcga_data/exp_rowRanges.RData")

#test_exp_rowRanges<-readRDS("./data/tcga_data/exp_rowRanges.RData")

#View(test_exp_rowRanges)

#### Collect DNA methylation data
query_mty <- GDCquery(project = "TCGA-SKCM",
                      data.category = "DNA Methylation",
                      platform = "Illumina Human Methylation 450",
                      data.type = "Methylation Beta Value")
#View(query_mty)

GDCdownload(query_mty)
data_mty <- GDCprepare(query_mty)

#View(data_mty)


#class(data_mty)
#dim(data_mty)
#View(data_mty)
#assays(data_mty)
#View(rowRanges(data_mty))
mty_assay <- as.data.frame(assay(data_mty)) # DNA methylation matrix
mty_colData <- as.data.frame(colData(data_mty)) # Patient annotation(475 patients)

mty_rowRanges <- as.data.frame(rowRanges(data_mty)) # cg probe annotation
#######################Save files##########################################
saveRDS(exp_colData,"./data/tcga_data/exp_colData.RData")
saveRDS(exp_rowRanges,"./data/tcga_data/exp_rowRanges.RData")


saveRDS(mty_colData, paste0("./data/tcga_data/",cancer_type,"_mty_colData.RData"))
saveRDS(mty_assay, paste0("./data/tcga_data/",cancer_type,"_mty_assay.RData"))
saveRDS(mty_rowRanges, paste0("./data/tcga_data/",cancer_type,"_mty_rowRanges.RData"))

#test_mty_rowRangess<-readRDS("./data/tcga_data/mty_rowRanges.RData")
#test_mty_assay<-readRDS("./data/tcga_data/mty_assay.RData")

#saveRDS(mty_assay,"./data/tcga_data/mty_assay.RData")

#### Collect gene mutation data
#maf <- GDCquery_Maf("SKCM", pipelines = "muse") # Mutation annotation
#saveRDS(maf, "./data/tcga_data/maf.RData")

#cancer_type<-"SKCM"

query <- GDCquery(
     project = paste0("TCGA-", cancer_type),
    data.category = "Simple Nucleotide Variation", 
     access = "open",
     data.type = "Masked Somatic Mutation", 
     #workflow.type = "MuSE Variant Aggregation and Masking"
      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
   )
GDCdownload(query)

maf <- GDCprepare(query)

filename1 <- paste0("./data/tcga_data/", cancer_type, "_maf.RData")

#saveRDS(maf, "./data/tcga_data/maf.RData")
saveRDS(maf, filename1)
#display_data<-readRDS(filename1)
#View(display_data)

#filecompare<-readRDS("./data/tcga_data/maf.RData")
#View(filecompare)





#### Collect copy number variation data
#query_cnv <- GDCquery(project = "TCGA-SKCM",
                      #data.category = "Copy Number Variation",
                     # data.type = "Gene Level Copy Number Scores")

#str(query_cnv)

#View(getResults(query_cnv))

#query_cnv1 <- GDCquery(project = "TCGA-SKCM",data.category = "Copy Number Variation",data.type="Gene Level Copy Number")
#
###########################################################
#query_cnv2 <- GDCquery(project = paste0("TCGA-", cancer_type),
                      #data.category = "Copy Number Variation")

# 获取查询结果
#results2 <- getResults(query_cnv2)
#View(results2)
# 从查询结果中获取 data_type 列的唯一值

#data_types2 <- unique(results2$cases)
#print(data_types2)
#View(data_types2)
# 使用正则表达式匹配
#matched_data_types <- data_types2[grepl("Gene Level Copy Number.*$", data_types2)]
#print(matched_data_types)

# 选择一个 data_type 值（例如，这里选择第一个）
#selected_data_type <- data_types2[1]
#print(selected_data_type)
# 创建具有所选数据类型的新查询
query_cnv <- GDCquery(project = paste0("TCGA-", cancer_type),
                      data.category = "Copy Number Variation",
                      data.type = "Gene Level Copy Number")

GDCdownload(query_cnv)
#class(query_cnv$results)
#str(query_cnv)

#class(query_cnv)
#colnames(query_cnv)
#colnames(getResults(query_cnv))
#View(getResults(query_cnv))
#query_results<-getResults(query_cnv)
#View(query_cnv)
#colnames(query_results)
#########################################################
#colnames(query_cnv$results)
#str(query_cnv$results)

potential_df <- query_cnv$results[[1]]

View(potential_df)
#str(potential_df)
#colnames(potential_df)
#class(potential_df)
#results_df <- as.data.frame(query_cnv$results)
potential_df$created_datetime <- ymd_hms(potential_df$created_datetime)

#colnames(results_df)
#results_df_unique <- results_df %>% 
  #distinct(cases, .keep_all = TRUE)
#unique_cases_df <- potential_df %>%
  #distinct(cases, .keep_all = TRUE)
potential_df$timestamp <- as.numeric(potential_df$created_datetime)

#unique_cases_df <- potential_df %>%
  #distinct(cases, .keep_all = TRUE) %>%
  #group_by(cases) %>%
  #slice_max(order_by = created_datetime) %>%
  #slice_max(order_by = as_datetime(created_datetime, tz = "GMT")) %>%
  #slice_max(order_by = c(year(created_datetime), month(created_datetime), day(created_datetime), hour(created_datetime), minute(created_datetime), second(created_datetime))) %>%
  #ungroup()

#######################################################################



#unique_cases_df <- potential_df %>%
 # distinct(cases, .keep_all = TRUE) %>%
 # group_by(cases) %>%
 # slice_max(order_by = timestamp) %>%
  #ungroup() %>%
 # select(-timestamp)
#####################################################
unique_cases_df <- potential_df %>%
  distinct(cases, .keep_all = TRUE) %>%
 group_by(cases) %>%
 mutate(datetime = ymd_hms(created_datetime)) %>%  # 将created_datetime转换为日期时间对象
 filter(datetime == max(datetime)) %>%            # 过滤出最新的日期时间
 ungroup() %>%
 select(-datetime)%>%
 select(-timestamp)


#############################################################

#####################################################
query_cnv$results[[1]] <- unique_cases_df
##########################################################################

#########################################################


saveRDS(query_cnv,"./data/tcga_data/query_cnv.RData")

View(query_cnv)
write.csv(query_cnv$results[[1]], file = "./data/tcga_data/query_results13.csv", row.names = TRUE)

############################################################
query_cnv_copy<-query_cnv
View(query_cnv_copy)
#colnames(query_cnv$results[[1]])
query_cnv_copy$results[[1]] <- query_cnv_copy$results[[1]] %>%
              distinct(sample.submitter_id, .keep_all = TRUE)
#colnames(query_cnv_copy)
#View(query_cnv_copy)
#dim(query_cnv_copy)
#View(query_cnv_copy)

View(query_cnv_copy)
dim(query_cnv_copy)
GDCdownload(query_cnv_copy)
#print(query_cnv)
#GDCdownload()
#rownames(query_cnv$results[[1]])
         
#View(query_cnv)

#str(query_cnv)
#colnames(query_cnv)

#dim(query_cnv_copy)
data_cnv <- GDCprepare(query_cnv_copy)

########################################################################################################################GDC
data_cnv_savage <- GDCprepare(query_cnv_copy)
#saveRDS(data_cnv_savage,"./data/tcga_data/data_cnv_savage_GDCprepare.RData")
#data_cnv_savage<-readRDS("./data/tcga_data/data_cnv_savage_GDCprepare.RData")
View(as.data.frame(colData(data_cnv_savage)))

View(as.data.frame(rowData(data_cnv_savage)))
#class(data_cnv_savage)
#View(data_cnv_savage)
#class(data_cnv_savage)
#rowdata_data_cnv<-rowData(data_cnv_savage)
#View(rowdata_data_cnv)
#head(rowdata_data_cnv)
#class(data_cnv)
#saveRDS(data_cnv,"./data/tcga_data/data_cnv_GDCprepare.RData")

#saveRDS(rowRanges(data_cnv),"./data/tcga_data/rowRanges_data_cnv_GDCprepare.RData")
#data_cnv<-readRDS("./data/tcga_data/rowRanges_data_cnv_GDCprepare.RData")
#data_cnv_savage<-saveRDS(data_cnv_savage,"./data/tcga_data/data_cnv_savage_GDCprepare.RData")
#View(data_cnv_savage)

#write.csv(data_cnv,"./data/tcga_data/data_cnv_GDCprepare.csv")
#write.csv(rowRanges(data_cnv),"./data/tcga_data/rowRanges_data_cnv_GDCprepare.csv")

#data_cnv_colData <- as.data.frame(colData(data_cnv))
#View(data_cnv_colData)
#data_cnv_copy<-data_cnv
#View(data_cnv_copy)
#View(data_cnv)
########################################################################################################################GDC
#class(data_cnv)
data_temp<-data_cnv_savage
#data_temp02<-data_cnv_savage
#View(data_temp)
data_temp <- as.data.frame(assay(data_temp))
saveRDS(data_temp,"./data/tcga_data/data_temp_savage.csv")
#data_temp02 <- as.data.frame((data_temp02))
#View(data_cnv)
#class(data_cnv)

#saveRDS(data_cnv,"./data/tcga_data/data_cnv_df.RData")
#write.csv(data_cnv,"./data/tcga_data/data_cnv_df.csv")
#data_cnv<-readRDS("./data/tcga_data/data_cnv.RData")

#View(data_cnv)
#### Collect clinical data
clinical <- GDCquery_clinic(project = "TCGA-SKCM", type = "clinical")

#### Collect clinical radiation and drug therapy data
query <- GDCquery(project = "TCGA-SKCM", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab"
)
GDCdownload(query)
sckm.tab.all <- GDCprepare(query)

therapy <- sckm.tab.all$clinical_drug_skcm # clinical drug therapy
saveRDS(therapy, paste0("./data/tcga_data/",cancer_type,"_therapy.RData"))
        #"./data/tcga_data/therapy.RData")
therapy$pharmaceutical_therapy_type # Therapy types


radiation <- sckm.tab.all$clinical_radiation_skcm # Clinical radiation therapy
#saveRDS(radiation, "./data/tcga_data/radiation.RData")

saveRDS(radiation, paste0("./data/tcga_data/",cancer_type,"_radiation.RData"))

########################################################### Process gene expression data#################################
#View(exp_assay)
colnames(exp_assay) <- substr(colnames(exp_assay), 1, 16)

#View(exp_assay)
#class(exp_assay)
#View(row.names(exp_assay))
#View(exp_rowRanges)
exp_assay$SYMBOL <- exp_rowRanges[row.names(exp_assay), "gene_name"] 
#View(exp_assay)
saveRDS(exp_assay,"./data/tcga_data/exp_assay.RData")
#colnames(exp_assay)
#rownames(exp_assay)
#dim(exp_assay)
#print(colnames(exp_assay)[1])

####################################################################################################################################################

exp_assay$SYMBOL <- exp_rowRanges[row.names(exp_assay), "gene_name"] 

num_rows <- nrow(exp_assay)
#print(num_rows)
# 获取列数
num_cols <- ncol(exp_assay)
dim(exp_assay)
#print(num_cols)
#dim(exp_assay)
#View(exp_assay)
exp_assay_test1 <- exp_assay[,c(num_cols,1:num_cols-1)]
#View(exp_assay_test1)
exp_assay_test2 <- as.matrix(exp_assay_test1)
#View(exp_assay_test2)
row.names(exp_assay_test2) <- exp_assay_test2[,1]
#View(exp_assay_test2)
#dim(exp_assay_test2)
exp_assay_test2 <- exp_assay_test2[,2:num_cols] # Convert Ensemble ID to corresponding gene symbols
#View(exp_assay_test2)
exp_assay_test2 <- matrix(as.numeric(exp_assay_test2), nrow = nrow(exp_assay_test2), dimnames = list(row.names(exp_assay_test2), colnames(exp_assay_test2)))
exp_assay_test2 <- avereps(exp_assay_test2) # If the gene corresponds to multiple gene expression values, take the average value,remove duplicated gene name

#View(exp_assay_test2)
#dim(exp_assay_test2)
exp_assay<-exp_assay_test2
dim(exp_assay)
####################################################################################################################################################
#exp_assay <- exp_assay[,c(473,1:472)]
#exp_assay <- as.matrix(exp_assay)
#row.names(exp_assay) <- exp_assay[,1]
#exp_assay <- exp_assay[,2:473] # Convert Ensemble ID to corresponding gene symbols
####################################################################################
#exp_assay_copy <- exp_assay[,c(474,1:473)]
#View(exp_assay_copy)
#exp_assay_copy <- as.matrix(exp_assay_copy)
#View(exp_assay_copy)
#exp_assay_copy <- exp_assay_copy[,2:474] # Convert Ensemble ID to corresponding gene symbols

#View(exp_assay_copy)
####################################################################################
#exp_assay <- matrix(as.numeric(exp_assay), nrow = nrow(exp_assay), dimnames = list(row.names(exp_assay), colnames(exp_assay)))
#exp_assay <- avereps(exp_assay) # If the gene corresponds to multiple gene expression values, take the average value
####################################################################################

#exp_assay_copy <- matrix(as.numeric(exp_assay_copy), nrow = nrow(exp_assay_copy), dimnames = list(row.names(exp_assay_copy), colnames(exp_assay_copy)))

#dim(exp_assay_copy)
#exp_assay_copy <- avereps(exp_assay_copy) # If the gene corresponds to multiple gene expression values, take the average value
####################################################################################

#### Process DNA methylation data
#class(data_mty)
#View(data_mty)
#dim(data_mty)
data_mty <- subset(data_mty, subset = (rowSums(is.na(assay(data_mty)))==0)) 
#View(getResults(data_mty))
#View(data_mty)
#dim(data_mty)

#mty_assay <- as.data.frame(assay(data_mty)) # DNA methylation matrix
#mty_colData <- as.data.frame(colData(data_mty)) # Patient annotation(475 patients)
#mty_rowRanges <- as.data.frame(rowRanges(data_mty)) # cg probe annotation

#View(mty_rowRanges)
#############################################################################################################################crisis
#data_mty_temp <-data_mty



#pre_subset_rowRanges <- rowRanges(data_mty_temp)
#pre_subset_rowRanges_df <- as.data.frame(pre_subset_rowRanges)
#View(pre_subset_rowRanges_df)

#############################################################################################################################crisis
#View(data_mty_temp)
#dim(data_mty_temp)
#data_mty<- subset(data_mty, subset = (as.data.frame(rowRanges(data_mty))$gene_HGNC!=NA))

data_mty<- subset(data_mty, subset =!is.na(as.data.frame(rowRanges(data_mty))$gene_HGNC))
#View(data_mty)
#colnames(data_mty)
#str(data_mty)
#class(data_mty)
#############################################################################################################################crisis
#View((as.data.frame(rowRanges(data_mty))$seqinfo[[1]]))

#if (!is.null(rowRanges(data_mty))) {
#  message("Row ranges are present.")
#}# else {
 # message("Row ranges are NULL.")
#}
#summary(rowRanges(data_mty))
#head(rowRanges(data_mty))


# Check for gene symbols directly within the rowRanges
#if (!is.null(mcols(rowRanges(data_mty))$Gene_Symbol)) {
 # message("Gene_Symbol column exists in row metadata.")
#} else {
 # message("Gene_Symbol column does not exist in row metadata.")
#}
###########################################################
#pre_subset_rowRanges <- rowRanges(data_mty)
#pre_subset_rowRanges_df <- as.data.frame(pre_subset_rowRanges)
#View(pre_subset_rowRanges_df)
#############################################################

mty_assay <- as.data.frame(assay(data_mty))
#View(mty_assay)
#str(mty_assay)
#rownames(mty_assay)
#colnames(mty_assay)
mty_rowRanges <- as.data.frame(rowRanges(data_mty))
#View(rowRanges(data_mty))
#row_ranges_df <- as.data.frame(rowRanges(data_mty))
#head(row_ranges_df)
#View(mty_rowRanges)

#dim(mty_rowRanges)

mty_symbol <- strsplit(mty_rowRanges$gene_HGNC, split = ";")
#dim(mty_symbol)
#length(mty_symbol)

#mty_symbol <- strsplit(mty_rowRanges$Gene_Symbol, split = ";")
#View(row.names(mty_rowRanges))
names(mty_symbol) <- row.names(mty_rowRanges)


mty_symbol <- lapply(mty_symbol, FUN = function(x){x<-unique(x)})


mty_symbol <- as.matrix(unlist(mty_symbol))

row.names(mty_symbol) <- substr(row.names(mty_symbol), 1, 10)
#View(mty_symbol)

mty_symbol <- data.frame("probe" = row.names(mty_symbol),"SYMBOL" = mty_symbol[,1])#cg14528386.1 ->H19_1
#View(mty_symbol)

#View(mty_assay)
mty_assay$probe <- row.names(mty_assay)
#View(mty_assay)
mty_assay <- merge(mty_assay, mty_symbol, by.x = "probe", by.y = "probe") 

#View(mty_assay)
#dim(mty_assay)


num_rows <- nrow(mty_assay)

# 获取列数
num_cols <- ncol(mty_assay)
#mty_mat<-mty_assay[,2:477]
#mty_mat<-as.matrix(mty_assay)
mty_mat <- as.matrix(mty_assay[,2:num_cols])

#dim(mty_assay)
#View(mty_assay)


colnames(mty_mat) <- substr(colnames(mty_mat), 1, 16)
row.names(mty_mat) <- mty_assay$SYMBOL # Convert probe ID to corresponding gene symbols
#dim(mty_mat)
#View(mty_mat)

mty_mat <- matrix(as.numeric(mty_mat), nrow = nrow(mty_mat), dimnames = list(row.names(mty_mat), colnames(mty_mat)))
mty_mat <- avereps(mty_mat) # If the gene corresponds to multiple methylation values, take the average value

View(mty_mat)
#### Process gene mutation data
maf <- maf[,c("Tumor_Sample_Barcode","Hugo_Symbol","Gene","Variant_Classification")]

#colnames(maf)
#View(maf)
rnames <- unique(maf$Hugo_Symbol)
cnames <- unique(maf$Tumor_Sample_Barcode)
snv_count <- matrix(data = 0, nrow = length(rnames), ncol = length(cnames), dimnames = list(rnames,cnames)) 
#save.image(file = 'my_work01.RData')
#View(snv_count)
#dim(snv_count)
# Calculate the frequency of genes' variants
for(i in 1:nrow(maf)){
  rname <- maf[i,]$Hugo_Symbol
  cname <- maf[i,]$Tumor_Sample_Barcode
  snv_count[rname, cname] <- snv_count[rname,cname] + 1
}
#View(snv_count)
#https://www.notion.so/57782b3efe144e91a93ae602f85a1cbd


colnames(snv_count) <- substr(colnames(snv_count), 1, 16)

#### Process copy number variation data
#data_cnv_tmp<-data_cnv
###############################################################################################
#row_names_df <- data.frame(row_names = row.names(data_cnv_copy))
#write.csv(row_names_df,file="./data/tcga_data/row_names_df.csv",row.names =TRUE)
# 在 RStudio 中查看数据框
#View(row_names_df)

# 获取行名
#row_names <- row.names(data_cnv_copy)

# 检测是否有重复的行名
#duplicated_row_names <- duplicated(row_names)

# 如果有重复的行名，打印出来
#if(any(duplicated_row_names)){
#  print("There are duplicated row names.")
#  print(row_names[duplicated_row_names])
#} else {
#  print("There are no duplicated row names.")
#}

#write.csv(data_cnv_copy,file="./data/tcga_data/data_cnv_copy.csv",row.names =TRUE)

#View(data_cnv_copy)

###############################################################################################!!!!!!!!!!!!!!!
#row.names(data_cnv) <- substr(data_cnv[,1], 1, 15)

data_temp$ensembl<-rownames(data_temp)

rownames(data_temp)<-seq_along(rownames(data_temp))

View(data_temp)

#ncol(data_temp)

#data_temp<data_temp[,c(ncol(data_temp),1:ncol(data_temp)-1)]

#View(data_temp)

data_temp$ensembl<-substr((data_temp$ensembl),1, 15)

na_count <- function(x) sum(is.na(x))

data_temp$na_count <- apply(data_temp, 1, na_count)

data_temp <- data_temp[order(data_temp$ensembl, data_temp$na_count), ]

data_temp <- data_temp[!duplicated(data_temp$ensembl), ]
#remove na_count
data_temp$na_count <- NULL

#dim(data_temp)
#View(data_temp)
#ncol(data_temp)

#data_temp <- unique(data_temp)
#ncol(data_temp)

row.names(data_temp)<-data_temp[,ncol(data_temp)]
View(data_temp)
#View(data_cnv)
#first_col_substr <- substr(data_cnv[,1], 1, 15)
#duplicated_values <- first_col_substr[duplicated(first_col_substr)]
#print(duplicated_values)
#print(row.names(data_cnv_copy))


data_cnv <- data_temp[,4:ncol(data_temp)]
colnames(data_cnv) <- substr(colnames(data_cnv), 1, 16)

#colnames(data_cnv)
###############################################################################################!!!!!!!!!!!!!!!!!!!!
col_data_cnv<-ncol(data_cnv)
#row.names(data_cnv) <- unique(substr(data_cnv[,col_data_cnv-1], 1, 15))

#View(as.data.frame(row.names(data_cnv)))

#############################################################################################################
#############################################################################################################
#new_data <- data.frame(gene_id = substr(data_cnv$gene_id, 1, 15), gene_name = data_cnv$gene_name)
#dim(new_data)
#new_data <- new_data[!duplicated(new_data$gene_id), ]
#dim(new_data)
#View(new_data)

#row.names(new_data)<-new_data[,1]
#new_data_tmp1<-new_data[,-1]
#############################################################################################################
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
View(data_cnv)
cnv_df <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                filters = c("ensembl_gene_id"),
                values = row.names(data_cnv),
                mart = ensembl)

#View(cnv_df)
#dim(cnv_df)

##############################################################
cnv_df <- cnv_df[which(cnv_df$hgnc_symbol !=''),]
#dim(cnv_df)
#View(cnv_df)

#saveRDS(cnv_df,"./data/tcga_data/cnv_df.RData")

#write.csv(cnv_df,"./data/tcga_data/cnv_df.csv")
View(cnv_df)
dim(cnv_df)
cnv_df<-unique(cnv_df)
dim(cnv_df)
#dim(cnv_df)
#View(data_cnv)
colnames(data_cnv)[colnames(data_cnv)=="ensembl"]<-"emsembl"
#View(as.data.frame(colnames(data_cnv)))
#View(data_cnv)
#length(colnames(data_cnv))
unique_cols<-!duplicated(t(data_cnv))
data_cnv<-data_cnv[,unique_cols]
#############################################################################################################
##########################################################################################################
#data_cnv$emsembl <- row.names(data_cnv)
data_cnv <- merge(x = data_cnv, y = cnv_df, by.x = "emsembl", by.y = "ensembl_gene_id")
data_cnv <- as.matrix(data_cnv)

#View(data_cnv)
#dim(data_cnv)
#############################################################################################################
##########################################################################################################
#View(data_cnv)
#new_data$emsembl <- row.names(new_data)
#new_data <- merge(x = new_data, y = cnv_df, by.x = "emsembl", by.y = "ensembl_gene_id")
#new_data <- as.matrix(new_data)
#dim(new_data)
#View(new_data)
#View(data_cnv)
row.names(data_cnv) <- data_cnv[,ncol(data_cnv)]

#dim(data_cnv)
data_cnv <- data_cnv[,2:ncol(data_cnv)-2]  # Convert Ensemble ID to corresponding gene symbols
data_cnv <- matrix(as.numeric(data_cnv), nrow = nrow(data_cnv), dimnames = list(row.names(data_cnv), colnames(data_cnv)))
data_cnv <- data_cnv[!duplicated(row.names(data_cnv)),] # Only one gene PRAMEF7 has duplicate copy number variation value, and the duplicate value is the same
data_cnv<-data_cnv[,-1]
#### Keep patient samples with four genomic profiles

#View(as.data.frame(colnames(exp_assay)))
#View(as.data.frame(colnames(mty_mat)))
#length(colnames(mty_mat))
samples <- intersect(colnames(exp_assay), colnames(mty_mat)) 
#View(as.data.frame(snv_count))
samples <- intersect(samples, colnames(snv_count))
length(samples)
samples <- intersect(samples, colnames(data_cnv))
saveRDS(samples, "./data/tcga_data_processed/test1_samples.RData")


samples_test<-readRDS("./data/tcga_data_processed/test1_samples.RData")

write.csv(samples,"./data/tcga_data_processed/test1_samples_test.csv")
#### Integrate the genomic data into the network
exp_assay <- exp_assay[,samples]
mty_mat <- mty_mat[,samples]
snv_count <- snv_count[,samples]
data_cnv <- data_cnv[,samples]
####################################################################added extra information####################
graph_comp<-readRDS("./data/network_processed/graph_comp.RData")
graph_comp_updated <- upgrade_graph(graph_comp)
genes <- as.character(unlist(vertex.attributes(graph_comp_updated))) # Nodes in the maximum connected subgraph

y1 <- which(row.names(exp_assay) %in% genes)
exp_intgr <- exp_assay[y1,] 
exp_intgr <- t(exp_intgr) # Gene expression data matrix
exp_intgr <- log10(exp_intgr + 1) # 1og10 conversion
saveRDS(exp_intgr, "./data/tcga_data_processed/test1_exp_intgr.RData")
write.csv(exp_intgr,"./data/python_related/data/test1_exp_intgr.csv")

y2 <- which(row.names(mty_mat) %in% genes)
mty_intgr <- mty_mat[y2,]
mty_intgr <- t(mty_intgr) # DNA methylation data matrix
saveRDS(mty_intgr, "./data/tcga_data_processed/test1_mty_intgr.RData")
write.csv(mty_intgr,"./data/python_related/data/test1_mty_intgr.csv")

y3 <- which(row.names(snv_count) %in% genes)
snv_intgr <- snv_count[y3,]
snv_intgr <- t(snv_intgr) # Gene mutation data matrix
saveRDS(snv_intgr, "./data/tcga_data_processed/test1_snv_intgr.RData")
write.csv(snv_intgr,"./data/python_related/data/test1_snv_intgr.csv")


y4 <- which(row.names(data_cnv) %in% genes)
cnv_intgr <- data_cnv[y4,]
cnv_intgr <- t(cnv_intgr) # Copy number variation data matrix
saveRDS(cnv_intgr, "./data/tcga_data_processed/test1_cnv_intgr.RData")
write.csv(clinicalInfo,"./data/python_related/data/test1_cnv_intgr.csv")


clinicalInfo <- mty_colData
row.names(clinicalInfo) <- clinicalInfo$sample
clinicalInfo <- clinicalInfo[samples,]
saveRDS(clinicalInfo, "./data/tcga_data_processed/test1_clinical_info.RData")
write.csv(clinicalInfo,"./data/python_related/data/test1_clinical_info.csv")





save.image(file = 'my_work02_20240410.RData')





