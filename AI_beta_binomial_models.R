# conda activate R_figures
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(RColorBrewer)
library(data.table)
library(ggpubr)
library(VGAM)
library(glmmTMB)
library(broom.mixed)
library(UpSetR)
library(scales)  
# library(GWASTools)
library(ggrepel)
library(purrr)
library(stringr)

##########################################
##########################################

# 1) AI table with additional information

##########################################
##########################################

# 1.1) TCGA AI table

# alt_ratio_all <- readRDS("/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_AI_final_tables/TCGA_ASE_table.RData")
TCGA_ASE_table_merged <- readRDS("/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_AI_final_tables/TCGA_ASE_table_merged.RData")
dim(TCGA_ASE_table_merged)
# 854787    721
# Remove unnecessary columns
TCGA_ASE_table_merged <- TCGA_ASE_table_merged[,-c(51:568)]
dim(TCGA_ASE_table_merged)

# 1.2) NMD classification for nonsense mutations

# TCGA somatic PTC NMD-triggering vs NMD-evading
path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/PTCs_dataset/somatic_PTCs_ASE_all_TCGA_seq.txt"
somatic_PTCs_ASE_all_TCGA <- read.table(file = path, header  = TRUE, sep = "\t")
# Keep only those withing our ASE data
TCGA_ASE_table_merged[,c("TSS_PTC_dist","X55_nt_last_exon","PTC_CDS_exon_length")] <- NULL
ASE_stopgain <- TCGA_ASE_table_merged[TCGA_ASE_table_merged$SNV_varity == "stopgain",]
dim(ASE_stopgain)
somatic_PTCs_ASE_all_TCGA_filt <- somatic_PTCs_ASE_all_TCGA %>%
    filter(start_pos %in% ASE_stopgain$pos)
dim(somatic_PTCs_ASE_all_TCGA_filt)
# PTCs can overlap different transcripts within the same sample
# For each same PTC start position and SAME sample let's take one transcript
ensembl_transcripts_MANE_path <- "/g/strcombio/fsupek_cancer1/gpalou/NMD_project/conversor_tables/ENSEMBL/ensembl_MANE_v107.txt"
ensembl_transcripts_MANE <- read.table(file = ensembl_transcripts_MANE_path, header = TRUE)
# First the MANE Select if possible, if not, then randomly
ensembl_transcripts_MANE$MANE_main_isoform <- "yes"
somatic_PTCs_ASE_all_TCGA_filt <- merge(somatic_PTCs_ASE_all_TCGA_filt,ensembl_transcripts_MANE[,c("ensembl_transcript_id","MANE_main_isoform")],
        by.x="transcript_id", by.y ="ensembl_transcript_id", all.x = TRUE)
somatic_PTCs_ASE_all_TCGA_filt$MANE_main_isoform <- ifelse(is.na(somatic_PTCs_ASE_all_TCGA_filt$MANE_main_isoform),"no","yes")
PTC_transcripts_all_TCGA_clean <- aggregate(NMD_efficiency_TPM ~ TCGA_cancer + TCGA_barcode + 
                    PTC_CDS_pos + start_pos + MANE_main_isoform, 
                    data = somatic_PTCs_ASE_all_TCGA_filt, FUN=function(x) x[sample(seq_along(x), 1)])
# Remove those transcripts with both MANE + random (keep onle the MANE ones) ~ 664 transcripts
PTC_transcripts_all_TCGA_clean <- PTC_transcripts_all_TCGA_clean[-which(duplicated(PTC_transcripts_all_TCGA_clean[,1:4], fromLast = TRUE)),]
# Merge again
common_col <- c("TCGA_cancer","TCGA_barcode","NMD_efficiency_TPM","PTC_CDS_pos","start_pos","MANE_main_isoform")
somatic_PTCs_ASE_all_TCGA_final <- merge(PTC_transcripts_all_TCGA_clean,somatic_PTCs_ASE_all_TCGA_filt, 
                    by.x=common_col, by.y=common_col, all.x = TRUE)
somatic_PTCs_ASE_all_TCGA_final <- somatic_PTCs_ASE_all_TCGA_final[!duplicated(somatic_PTCs_ASE_all_TCGA_final),]
dim(somatic_PTCs_ASE_all_TCGA_final)
somatic_PTCs_ASE_all_TCGA_final <- somatic_PTCs_ASE_all_TCGA_final[somatic_PTCs_ASE_all_TCGA_final$MANE_main_isoform == "yes",]
# Merge with ASE
cols <- c("chr","start_pos","Ref","Alt","gene_id",
        "TCGA_barcode","TSS_PTC_dist","X55_nt_last_exon","PTC_CDS_exon_length")
somatic_PTCs_ASE_all_TCGA_final <- somatic_PTCs_ASE_all_TCGA_final[,cols]
somatic_PTCs_ASE_all_TCGA_final <- somatic_PTCs_ASE_all_TCGA_final[!duplicated(somatic_PTCs_ASE_all_TCGA_final),]
dim(somatic_PTCs_ASE_all_TCGA_final)

dim(ASE_stopgain)
ASE_stopgain <- merge(ASE_stopgain,somatic_PTCs_ASE_all_TCGA_final, all.x = TRUE,
        by.x = c("chrom","pos","ref","alt","gene_id","sample"),
        by.y = c("chr","start_pos","Ref","Alt","gene_id","TCGA_barcode"))
dim(ASE_stopgain)

dim(TCGA_ASE_table_merged)
TCGA_ASE_table <- TCGA_ASE_table_merged %>%
  filter(SNV_varity != "stopgain") %>%
  mutate(TSS_PTC_dist = NA) %>%
  mutate(X55_nt_last_exon = NA) %>%
  mutate(PTC_CDS_exon_length = NA)
TCGA_ASE_table <- rbind(TCGA_ASE_table,ASE_stopgain)
dim(TCGA_ASE_table)

TCGA_ASE_table <- TCGA_ASE_table[!duplicated(TCGA_ASE_table[, c("sample", "cancer", "chrom", "pos")]),]
dim(TCGA_ASE_table)

# 1.3) Add CNA ASCAT values

TCGA_CNA_ASCAT <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_AI_final_tables/CNA_ASCAT/TCGA_CNA_ASCAT.txt",
              header = TRUE, sep = "\t")
print(dim(TCGA_CNA_ASCAT))
# Convert to data.table
setDT(TCGA_ASE_table)
setDT(TCGA_CNA_ASCAT)
# Rename columns for easier matching
colnames(TCGA_CNA_ASCAT) <- c("GDC_Aliquot","chrom", "DNA_segment_start", "DNA_segment_end", "Copy_Number", "Major_Copy_Number", 
                            "Minor_Copy_Number", "sample", "cancer")

# Perform the range-based join
TCGA_ASE_table$SNV_pos <- TCGA_ASE_table$pos
merged_data <- TCGA_ASE_table[TCGA_CNA_ASCAT, 
                             on = .(sample, cancer, chrom, SNV_pos >= DNA_segment_start, SNV_pos <= DNA_segment_end), 
                             nomatch = 0, 
                             .(DNA_segment_start, DNA_segment_end, Copy_Number, Major_Copy_Number, Minor_Copy_Number, .SD)]
# Rename columns by removing the ".SD." prefix
setnames(merged_data, old = names(merged_data), new = gsub("^\\.SD\\.", "", names(merged_data)))
TCGA_ASE_table$SNV_pos <- NULL
merged_data$SNV_pos <- NULL

# 1.3.1) Add SNVs with no ASCAT values

setDT(TCGA_ASE_table)
setDT(merged_data)
# Select relevant columns from both data.tables
all_rows <- unique(TCGA_ASE_table[, .(sample, cancer, chrom, pos)])
merged_rows <- unique(merged_data[, .(sample, cancer, chrom, pos)])
merged_rows$new <- TRUE

# Perform a left join to keep all rows from TCGA_ASE_table
merged_result <- merge(all_rows, merged_rows, by = c("sample", "cancer", "chrom", "pos"), all.x = TRUE)
# Order as TCGA_ASE_table

# Order merged_result by sample, cancer, chrom, and pos
setorder(merged_result, sample, cancer, chrom, pos)

# Check the order of TCGA_ASE_table using `order()`
setorder(TCGA_ASE_table, sample, cancer, chrom, pos)

# Check
table(rowSums(merged_result[,c("pos","sample","cancer","chrom")] == TCGA_ASE_table[,c("pos","sample","cancer","chrom")]) == 4)

# Add the left rows to our dataset
left_ind <- which(is.na(merged_result$new))
# left_rows <- merged_result[ind,c("sample", "cancer", "chrom", "pos")]
TCGA_ASE_table_filt <- TCGA_ASE_table[left_ind,]
TCGA_ASE_table_filt[,c("DNA_segment_start","DNA_segment_end","Copy_Number","Major_Copy_Number","Minor_Copy_Number")] <- NA
TCGA_ASE_table_ASCAT <- rbind(TCGA_ASE_table_filt,merged_data)
setorder(TCGA_ASE_table_ASCAT, sample, cancer, chrom, pos)
dim(TCGA_ASE_table_ASCAT)
dim(TCGA_ASE_table)

# Check
# TCGA_ASE_table[75350,c("sample", "cancer", "chrom", "pos","RNA_ASE")]
# TCGA_ASE_table_ASCAT[75350,c("sample", "cancer", "chrom", "pos","RNA_ASE","Copy_Number","Major_Copy_Number","Minor_Copy_Number")]
# TCGA_CNA_ASCAT %>%
#   filter(chrom == "chr5",
#         cancer == "TCGA-LUSC",
#         sample == "TCGA-51-4079")

# Create WT_CNA and MUT_CNA columns based on RNA_ref and RNA_alt counts
TCGA_ASE_table_ASCAT[, `:=`(
  WT_CNA = ifelse( (DNA_ref / ( 2*(1-purity) ) ) >= (DNA_alt / (purity) ), Major_Copy_Number, Minor_Copy_Number),
  MUT_CNA = ifelse( (DNA_alt / (purity) ) > (DNA_ref / (2*(1-purity) ) ), Major_Copy_Number, Minor_Copy_Number)
)]
dim(TCGA_ASE_table_ASCAT)
# Print result
print(TCGA_ASE_table_ASCAT[50000, .(RNA_ref, RNA_alt, Major_Copy_Number, Minor_Copy_Number, WT_CNA, MUT_CNA)])
# Print cleaned column names
TCGA_ASE_table_ASCAT[1:2,]

# 1.4) Check for low coverage regions

# Define window size (adjust based on sequencing depth)
window_size <- 100000  # 100 kb

# Create windows
TCGA_ASE_table_ASCAT[, window_start := floor(pos / window_size) * window_size]
TCGA_ASE_table_ASCAT[, window_end := window_start + window_size]

# Aggregate coverage per window
window_coverage <- TCGA_ASE_table_ASCAT[, .(
  avg_coverage = mean(DNA_ref + DNA_alt, na.rm = TRUE),
  median_coverage = median(DNA_ref + DNA_alt, na.rm = TRUE),
  snv_count = .N
), by = .(chrom, sample, window_start, window_end)]

# Identify low-coverage regions (adjust threshold)
low_coverage_windows <- window_coverage[avg_coverage < 10]  # Set threshold based on expected coverage

# Print low-coverage regions
print(low_coverage_windows)

### Remove those SNVs within low-coverage regions

# Ensure TCGA_ASE_table_ASCAT and low_coverage_windows are data.tables
setDT(TCGA_ASE_table_ASCAT)
setDT(low_coverage_windows)

# Prepare low_coverage_windows for `foverlaps()` by renaming
setnames(low_coverage_windows, c("window_start", "window_end"), c("start", "end"))

# Ensure TCGA_ASE_table_ASCAT has the correct keys for range-based matching
TCGA_ASE_table_ASCAT[, start := pos]
TCGA_ASE_table_ASCAT[, end := pos]

# Set keys for foverlaps (range-based join)
setkey(TCGA_ASE_table_ASCAT, chrom, sample, start, end)
setkey(low_coverage_windows, chrom, sample, start, end)

# Perform an anti-join (keep only variants NOT in low-coverage regions)
# Find SNVs that fall within low-coverage regions
snvs_to_remove <- foverlaps(TCGA_ASE_table_ASCAT, low_coverage_windows, nomatch = 0, type = "within", which = FALSE, mult = "all")

# Ensure unique IDs for removal
snvs_to_remove <- snvs_to_remove[, .(chrom, sample, pos)]

# Perform anti-join to exclude variants in low-coverage windows
TCGA_ASE_table_ASCAT_filt <- TCGA_ASE_table_ASCAT[!snvs_to_remove, on = .(chrom, sample, pos)]

# Print cleaned dataset
print(TCGA_ASE_table_ASCAT_filt)
dim(TCGA_ASE_table_ASCAT_filt)

# 1.5) Add new purity estimates based on UCDeconvolve

UCDeconvolve_ordered <- readRDS("/g/strcombio/fsupek_cancer1/gpalou/TCGA_purity/TCGA_UCDeconvolve.RData")

TCGA_ASE_table_ASCAT_filt <- merge(TCGA_ASE_table_ASCAT_filt,UCDeconvolve_ordered[,c("TCGA_sample","hematopoietic.cell", "leukocyte","lymphocyte" ,"b.cell","t.cell","macrophage",
        "endo.epithelial.cell", "endothelial.cell","fibroblast")], 
        all.x = TRUE, by.x = "sample", by.y = "TCGA_sample")
dim(TCGA_ASE_table_ASCAT_filt)
TCGA_ASE_table_ASCAT_filt[1,]

# 1.6) Add gene-level TPM expression values

TCGA_RNAseq_TPM <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/primary_tumor/TCGA_RNAseq_matrix_TPM_gene.txt",
        header = TRUE, sep = "\t")
colnames(TCGA_RNAseq_TPM) <- gsub("\\.","-",colnames(TCGA_RNAseq_TPM))

# Ensure TCGA_RNAseq_TPM is a dataframe and store ENSG gene IDs
TCGA_RNAseq_TPM <- as.data.frame(TCGA_RNAseq_TPM)  
TCGA_RNAseq_TPM$gene_id <- rownames(TCGA_RNAseq_TPM)  # Add gene_id column

# Extract relevant gene IDs and sample barcodes from TCGA_ASE_table_ASCAT_filt
gene_list <- unique(TCGA_ASE_table_ASCAT_filt$gene_id)  # Get unique gene IDs
sample_list <- unique(substr(TCGA_ASE_table_ASCAT_filt$sample, 1, 12))  # Extract first 12 chars (short TCGA barcodes)

# Pre-filter TCGA_RNAseq_TPM to keep only relevant genes and samples
TCGA_RNAseq_TPM_filtered <- TCGA_RNAseq_TPM %>%
  filter(gene_id %in% gene_list) # Keep only genes present in TCGA_ASE_table_ASCAT_filt
TCGA_RNAseq_TPM_filtered <- TCGA_RNAseq_TPM_filtered[,colnames(TCGA_RNAseq_TPM_filtered) %in% c(sample_list,"gene_id")]
dim(TCGA_RNAseq_TPM_filtered)

# Convert TCGA_RNAseq_TPM to long format for easier merging
TCGA_RNAseq_TPM_long <- TCGA_RNAseq_TPM_filtered %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "TPM")
dim(TCGA_RNAseq_TPM_long)
head(TCGA_RNAseq_TPM_long)

# Merge TCGA_ASE_table_ASCAT_filt with filtered RNA-seq data
TCGA_ASE_table_ASCAT_filt <- TCGA_ASE_table_ASCAT_filt %>%
  left_join(TCGA_RNAseq_TPM_long, by = c("gene_id", "sample"))

# Print the first rows to confirm
dim(TCGA_ASE_table_ASCAT_filt)
TCGA_ASE_table_ASCAT_filt[1,]
cor.test(TCGA_ASE_table_ASCAT_filt$log_odds,TCGA_ASE_table_ASCAT_filt$TPM)

# 1.7) # Add imprinting genes from the paper used in PCAWG

imprinted_genes_list <- read.table(file = "/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_AI_final_tables/conversor_tables/imprinted_genes_Morison_2006.txt",
            header = FALSE, sep = "\t", stringsAsFactors = FALSE)$V1
TCGA_ASE_table_ASCAT_filt$imprinting_2 <- ifelse(TCGA_ASE_table_ASCAT_filt$Gene.Symbol %in% imprinted_genes_list, TRUE,FALSE)
dim(TCGA_ASE_table_ASCAT_filt)

# 1.8) Mutant CNA fraction + filters

# df <- TCGA_ASE_table_ASCAT_filt[,c("pos","RNA_alt","RNA_ref","RNA_ASE", "DNA_alt","DNA_ref","DNA_ASE", "TSS_PTC_dist", "X55_nt_last_exon",
#             "purity","Copy_Number","WT_CNA","MUT_CNA", "TPM","lymphocyte","fibroblast","endothelial.cell","pvalue",
#             "sample","log_odds","log_odds_cna","log_odds_total","gene","gene_cognate", "SNV_varity","SNV_type", "Gene.Symbol")]
TCGA_ASE_table_ASCAT_filt$RNA_alt <- round(TCGA_ASE_table_ASCAT_filt$RNA_alt)
TCGA_ASE_table_ASCAT_filt$RNA_ref <- round(TCGA_ASE_table_ASCAT_filt$RNA_ref)
TCGA_ASE_table_ASCAT_filt$RNA_total_reads <- TCGA_ASE_table_ASCAT_filt$RNA_alt + TCGA_ASE_table_ASCAT_filt$RNA_ref
TCGA_ASE_table_ASCAT_filt <- data.frame(TCGA_ASE_table_ASCAT_filt)
TCGA_ASE_table_ASCAT_filt$DNA_alt <- round(as.numeric(TCGA_ASE_table_ASCAT_filt$DNA_alt))
TCGA_ASE_table_ASCAT_filt$DNA_ref <- round(as.numeric(TCGA_ASE_table_ASCAT_filt$DNA_ref))
TCGA_ASE_table_ASCAT_filt$DNA_total_reads <- TCGA_ASE_table_ASCAT_filt$DNA_alt + TCGA_ASE_table_ASCAT_filt$DNA_ref

# Compute Estimated Mutant Copy Number (MUT_CNA)
TCGA_ASE_table_ASCAT_filt$estimate_MUT_CNA <- (TCGA_ASE_table_ASCAT_filt$DNA_ASE / TCGA_ASE_table_ASCAT_filt$purity) * 
                    ((TCGA_ASE_table_ASCAT_filt$purity * (TCGA_ASE_table_ASCAT_filt$Copy_Number)) + 
                    (2 * (1 - TCGA_ASE_table_ASCAT_filt$purity)))

# Copy_Number / 2 ?? --> Check paper

# Compute the fraction of mutant copy number relative to total CN, capped at 1
TCGA_ASE_table_ASCAT_filt$MUT_CNA_fraction <- pmin(TCGA_ASE_table_ASCAT_filt$estimate_MUT_CNA / (TCGA_ASE_table_ASCAT_filt$Copy_Number), 1)

TCGA_ASE_table_ASCAT_filt <- data.frame(TCGA_ASE_table_ASCAT_filt)
dim(TCGA_ASE_table_ASCAT_filt)
# cols <- c("MUT_CNA_fraction", "purity", "DNA_ref", "DNA_alt", "MUT_CNA", "WT_CNA",
#                                   "RNA_ref", "RNA_alt", "sample", "endothelial.cell",
#                                   "TPM", "fibroblast", "lymphocyte")
# na_counts <- colSums(is.na(df[, cols]))
# print(na_counts)

# rem_rows <- which(rowSums(is.na(df[,cols])) > 0)
# df <- df[-rem_rows,]
# dim(df)

saveRDS(TCGA_ASE_table_ASCAT_filt, "/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_AI_final_tables/TCGA_AI_table_bb_models_tmp1.RData")
TCGA_ASE_table_ASCAT_filt <- readRDS("/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_AI_final_tables/TCGA_AI_table_bb_models_tmp1.RData")
dim(TCGA_ASE_table_ASCAT_filt)

##############################################
##############################################
# 2) RNA-AI modelling with beta-binomial model
##############################################
##############################################

# 2.1) RNA Beta-Binomial model

# Stage 1: Adjust RNA counts for CNA and impurity effects using beta-binomial

RNA_bb_model1 <- VGAM::vglm(
  cbind(RNA_alt, RNA_ref) ~ MUT_CNA_fraction + sqrt(TPM),
  family = VGAM::betabinomial,
  na.action = na.exclude,
  data = TCGA_ASE_table_ASCAT_filt
)

RNA_bb_model2 <- VGAM::vglm(
  cbind(RNA_alt, RNA_ref) ~ MUT_CNA_fraction + lymphocyte + fibroblast + endothelial.cell + 
        sqrt(TPM) + lymphocyte:sqrt(TPM) + fibroblast:sqrt(TPM) + endothelial.cell:sqrt(TPM),
  family = VGAM::betabinomial,
  na.action = na.exclude,
  data = TCGA_ASE_table_ASCAT_filt
)

RNA_bb_model3 <- VGAM::vglm(
  cbind(RNA_alt, RNA_ref) ~ Copy_Number + purity + Copy_Number:purity + sqrt(TPM) + purity:sqrt(TPM),
  family = VGAM::betabinomial,
  na.action = na.exclude,
  data = TCGA_ASE_table_ASCAT_filt
)

RNA_bb_model4 <- VGAM::vglm(
  cbind(RNA_alt, RNA_ref) ~ MUT_CNA + WT_CNA + MUT_CNA:WT_CNA + purity + MUT_CNA:purity + 
              WT_CNA:purity + sqrt(TPM) + purity:sqrt(TPM),
  family = VGAM::betabinomial,
  na.action = na.exclude,
  data = TCGA_ASE_table_ASCAT_filt
)

RNA_null_model <- VGAM::vglm(
  cbind(RNA_alt, RNA_ref) ~ purity + sqrt(TPM) + purity:sqrt(TPM),
  family = VGAM::betabinomial,
  na.action = na.exclude,
  data = TCGA_ASE_table_ASCAT_filt
)

# Save
saveRDS(RNA_bb_model1, "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/RNA_AI_bb_model1.RData")
saveRDS(RNA_bb_model2, "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/RNA_AI_bb_model2.RData")
saveRDS(RNA_bb_model3, "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/RNA_AI_bb_model3.RData")
saveRDS(RNA_bb_model4, "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/RNA_AI_bb_model4.RData")
# saveRDS(RNA_bb_model5, "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/RNA_AI_bb_model5.RData")
saveRDS(RNA_null_model, "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/RNA_AI_null_model.RData")

RNA_bb_model1 <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/RNA_AI_bb_model1.RData")
RNA_bb_model2 <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/RNA_AI_bb_model2.RData")
RNA_bb_model3 <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/RNA_AI_bb_model3.RData")
RNA_bb_model4 <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/RNA_AI_bb_model4.RData")
# RNA_bb_model5 <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/RNA_AI_bb_model5.RData")
RNA_null_model <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/RNA_AI_null_model.RData")

# 2.2) DNA modelling with beta-binomial model

DNA_bb_model1 <- VGAM::vglm(
  cbind(DNA_alt, DNA_ref) ~ MUT_CNA_fraction,
  family = VGAM::betabinomial,
  na.action = na.exclude,
  data = TCGA_ASE_table_ASCAT_filt
)

DNA_bb_model2 <- VGAM::vglm(
  cbind(DNA_alt, DNA_ref) ~ Copy_Number + purity + Copy_Number:purity,
  family = VGAM::betabinomial,
  na.action = na.exclude,
  data = TCGA_ASE_table_ASCAT_filt
)

DNA_bb_model3 <- VGAM::vglm(
  cbind(DNA_alt, DNA_ref) ~ MUT_CNA + purity + MUT_CNA:purity,
  family = VGAM::betabinomial,
  na.action = na.exclude,
  data = TCGA_ASE_table_ASCAT_filt
)

DNA_bb_model4 <- VGAM::vglm(
  cbind(DNA_alt, DNA_ref) ~ MUT_CNA + WT_CNA + MUT_CNA*WT_CNA + purity,
  family = VGAM::betabinomial,
  na.action = na.exclude,
  data = TCGA_ASE_table_ASCAT_filt
)

DNA_null_model <- VGAM::vglm(
  cbind(DNA_alt, DNA_ref) ~ purity,  # Removing CNA
  family = VGAM::betabinomial,
  na.action = na.exclude,
  data = TCGA_ASE_table_ASCAT_filt
)

saveRDS(DNA_bb_model1, "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/DNA_AI_bb_model1.RData")
saveRDS(DNA_bb_model2, "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/DNA_AI_bb_model2.RData")
saveRDS(DNA_bb_model3, "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/DNA_AI_bb_model3.RData")
saveRDS(DNA_bb_model4, "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/DNA_AI_bb_model4.RData")
saveRDS(DNA_null_model, "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/DNA_AI_bb_null_model.RData")

DNA_bb_model1 <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/DNA_AI_bb_model1.RData")
DNA_bb_model2 <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/DNA_AI_bb_model2.RData")
DNA_bb_model3 <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/DNA_AI_bb_model3.RData")
DNA_bb_model4 <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/DNA_AI_bb_model4.RData")
DNA_null_model <- readRDS("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/DNA_AI_bb_null_model.RData")

# 2.3) Significant p-values

TCGA_ASE_table_ASCAT_filt <- readRDS("/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_AI_final_tables/TCGA_AI_table_bb_models_tmp1.RData")
dim(TCGA_ASE_table_ASCAT_filt)

analyze_allelic_imbalance <- function(df_merged, model, model_type = "DNA_bb_model1") {
  # Extract prefix (DNA/RNA) and model number from model_type
  parts <- strsplit(model_type, "_")[[1]]
  prefix <- parts[1]  # DNA or RNA
  model_num <- gsub("model", "", parts[length(parts)])  # Extract model number
  
  # Generate column suffix based on model type
  suffix <- paste0("_", parts[2], "_model", model_num)  # e.g., "_bb_model1"
  
  # Make a copy of the input dataframe to avoid modifying the original
  results_df <- df_merged
  
  # Get expected proportions from the model
  expected_prop <- fitted(model)[, 1]
  results_df[[paste0(prefix, "_expected_ratio", suffix)]] <- expected_prop
  
  # Calculate observed ratios if not already present
  ratio_col <- paste0(prefix, "_observed_ratio")
  if (!ratio_col %in% names(results_df)) {
    results_df[[ratio_col]] <- results_df[[paste0(prefix, "_alt")]] / 
                               (results_df[[paste0(prefix, "_alt")]] + results_df[[paste0(prefix, "_ref")]])
  }
  
  # Get dispersion parameter from model
  model_coefs <- coef(model)
  # Use plogis for dispersion as it seems to be the approach you're using
  model_dispersion <- plogis(model_coefs["(Intercept):2"])
  
  # Calculate p-values for model (both tails)
  alt_counts <- results_df[[paste0(prefix, "_alt")]]
  total_counts <- results_df[[paste0(prefix, "_alt")]] + results_df[[paste0(prefix, "_ref")]]
  
  # Lower tail: observed is lower than expected (favors reference)
  results_df[[paste0(prefix, "_p_lower", suffix)]] <- VGAM::pbetabinom(
    alt_counts, size = total_counts, prob = expected_prop, rho = model_dispersion
  )
  
  # Upper tail: observed is higher than expected (favors alternate)
  results_df[[paste0(prefix, "_p_upper", suffix)]] <- 1 - VGAM::pbetabinom(
    alt_counts - 1, size = total_counts, prob = expected_prop, rho = model_dispersion
  )

  # FDR adjust
  results_df[[paste0(prefix, "_FDR_lower", suffix)]] <- p.adjust(results_df[[paste0(prefix, "_p_lower", suffix)]], method = "fdr")
  results_df[[paste0(prefix, "_FDR_upper", suffix)]] <- p.adjust(results_df[[paste0(prefix, "_p_upper", suffix)]], method = "fdr")
  
  # Return the augmented dataframe
  return(results_df)
}

analyze_all_allelic_imbalance_models <- function(df_merged) {
  # Make a copy of the input dataframe
  results_df <- df_merged
  
  # List of all DNA models with their model_type identifiers
  dna_models <- list(
    list(model = DNA_bb_model1, type = "DNA_bb_model1"),
    list(model = DNA_bb_model2, type = "DNA_bb_model2"),
    list(model = DNA_bb_model3, type = "DNA_bb_model3"),
    list(model = DNA_bb_model4, type = "DNA_bb_model4"),
    list(model = DNA_null_model, type = "DNA_null_model")
  )
  
  # List of all RNA models with their model_type identifiers
  rna_models <- list(
    list(model = RNA_bb_model1, type = "RNA_bb_model1"),
    list(model = RNA_bb_model2, type = "RNA_bb_model2"),
    list(model = RNA_bb_model3, type = "RNA_bb_model3"),
    list(model = RNA_bb_model4, type = "RNA_bb_model4"),
    list(model = RNA_null_model, type = "RNA_null_model")
  )
  
  # Process all DNA models
  cat("Processing DNA models...\n")
  for (model_info in dna_models) {
    cat("  Processing", model_info$type, "...\n")
    results_df <- analyze_allelic_imbalance(
      results_df, 
      model = model_info$model,
      model_type = model_info$type
    )
  }
  
  # Process all RNA models
  cat("Processing RNA models...\n")
  for (model_info in rna_models) {
    cat("  Processing", model_info$type, "...\n")
    results_df <- analyze_allelic_imbalance(
      results_df, 
      model = model_info$model,
      model_type = model_info$type
    )
  }
  
  # Return the dataframe with all model results
  return(results_df)
}

# Process all models and get comprehensive results
TCGA_ASE_table_bb_models <- analyze_all_allelic_imbalance_models(df_merged = TCGA_ASE_table_ASCAT_filt)

dim(TCGA_ASE_table_ASCAT_filt)
dim(TCGA_ASE_table_bb_models)

saveRDS(TCGA_ASE_table_bb_models, "/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_AI_final_tables/TCGA_AI_table_bb_models_tmp2.RData")

##############################################
##############################################
# 3) Plots and analysis for quality control
##############################################
##############################################

TCGA_ASE_table_bb_models <- readRDS("/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_AI_final_tables/TCGA_AI_table_bb_models_tmp2.RData")

# 3.1) Observed vs Expected RNA/DNA-AI for each bb-model, including old method

# Plot RNA models including log_odds model with percentages of significance
plot_RNA_models <- function(df) {
  # Get all RNA model columns for expected ratio
  rna_models <- grep("RNA_expected_ratio_bb_model", names(df), value = TRUE)
  # Add null model if it exists
  rna_null <- grep("RNA_expected_ratio_null", names(df), value = TRUE)
  rna_models <- c(rna_models, rna_null)
  
  # Create a long-format dataframe for plotting
  plot_data <- data.frame(
    SNV_ID = rep(1:nrow(df), length(rna_models) + 1),  # +1 for log_odds
    observed_ratio = rep(df$RNA_observed_ratio, length(rna_models) + 1),
    model = rep(NA, nrow(df) * (length(rna_models) + 1)),
    expected_ratio = rep(NA, nrow(df) * (length(rna_models) + 1)),
    significant = rep(FALSE, nrow(df) * (length(rna_models) + 1)),
    sig_type = rep(NA, nrow(df) * (length(rna_models) + 1))  # Add significance type
  )
  
  # Store significance percentages for annotation
  sig_stats <- data.frame(
    model = character(),
    lower_pct = numeric(),
    upper_pct = numeric(),
    total_sig_pct = numeric(),
    log_neg_pct = numeric(),
    log_pos_pct = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Fill in the data for each beta-binomial model
  row_idx <- 1
  for (model_col in rna_models) {
    # Extract model name from column name
    model_name <- gsub("RNA_expected_ratio_", "", model_col)
    
    # Fill in this model's data
    idx <- row_idx:(row_idx + nrow(df) - 1)
    plot_data$model[idx] <- model_name
    plot_data$expected_ratio[idx] <- df[[model_col]]
    
    # Determine significance based on p-values and categorize by type
    p_lower_col <- paste0("RNA_p_lower_", model_name)
    p_upper_col <- paste0("RNA_p_upper_", model_name)
    
    if (p_lower_col %in% names(df) && p_upper_col %in% names(df)) {
      sig_lower <- !is.na(df[[p_lower_col]]) & (df[[p_lower_col]] < 0.05)
      sig_upper <- !is.na(df[[p_upper_col]]) & (df[[p_upper_col]] < 0.05)
      
      plot_data$significant[idx] <- sig_lower | sig_upper
      
      # Categorize significance type
      plot_data$sig_type[idx] <- "Not Significant"
      plot_data$sig_type[idx][sig_lower & !sig_upper] <- "Lower"
      plot_data$sig_type[idx][!sig_lower & sig_upper] <- "Upper"
      
      # Calculate percentages for annotation
      n_valid <- sum(!is.na(df[[p_lower_col]]) & !is.na(df[[p_upper_col]]))
      
      if (n_valid > 0) {
        lower_pct <- sum(sig_lower & !sig_upper) / n_valid * 100
        upper_pct <- sum(!sig_lower & sig_upper) / n_valid * 100
        total_sig_pct <- sum(sig_lower | sig_upper) / n_valid * 100
        
        sig_stats <- rbind(sig_stats, data.frame(
          model = model_name,
          lower_pct = lower_pct,
          upper_pct = upper_pct,
          total_sig_pct = total_sig_pct,
          log_neg_pct = NA,
          log_pos_pct = NA
        ))
      }
    }
    
    row_idx <- row_idx + nrow(df)
  }
  
  # Add log_odds model
  if ("log_odds" %in% names(df)) {
    idx <- row_idx:(row_idx + nrow(df) - 1)
    plot_data$model[idx] <- "log_odds"
    # For log_odds, we'll use observed ratio for x-axis (since we don't have an expected value)
    plot_data$expected_ratio[idx] <- df$log_odds  # neutral expectation
    
    # Determine significance based on log_odds and categorize
    log_neg <- !is.na(df$log_odds) & (df$log_odds < -1)
    log_pos <- !is.na(df$log_odds) & (df$log_odds > 1)
    
    plot_data$significant[idx] <- log_neg | log_pos
    
    # Categorize log_odds significance
    plot_data$sig_type[idx] <- "Not Significant"
    plot_data$sig_type[idx][log_neg] <- "Negative"
    plot_data$sig_type[idx][log_pos] <- "Positive"
    
    # Calculate percentages for log_odds
    n_valid <- sum(!is.na(df$log_odds))
    
    if (n_valid > 0) {
      log_neg_pct <- sum(log_neg) / n_valid * 100
      log_pos_pct <- sum(log_pos) / n_valid * 100
      total_log_pct <- sum(log_neg | log_pos) / n_valid * 100
      
      sig_stats <- rbind(sig_stats, data.frame(
        model = "log_odds",
        lower_pct = NA,
        upper_pct = NA,
        total_sig_pct = total_log_pct,
        log_neg_pct = log_neg_pct,
        log_pos_pct = log_pos_pct
      ))
    }
  }
  
  # Create annotation text for each model
  annotations <- lapply(1:nrow(sig_stats), function(i) {
    model <- sig_stats$model[i]
    
    if (model != "log_odds") {
      # For beta-binomial models
      sprintf(
        "%s\nLower: %.1f%%\nUpper: %.1f%%\nTotal: %.1f%%",
        model,
        sig_stats$lower_pct[i],
        sig_stats$upper_pct[i],
        sig_stats$total_sig_pct[i]
      )
    } else {
      # For log_odds model
      sprintf(
        "%s\nNeg: %.1f%%\nPos: %.1f%%\nTotal: %.1f%%",
        model,
        sig_stats$log_neg_pct[i],
        sig_stats$log_pos_pct[i],
        sig_stats$total_sig_pct[i]
      )
    }
  })
  names(annotations) <- sig_stats$model
  
  # Create plot with custom significance colors
  ggplot(plot_data, aes(x = expected_ratio, y = observed_ratio, color = sig_type)) +
    geom_point(alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    facet_wrap(~ model, scales = "free") +
    scale_color_manual(values = c(
      "Not Significant" = "grey60", 
      "Lower" = "blue", 
      "Upper" = "red", 
      "Negative" = "blue",
      "Positive" = "red"
    )) +
    # Add annotations with percentages
    geom_text(
      data = data.frame(
        model = sig_stats$model,
        x = rep(0.8, nrow(sig_stats)),  # Position for annotation
        y = rep(0.2, nrow(sig_stats)),  # Position for annotation
        label = unlist(annotations)
      ),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 0,
      size = 3
    ) +
    labs(x = "RNA-AI Expected", 
         y = "RNA-AI Observed",
         title = "RNA-AI Models",
         subtitle = "Blue: Lower tail sig, Red: Upper tail sig. (p<0.05)",
         color = "Significance") +
    theme_classic() +
    theme(legend.position = "bottom")
}

# Plot DNA models with more complex criteria, including log_odds model
plot_DNA_models <- function(df) {
  # Get all DNA model columns for expected ratio
  dna_models <- grep("DNA_expected_ratio_bb_model", names(df), value = TRUE)
  # Add null model if it exists
  dna_null <- grep("DNA_expected_ratio_null", names(df), value = TRUE)
  dna_models <- c(dna_models, dna_null)
  
  # Create a long-format dataframe for plotting
  plot_data <- data.frame(
    SNV_ID = rep(1:nrow(df), length(dna_models) + 1),  # +1 for log_odds
    observed_ratio = rep(df$DNA_observed_ratio, length(dna_models) + 1),
    model = rep(NA, nrow(df) * (length(dna_models) + 1)),
    expected_ratio = rep(NA, nrow(df) * (length(dna_models) + 1)),
    significant = rep(FALSE, nrow(df) * (length(dna_models) + 1)),
    sig_type = rep(NA, nrow(df) * (length(dna_models) + 1))  # Add significance type
  )
  
  # Store significance percentages for annotation
  sig_stats <- data.frame(
    model = character(),
    lower_pct = numeric(),
    upper_pct = numeric(),
    sig_null_not_full_pct = numeric(),
    log_neg_pct = numeric(),
    log_pos_pct = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Fill in the data for each beta-binomial model
  row_idx <- 1
  for (model_col in dna_models) {
    # Extract model name from column name
    model_name <- gsub("DNA_expected_ratio_", "", model_col)
    
    # Fill in this model's data
    idx <- row_idx:(row_idx + nrow(df) - 1)
    plot_data$model[idx] <- model_name
    plot_data$expected_ratio[idx] <- df[[model_col]]
    
    # Determine significance based on p-values
    p_lower_col <- paste0("DNA_p_lower_", model_name)
    p_upper_col <- paste0("DNA_p_upper_", model_name)
    
    if (p_lower_col %in% names(df) && p_upper_col %in% names(df)) {
      sig_lower <- !is.na(df[[p_lower_col]]) & (df[[p_lower_col]] < 0.05)
      sig_upper <- !is.na(df[[p_upper_col]]) & (df[[p_upper_col]] < 0.05)
      sig_any <- sig_lower | sig_upper
      
      # Basic significance for coloring points
      plot_data$sig_type[idx] <- "Not Significant"
      plot_data$sig_type[idx][sig_lower & !sig_upper] <- "Lower"
      plot_data$sig_type[idx][!sig_lower & sig_upper] <- "Upper"
      
      # Calculate percentages for basic significance
      n_valid <- sum(!is.na(df[[p_lower_col]]) & !is.na(df[[p_upper_col]]))
      
      if (n_valid > 0) {
        lower_pct <- sum(sig_lower & !sig_upper) / n_valid * 100
        upper_pct <- sum(!sig_lower & sig_upper) / n_valid * 100
        
        # For the special "significant in null but not in full model" criteria
        if (model_name != "null_model" && "DNA_p_lower_null_model" %in% names(df) && "DNA_p_upper_null_model" %in% names(df)) {
          # Logic: NOT significant in full model (p_lower >= 0.05 AND p_upper >= 0.05)
          # AND significant in null model (p_lower_null < 0.05 OR p_upper_null < 0.05)
          p_lower_null <- df$DNA_p_lower_null_model
          p_upper_null <- df$DNA_p_upper_null_model
          
          not_sig_full <- !is.na(df[[p_lower_col]]) & !is.na(df[[p_upper_col]]) &
                          (df[[p_lower_col]] >= 0.05 & df[[p_upper_col]] >= 0.05)
                          
          sig_null <- !is.na(p_lower_null) & !is.na(p_upper_null) &
                    (p_lower_null < 0.05 | p_upper_null < 0.05)
                    
          special_sig <- not_sig_full & sig_null
          
          # Update the significance field for this special case
          plot_data$significant[idx] <- special_sig
          
          # Calculate percentage for special significance
          sig_null_not_full_pct <- sum(special_sig, na.rm = TRUE) / n_valid * 100
        } else {
          # Just use regular significance
          plot_data$significant[idx] <- sig_any
          sig_null_not_full_pct <- NA
        }
        
        sig_stats <- rbind(sig_stats, data.frame(
          model = model_name,
          lower_pct = lower_pct,
          upper_pct = upper_pct,
          sig_null_not_full_pct = sig_null_not_full_pct,
          log_neg_pct = NA,
          log_pos_pct = NA
        ))
      }
    }
    
    row_idx <- row_idx + nrow(df)
  }
  
  # Add log_odds model
  if ("log_odds_cna" %in% names(df)) {
    idx <- row_idx:(row_idx + nrow(df) - 1)
    plot_data$model[idx] <- "log_odds_cna"
    # For log_odds, we'll use observed ratio for x-axis (since we don't have an expected value)
    plot_data$expected_ratio[idx] <- df$log_odds_cna  # neutral expectation
    
    # Determine significance based on log_odds and categorize
    log_neg <- !is.na(df$log_odds_cna) & (df$log_odds_cna < -1)
    log_pos <- !is.na(df$log_odds_cna) & (df$log_odds_cna > 1)
    
    plot_data$significant[idx] <- log_neg | log_pos
    
    # Categorize log_odds significance
    plot_data$sig_type[idx] <- "Not Significant"
    plot_data$sig_type[idx][log_neg] <- "Negative"
    plot_data$sig_type[idx][log_pos] <- "Positive"
    
    # Calculate percentages for log_odds
    n_valid <- sum(!is.na(df$log_odds_cna))
    
    if (n_valid > 0) {
      log_neg_pct <- sum(log_neg) / n_valid * 100
      log_pos_pct <- sum(log_pos) / n_valid * 100
      
      sig_stats <- rbind(sig_stats, data.frame(
        model = "log_odds_cna",
        lower_pct = NA,
        upper_pct = NA,
        sig_null_not_full_pct = NA,
        log_neg_pct = log_neg_pct,
        log_pos_pct = log_pos_pct
      ))
    }
  }
  
  # Create annotation text for each model
  annotations <- lapply(1:nrow(sig_stats), function(i) {
    model <- sig_stats$model[i]
    
    if (model == "null_model") {
      # For null model
      sprintf(
        "%s\nLower: %.1f%%\nUpper: %.1f%%",
        model,
        sig_stats$lower_pct[i],
        sig_stats$upper_pct[i]
      )
    } else if (model == "log_odds") {
      # For log_odds model
      sprintf(
        "%s\nNeg: %.1f%%\nPos: %.1f%%",
        model,
        sig_stats$log_neg_pct[i],
        sig_stats$log_pos_pct[i]
      )
    } else {
      # For other beta-binomial models
      sprintf(
        "%s\nLower: %.1f%%\nUpper: %.1f%%\nSig in Null\nnot in Full: %.1f%%",
        model,
        sig_stats$lower_pct[i],
        sig_stats$upper_pct[i],
        sig_stats$sig_null_not_full_pct[i]
      )
    }
  })
  names(annotations) <- sig_stats$model
  
  # Create plot
  ggplot(plot_data, aes(x = expected_ratio, y = observed_ratio)) +
    # Plot all points in gray first
    geom_point(alpha = 0.5, color = "grey80") +
    # Then plot colored points by significance type
    geom_point(
      data = subset(plot_data, sig_type != "Not Significant"),
      aes(color = sig_type),
      alpha = 0.7
    ) +
    # Special case for "significant in null but not in full" criteria
    geom_point(
      data = subset(plot_data, significant & model != "log_odds_cna" & model != "null_model"),
      color = "green",
      shape = 1,
      size = 3,
      alpha = 0.8
    ) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    facet_wrap(~ model, scales = "free") +
    scale_color_manual(values = c(
      "Not Significant" = "grey60", 
      "Lower" = "blue", 
      "Upper" = "red", 
      "Negative" = "blue",
      "Positive" = "red"
    )) +
    # Add annotations with percentages
    geom_text(
      data = data.frame(
        model = sig_stats$model,
        x = rep(0.8, nrow(sig_stats)),  # Position for annotation
        y = rep(0.2, nrow(sig_stats)),  # Position for annotation
        label = unlist(annotations)
      ),
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 0,
      size = 3
    ) +
    labs(x = "DNA-AI Expected", 
         y = "DNA-AI Observed",
         title = "DNA-AI Models",
         subtitle = "Blue: Lower tail, Red: Upper tail, Green circles: Significant in Null but not in Full model",
         color = "Significance Type") +
    theme_classic() +
    theme(legend.position = "bottom")
}

# Create plots for RNA models
rna_plot <- plot_RNA_models(TCGA_ASE_table_bb_models)
print(rna_plot)

# Create plots for DNA models
dna_plot <- plot_DNA_models(TCGA_ASE_table_bb_models)
print(dna_plot)

# Save the plots if needed
ggsave("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/bb_models/obs_vs_pred_RNA_AI.png", 
       rna_plot, width = 200, height = 180, units = "mm")
# Save the plots if needed
ggsave("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/bb_models/obs_vs_pred_DNA_AI.png", 
       dna_plot, width = 200, height = 180, units = "mm")

# 3.2) Absolute number of significant hits for bb-models & old method

# A) RNA-AI

# Define significance criteria for all models
df <- TCGA_ASE_table_bb_models %>%
  dplyr::mutate(
    # Original criteria
    log_odds_sig = abs(log_odds) > 1,
    # RNA beta-binomial model significance criteria
    RNA_sig_bb_model1 = RNA_p_upper_bb_model1 < 0.05 | RNA_p_lower_bb_model1 < 0.05,
    RNA_sig_bb_model2 = RNA_p_upper_bb_model2 < 0.05 | RNA_p_lower_bb_model2 < 0.05,
    RNA_sig_bb_model3 = RNA_p_upper_bb_model3 < 0.05 | RNA_p_lower_bb_model3 < 0.05,
    RNA_sig_bb_model4 = RNA_p_upper_bb_model4 < 0.05 | RNA_p_lower_bb_model4 < 0.05,
    RNA_sig_null_model = RNA_p_upper_null_model < 0.05 | RNA_p_lower_null_model < 0.05
  )

# Select only significant calls for overlap analysis
for (gene_type in c("random", "og", "tsg", "TP53")) {
  print(gene_type)

  if (gene_type == "random") {
    df_filt <- df %>% dplyr::filter(gene_cognate == "random")
  } else if (gene_type == "og") {
    df_filt <- df %>% dplyr::filter(gene_cognate == "og_cognate")
  } else if (gene_type == "tsg") {
    df_filt <- df %>% dplyr::filter(gene_cognate %in% c("tsg_cognate"))
  } else if (gene_type == "TP53") {
    df_filt <- df %>% dplyr::filter(gene_cognate %in% c("TP53"))
  }

  # Select all model significance columns for UpSet plot
  sig_df <- df_filt %>%
    dplyr::select(
      RNA_sig_bb_model1,
      RNA_sig_bb_model2,
      RNA_sig_bb_model3,
      RNA_sig_bb_model4,
      RNA_sig_null_model,
      log_odds_sig
    )

  # Convert TRUE/FALSE to binary format for UpSet plot
  sig_df <- sig_df %>%
    dplyr::mutate_all(~as.integer(.))

  print(dim(sig_df))

  sig_df <- na.omit(sig_df)

  # Define output path
  final_figure_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/bb_models/",gene_type, "_prop_sig_RNA_AI_by_bb_models.png")

  # Save the UpSet plot
  png(final_figure_path, width = 300, height = 200, units = "mm", res = 300)  # Increased width for more sets
  plot <- upset(
    sig_df,
    sets = c(
      "RNA_sig_bb_model1",
      "RNA_sig_bb_model2",
      "RNA_sig_bb_model3",
      "RNA_sig_bb_model4",
      "RNA_sig_null_model",
      "log_odds_sig"
    ),
    order.by = "freq",
    keep.order = TRUE,
    mainbar.y.label = "Significant SNV Overlaps",
    sets.x.label = "Number of Significant SNVs per Method",
    sets.bar.color = "darkblue",
    main.bar.color = "darkred",
    text.scale = 1.2  # Slightly larger text for better readability
  )
  print(plot)
  dev.off()  # Close the plotting device
}

# B) DNA

# Define significance criteria for all DNA models
df <- TCGA_ASE_table_bb_models %>%
  dplyr::mutate(
    # Original criteria
    log_odds_cna_sig = abs(log_odds_cna) > 1,
    
    # DNA beta-binomial model significance criteria - significance defined as:
    # (p_upper OR p_lower < 0.05 in bb model) AND (p_upper OR p_lower < 0.05 in null model)
    # DNA_sig_bb_model1 = ( DNA_p_upper_bb_model1 >= 0.05 & DNA_p_upper_null_model < 0.05 ) |
    #                      ( DNA_p_lower_bb_model1 >= 0.05 & DNA_p_lower_null_model < 0.05),
    # DNA_sig_bb_model2 = ( DNA_p_upper_bb_model2 >= 0.05 & DNA_p_upper_null_model < 0.05 ) |
    #                      ( DNA_p_lower_bb_model2 >= 0.05 & DNA_p_lower_null_model < 0.05),
    # DNA_sig_bb_model3 = ( DNA_p_upper_bb_model3 >= 0.05 & DNA_p_upper_null_model < 0.05 ) |
    #                      ( DNA_p_lower_bb_model3 >= 0.05 & DNA_p_lower_null_model < 0.05),
    # DNA_sig_bb_model4 = ( DNA_p_upper_bb_model4 >= 0.05 & DNA_p_upper_null_model < 0.05 ) |
    #                      ( DNA_p_lower_bb_model4 >= 0.05 & DNA_p_lower_null_model < 0.05),
    DNA_sig_bb_model1 = ( DNA_p_upper_bb_model1 < 0.05) |
                         ( DNA_p_lower_bb_model1 < 0.05),
    DNA_sig_bb_model2 = ( DNA_p_upper_bb_model2 < 0.05) |
                         ( DNA_p_lower_bb_model2 < 0.05),
    DNA_sig_bb_model3 = ( DNA_p_upper_bb_model3 < 0.05) |
                         ( DNA_p_lower_bb_model3 < 0.05),
    DNA_sig_bb_model4 = ( DNA_p_upper_bb_model4 < 0.05) |
                         ( DNA_p_lower_bb_model4 < 0.05),
    # Null model by itself - significance defined as p_upper OR p_lower < 0.05
    DNA_sig_null_model = DNA_p_upper_null_model < 0.05 | DNA_p_lower_null_model < 0.05
  )

# Select only significant calls for overlap analysis
for (gene_type in c("random", "og", "tsg", "TP53")) {
  print(gene_type)

  if (gene_type == "random") {
    df_filt <- df %>% dplyr::filter(gene_cognate == "random")
  } else if (gene_type == "og") {
    df_filt <- df %>% dplyr::filter(gene_cognate == "og_cognate")
  } else if (gene_type == "tsg") {
    df_filt <- df %>% dplyr::filter(gene_cognate %in% c("tsg_cognate"))
  } else if (gene_type == "TP53") {
    df_filt <- df %>% dplyr::filter(gene_cognate %in% c("TP53"))
  }

  # Select all model significance columns for UpSet plot
  sig_df <- df_filt %>%
    dplyr::select(
      DNA_sig_bb_model1,
      DNA_sig_bb_model2,
      DNA_sig_bb_model3,
      DNA_sig_bb_model4,
      DNA_sig_null_model,
      log_odds_cna_sig
    )

  # Convert TRUE/FALSE to binary format for UpSet plot
  sig_df <- sig_df %>%
    dplyr::mutate_all(~as.integer(.))

  sig_df <- na.omit(sig_df)

  print(dim(sig_df))

  # Define output path
  final_figure_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/bb_models/", gene_type, "_prop_sig_DNA_AI_by_bb_models.png")

  # Save the UpSet plot
  png(final_figure_path, width = 300, height = 200, units = "mm", res = 300)  # Increased width for more sets
  plot <- upset(
    sig_df,
    sets = c(
      "DNA_sig_bb_model1",
      "DNA_sig_bb_model2",
      "DNA_sig_bb_model3",
      "DNA_sig_bb_model4",
      "DNA_sig_null_model",
      "log_odds_cna_sig"
    ),
    order.by = "freq",
    keep.order = TRUE,
    mainbar.y.label = "Significant SNV Overlaps",
    sets.x.label = "Number of Significant SNVs per Method",
    sets.bar.color = "darkblue",
    main.bar.color = "darkred",
    text.scale = 1.2  # Slightly larger text for better readability
  )
  print(plot)
  dev.off()  # Close the plotting device
}

# 3.3) % DNA-AI vs RNA-AI

df <- TCGA_ASE_table_bb_models %>%
  # a) restrict to mRNA‐AI (null p < 0.05)

  # b) compute CNA‐AI & reg‐AI for this pairing
  mutate(
    mRNA = ( RNA_p_upper_null_model < 0.05 | RNA_p_lower_null_model < 0.05 ),
    CNA = (DNA_p_upper_null_model < 0.05 &
            DNA_p_upper_bb_model1 >= 0.05) |
          (DNA_p_lower_null_model < 0.05 &
            DNA_p_lower_bb_model1 >= 0.05),
    
    reg = ( RNA_p_upper_bb_model1 < 0.05 &
            DNA_p_upper_null_model >= 0.05) |
          ( RNA_p_lower_bb_model1 < 0.05 &
            DNA_p_lower_null_model >= 0.05)
  )

prop.table(table(df$mRNA) )
prop.table(table(df$reg) )
table(df$reg)
prop.table(table(df$CNA) )
table(df$CNA)
table(abs(df$log_odds) > 1)
prop.table(table(abs(df$log_odds) > 1))

# A) % of significant reg vs CNA-AI from significant mRNA-AI SNVs for each BB-model --> All gene_cognate

# Define your five RNA–DNA model pairings
model_pairs <- tibble::tribble(
  ~model_pair,   ~rna_model, ~dna_model,
  "RNA1-DNA1",   "1",        "1",
  "RNA2-DNA1",   "2",        "1",
  "RNA3-DNA2",   "3",        "2",
  "RNA4-DNA3",   "4",        "3",
  "RNA4-DNA4",   "4",        "4"
)

# Define your gene-type filters
gene_types <- c("random", "og", "tsg", "TP53")
gene_filters <- list(
  random = expr(gene_cognate == "random"),
  og     = expr(gene_cognate == "og_cognate"),
  tsg    = expr(gene_cognate %in% c("tsg_cognate")),
  TP53   = expr(gene_cognate == "TP53")
)

# Build one big results table
results_all <- imap_dfr(gene_filters, function(filter_expr, gene_type) {
  df_filt <- TCGA_ASE_table_bb_models %>%
    filter(!!filter_expr) #%>%
    # filter(SNV_varity == "stopgain")
  
  # skip if no SNVs for this gene set
  if(nrow(df_filt) == 0) return(tibble())
  
  model_pairs %>%
    pmap_dfr(function(model_pair, rna_model, dna_model) {
      df_filt %>%
        # a) restrict to mRNA‐AI (null p < 0.05)
        dplyr::filter(RNA_p_upper_null_model < 0.05 | RNA_p_lower_null_model < 0.05) %>%
        
        # b) compute CNA‐AI & reg‐AI for this pairing
        dplyr::mutate(
          CNA = (DNA_p_upper_null_model < 0.05 &
                 !!sym(paste0("DNA_p_upper_bb_model", dna_model)) >= 0.05) |
                (DNA_p_lower_null_model < 0.05 &
                 !!sym(paste0("DNA_p_lower_bb_model", dna_model)) >= 0.05),
          
          reg = ( !!sym(paste0("RNA_p_upper_bb_model", rna_model)) < 0.05 &
                  DNA_p_upper_null_model >= 0.05) |
                ( !!sym(paste0("RNA_p_lower_bb_model", rna_model)) < 0.05 &
                  DNA_p_lower_null_model >= 0.05)
        ) %>%
        
        # c) drop NAs, tally & pct
        tidyr::drop_na(CNA, reg) %>%
        dplyr::count(CNA, reg, name = "n") %>%
        dplyr::mutate(
          percent = n / sum(n) * 100,
          category = case_when(
            CNA &  reg  ~ "both",
            CNA & !reg  ~ "CNA-AI",
           !CNA &  reg  ~ "reg-AI",
            TRUE        ~ "neither"
          ),
          model_pair = model_pair,
          gene_type  = gene_type
        ) %>%
        dplyr::select(gene_type, model_pair, category, n, percent)
    })
})

results_all <- results_all %>%
  mutate(gene_type = recode(gene_type,
                            random = "Passenger genes"))

# Plot: stacked bars of % by model_pair, facetted over gene_type
plot <- ggplot(results_all, aes(x = model_pair, y = percent, fill = category)) +
  geom_col(position = "stack") +
  # add percent labels
  geom_text(aes(label = ifelse(percent >= 2, sprintf("%.0f%%", percent), "")),
            position = position_stack(vjust = 0.5),
            size = 5,         # bump text-in-bar size
            color = "white",
            fontface = "bold") +
  facet_wrap(~ gene_type, scales = "free_y") +
  labs(
    x    = "RNA–DNA bb‑model pair",
    y    = "Percent of significant mRNA‑AI",
    fill = "AI Source"
  ) +
  theme_classic(base_size = 16) +     # increase all base text
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.x       = element_text(size = 18, face = "bold"),
    axis.title.y       = element_text(size = 18, face = "bold"),
    strip.text         = element_text(size = 17, face = "bold"),  # facet titles
    legend.position    = "top",
    legend.title       = element_text(size = 16, face = "bold"),
    legend.text        = element_text(size = 15)
  ) +
  scale_fill_brewer(palette = "Set1", name = "AI source")

 # Save the plot
ggsave(
  "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/bb_models/prop_sig_DNA_vs_RNA_from_total_all_gene_types.png",
  plot,
  width  = 200,
  height = 200,
  units  = "mm",
  dpi    = 300
)

# B) % of significant reg vs CNA-AI from significant mRNA-AI SNVs for each BB-model --> Passenger genes & nonsense mut only

# Annotate NMD_status and filter to stopgains on random genes
df_nmd <- TCGA_ASE_table_bb_models %>%
  mutate(
    NMD_status = case_when(
      TSS_PTC_dist < 150 | X55_nt_last_exon == "NMD-evading" ~ "NMD-evading",
      TRUE                                             ~ "NMD-triggering"
    )
  ) %>%
  filter(
    SNV_varity == "stopgain",
    gene_cognate  == "random"
  )

# Your five RNA–DNA bb‑model pairings
model_pairs <- tibble::tribble(
  ~model_pair,   ~rna_model, ~dna_model,
  "RNA1-DNA1",   "1",        "1",
  "RNA2-DNA1",   "2",        "1",
  "RNA3-DNA2",   "3",        "2",
  "RNA4-DNA3",   "4",        "3",
  "RNA4-DNA4",   "4",        "4"
)

# Compute counts & percentages by NMD_status for each pairing
results_nmd <- model_pairs %>%
  pmap_dfr(function(model_pair, rna_model, dna_model) {
    df_nmd %>%
      # mRNA AI (null model p < 0.05)
      filter(RNA_p_upper_null_model < 0.05 |
             RNA_p_lower_null_model < 0.05) %>%
      # compute CNA & reg flags
      mutate(
        CNA = (DNA_p_upper_null_model < 0.05 &
               !!sym(paste0("DNA_p_upper_bb_model", dna_model)) >= 0.05) |
              (DNA_p_lower_null_model < 0.05 &
               !!sym(paste0("DNA_p_lower_bb_model", dna_model)) >= 0.05),
        reg = ( !!sym(paste0("RNA_p_upper_bb_model", rna_model)) < 0.05 &
                DNA_p_upper_null_model >= 0.05) |
              ( !!sym(paste0("RNA_p_lower_bb_model", rna_model)) < 0.05 &
                DNA_p_lower_null_model >= 0.05)
      ) %>%
      tidyr::drop_na(CNA, reg) %>%
      # count within each NMD_status
      dplyr::group_by(NMD_status) %>%
      dplyr::count(CNA, reg, name = "n") %>%
      dplyr::mutate(
        percent  = n / sum(n) * 100,
        category = case_when(
          CNA &  reg  ~ "both",
          CNA & !reg  ~ "CNA-AI",
         !CNA &  reg  ~ "reg-AI",
          TRUE        ~ "neither"
        ),
        model_pair = model_pair
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(model_pair, NMD_status, category, n, percent)
  })

# Plot
plot <- ggplot(results_nmd, aes(x = NMD_status, y = percent, fill = category)) +
  geom_col(position = "stack") +
  # only label bars ≥2%, others get empty label
  geom_text(aes(label = ifelse(percent >= 2, sprintf("%.0f%%", percent), "")),
            position = position_stack(vjust = 0.5),
            size = 5,         # bump text-in-bar size
            color = "white",
            fontface = "bold") +
  facet_wrap(~ model_pair) +
  labs(
    x    = "NMD status",
    y    = "% of nonsense SNVs with sig mRNA-AI",
    fill = "AI source"
  ) +
  theme_classic(base_size = 18) +     # increase all base text
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.x       = element_text(size = 18, face = "bold"),
    axis.title.y       = element_text(size = 18, face = "bold"),
    strip.text         = element_text(size = 17, face = "bold"),  # facet titles
    legend.position    = "top",
    legend.title       = element_text(size = 16, face = "bold"),
    legend.text        = element_text(size = 15)
  ) +
  # pick a “cool” palette: viridis or one from RColorBrewer
  # scale_fill_viridis_d(option = "C", name = "AI source")
  # OR, for RColorBrewer:
  scale_fill_brewer(palette = "Set1", name = "AI source")

# Save with metric dimensions
ggsave(
  "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/bb_models/prop_sig_DNA_vs_RNA_from_total_passengers_nonsense_mut.png",
  plot,
  width  = 200,
  height = 200,
  units  = "mm",
  dpi    = 300
)

# C) ORs of AI in NMD-triggering vs NMD-evading

# 1) Build the base NMD dataframe: annotate & filter stopgain/random
df_base_nmd <- TCGA_ASE_table_bb_models %>%
  mutate(
    NMD_status = case_when(
      TSS_PTC_dist < 150 | X55_nt_last_exon == "NMD-evading" ~ "NMD-evading",
      TRUE                                                  ~ "NMD-triggering"
    )
  ) %>%
  filter(
    SNV_varity == "stopgain",
    gene_cognate  == "random"
  )

# 2) Define your four bb‐model pairs + the new log_odds class
model_defs <- tibble::tribble(
  ~model_pair,   ~rna_model, ~dna_model,
  "RNA1-DNA1",   "1",        "1",
  "RNA2-DNA1",   "2",        "1",
  "RNA3-DNA2",   "3",        "2",
  "RNA4-DNA3",   "4",        "3",
  "log_odds",    NA,         NA
)

# 3) Compute ORs, CIs, p‐values & 2×2 counts for each “model_pair”
odds_results <- model_defs %>%
  pmap_dfr(function(model_pair, rna_model, dna_model) {
    df_tmp <- df_base_nmd
    
    if (model_pair != "log_odds") {
      # for the bb models, keep only mRNA-AI in the NULL model
      df_tmp <- df_tmp %>%
        filter(RNA_p_upper_null_model < 0.05 | RNA_p_lower_null_model < 0.05) %>%
        mutate(
          AI_sig = ( !!sym(paste0("RNA_p_upper_bb_model", rna_model)) < 0.05 &
                     DNA_p_upper_null_model >= 0.05 ) |
                   ( !!sym(paste0("RNA_p_lower_bb_model", rna_model)) < 0.05 &
                     DNA_p_lower_null_model >= 0.05 )
        )
    } else {
      # for the log_odds class, significance = abs(log_odds)>1
      df_tmp <- df_tmp %>% mutate(
        AI_sig = abs(log_odds) > 1
      )
    }
    
    df_tmp <- df_tmp %>% drop_na(AI_sig, NMD_status)
    
    # build the 2×2 table and run Fisher's test
    tab <- table(df_tmp$AI_sig, df_tmp$NMD_status)
    ft  <- fisher.test(tab)
    
    # return one row per model_pair
    tibble(
      model_pair   = model_pair,
      OR           = as.numeric(ft$estimate),
      lower_CI     = ft$conf.int[1],
      upper_CI     = ft$conf.int[2],
      p_value      = ft$p.value,
      n_evade_not = tab["FALSE", "NMD-evading"],
      n_evade_sig = tab["TRUE",  "NMD-evading"],
      n_trig_not  = tab["FALSE", "NMD-triggering"],
      n_trig_sig  = tab["TRUE",  "NMD-triggering"]
    )
  })

# 4) Pivot counts into long form for the grouped‐bar plot
counts_long <- odds_results %>%
  select(model_pair, starts_with("n_")) %>%
  pivot_longer(
    cols         = starts_with("n_"),
    names_to     = c("NMD_status","AI_sig"),
    names_pattern= "n_(evade|trig)_(not|sig)",
    values_to    = "count"
  ) %>%
  mutate(
    NMD_status = recode(NMD_status,
                        evade = "NMD-evading",
                        trig  = "NMD-triggering"),
    AI_sig     = recode(AI_sig,
                        not = "not-significant",
                        sig = "AI-significant")
  )

# 5) Build annotation df for OR + CI labels
annotation_df <- odds_results %>%
  mutate(
    label = paste0(
      "OR=", sprintf("%.2f", OR),
      "\n95% CI\n [", sprintf("%.2f", lower_CI), "-", sprintf("%.2f", upper_CI), "]"
    )
  ) %>%
  left_join(
    counts_long %>%
      group_by(model_pair) %>%
      summarize(y_max = max(count), .groups = "drop"),
    by = "model_pair"
  ) %>%
  distinct(model_pair, label, y_max)

# 6) Plot 2×2 counts with OR+CI annotations
plot <- ggplot(counts_long, aes(x = NMD_status, y = count, fill = AI_sig)) +
  # bars
  geom_col(position = position_dodge(0.8), width = 0.7) +
  
  # counts above bars, bumped up, larger text
  geom_text(aes(label = count),
            position = position_dodge(0.8),
            vjust    = 0,
            size     = 5) +
  
  # OR + CI annotation: italic, slightly lower-right, with newline
  geom_text(
    data        = annotation_df,
    aes(x        = 1, 
        y        = y_max * 0.8, 
        label    = label),
    inherit.aes = FALSE,
    size        = 5,
    fontface    = "italic",
    # hjust       = 0,              # anchor left of the point
    # vjust       = 1,              # anchor top of the text at the point
    # nudge_x     = 0.2,            # shift right
    # nudge_y     = -0.05 * annotation_df$y_max  # shift downward
  ) +
  
  # facets with bold titles
  facet_wrap(~ model_pair, ncol = 3, scales = "free_y") +
  
  # legend relabeled & moved
  scale_fill_manual(
    values = c("not-significant" = "grey70", "AI-significant" = "steelblue"),
    name   = "Significant AI",
    labels = c("Yes", "No")
  ) +
  
  labs(
    x     = "NMD status",
    y     = "Count of nonsense SNVs",
    title = "2×2 Contingency Counts of Significance vs NMD Status"
  ) +
  
  # bump all text sizes, bold y-axis title & facet titles
  theme_classic(base_size = 16) +
  theme(
    axis.text.x        = element_text(angle = 25, hjust = 1, size = 15),
    axis.title.x       = element_text(size = 18),
    axis.title.y       = element_text(size = 18, face = "bold"),
    strip.text         = element_text(size = 17, face = "bold"),
    legend.position    = "top",
    legend.title       = element_text(size = 16),
    legend.text        = element_text(size = 15),
    panel.grid.major.x = element_blank(),
    plot.title         = element_blank()
  )

ggsave(
  "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/bb_models/prop_sig_reg_AI_ORs_nonsense.png",
  plot,
  width  = 275,
  height = 225,
  units  = "mm",
  dpi    = 300
)

# 7) Plot ORs + CIs alone on a log scale
plot <- ggplot(odds_results, aes(x = model_pair, y = OR)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_log10() +
  labs(
    x     = "Model",
    y     = "Odds Ratio (log scale)",
    title = "ORs of Significance in NMD-Triggering vs NMD-Evading\n(including log_odds class)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/bb_models/prop_sig_reg_AI_ORs_nonsense_2.png",
  plot,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300
)
