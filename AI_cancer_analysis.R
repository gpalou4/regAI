################################################################################################
########################################## LIBRARIES ###########################################
################################################################################################

# conda activate R_figures
source("/home/gpalou/projects/ASE/Figures/ggplot_themes.R")
library(ggsignif)
library(GWASTools)
library(ggfortify)
library(VennDiagram)
library(survival)
library(survminer)
library(forestmodel)
library(ggvenn)
library(readxl)
library(tibble)

# setwd("F:/ASE_project/") 
setwd("/home/gpalou/projects/ASE")

################################################################################################
########################################## FUNCTIONS ###########################################
################################################################################################

ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g,
                                                               surv.plot.height = NULL,
                                                               risk.table.height = NULL,
                                                               ncensor.plot.height = NULL)}  

################################################################################################
########################################## SCRIPT ##############################################
################################################################################################

data_path <- "/g/strcombio/fsupek_cancer1/gpalou/ASE_project"

# 1) Open data
TCGA_AI_table_final <- readRDS("/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_AI_final_tables/TCGA_AI_table_bb_models_final.RData")
dim(TCGA_AI_table_final)
#854785    380

# 1.1) Significance thresholds

p_value <- 0.05
TCGA_AI_table_final <- TCGA_AI_table_final %>%
  dplyr::mutate(
    pos_mRNA_AI = ( RNA_p_upper_null_model < p_value ),
    neg_mRNA_AI = ( RNA_p_lower_null_model < p_value ),
    no_mRNA_AI = !(pos_mRNA_AI | neg_mRNA_AI),
    pos_CNA_AI = (DNA_p_upper_null_model < p_value &
            DNA_p_upper_bb_model1 >= p_value),
    neg_CNA_AI = (DNA_p_lower_null_model < p_value &
            DNA_p_lower_bb_model1 >= p_value),
    no_CNA_AI = !(pos_CNA_AI | neg_CNA_AI),
    pos_reg_AI = ( RNA_p_upper_bb_model1 < p_value &
            DNA_p_upper_null_model >= p_value),
    neg_reg_AI = ( RNA_p_lower_bb_model1 < p_value &
            DNA_p_lower_null_model >= p_value),
    no_reg_AI = !(pos_reg_AI | neg_reg_AI),
    all_mRNA_AI = (pos_mRNA_AI | neg_mRNA_AI),
    all_CNA_AI = (pos_CNA_AI | neg_CNA_AI),
    all_reg_AI = (pos_reg_AI | neg_reg_AI)
  ) %>%
  filter(!is.na(pos_reg_AI)) %>%
  # define the three reg-AI classes
  mutate(
    reg_AI_type = case_when(
      pos_reg_AI ~ "Pos reg-AI",
      neg_reg_AI ~ "Neg reg-AI",
      TRUE       ~ "No reg-AI"
    ),
    # make sure it's a factor in the order you want
    reg_AI_type = factor(
      reg_AI_type,
      levels = c( "Pos reg-AI", "No reg-AI", "Neg reg-AI"),
      labels = c("Pos reg-AI", "No reg-AI", "Neg reg-AI")
    )
  ) %>%
  # define the three CNA-AI classes
  mutate(
    CNA_AI_type = case_when(
      pos_CNA_AI ~ "Pos CNA-AI",
      neg_CNA_AI ~ "Neg CNA-AI",
      TRUE       ~ "No CNA-AI"
    ),
    # make sure it's a factor in the order you want
    CNA_AI_type = factor(
      CNA_AI_type,
      levels = c( "Pos CNA-AI", "No CNA-AI", "Neg CNA-AI"),
      labels = c("Pos CNA-AI", "No CNA-AI", "Neg CNA-AI")
    )
  ) %>%
  # define the three mRNA-AI classes
  mutate(
    mRNA_AI_type = case_when(
      pos_mRNA_AI ~ "Pos mRNA-AI",
      neg_mRNA_AI ~ "Neg mRNA-AI",
      TRUE       ~ "No mRNA-AI"
    ),
    # make sure it's a factor in the order you want
    mRNA_AI_type = factor(
      mRNA_AI_type,
      levels = c( "Pos mRNA-AI", "No mRNA-AI", "Neg mRNA-AI"),
      labels = c("Pos mRNA-AI", "No mRNA-AI", "Neg mRNA-AI")
    )
  ) %>%
  # switch to rowwise for c_across()
  rowwise() %>%
  mutate(
    Puffin_MaxDelta = {
      values <- c_across(c(
        PRO_CAP_50bp_TSS_sum_all_abs,
        FANTOM_CAGE_50bp_TSS_sum_all_abs,
        GRO_CAP_50bp_TSS_sum_all_abs
      ))
      idx <- which.max(abs(values))
      if (length(idx) == 0) 
        NA_real_ 
      else 
        values[idx]
    },
    Sei_MaxDelta = {
      values <- c_across(c(
        TN1.Transcription,
        TN2.Transcription,
        TN3.Transcription,
        TN4.Transcription
      ))
      max_abs_value <- max(abs(values), na.rm = TRUE)
      max_abs_value
    }

  ) %>%
  ungroup()

setDT(TCGA_AI_table_final)

TCGA_AI_table_final[, `:=`(
  # assign major/minor based on whether mutant fraction > 0.5
  MUT_CNA2 = ifelse(
    MUT_CNA_fraction > 0.5,
    # estimate_MUT_CNA >= 2,
    Major_Copy_Number,
    Minor_Copy_Number
  ),
  WT_CNA2  = ifelse(
    MUT_CNA_fraction > 0.5,
    # estimate_MUT_CNA >= 2,
    Minor_Copy_Number,
    Major_Copy_Number
  )
)]

TCGA_AI_table_final <- data.frame(TCGA_AI_table_final)

dim(TCGA_AI_table_final)
table(TCGA_AI_table_final$mRNA_AI_type,TCGA_AI_table_final$CNA_AI_type)
table(TCGA_AI_table_final$mRNA_AI_type,TCGA_AI_table_final$reg_AI_type)
table(TCGA_AI_table_final$mRNA_AI_type)
summary(TCGA_AI_table_final$Sei_MaxDelta)
summary(TCGA_AI_table_final$Puffin_MaxDelta)
summary(TCGA_AI_table_final$promoterAI)

# 1.2) Reclassify cancer gene sets based on CGC+Mutpanning

# CGC genes
CGC_genes <- read.csv("/g/strcombio/fsupek_home/ppericot/ASE_project/cancer_gene_census_updated_ensembl.tsv", header = T, sep = "\t") 
# Mutpanning genes
Mutpanning_genes <- read.csv("/g/strcombio/fsupek_home/ppericot/ASE_project/MutPanningGeneTumorPairs.csv")
colnames(Mutpanning_genes)[1] <- "Gene.Symbol"
table(unique(CGC_genes$Gene.Symbol) %in% unique(Mutpanning_genes$Gene.Symbol))
table(unique(Mutpanning_genes$Gene.Symbol) %in% unique(CGC_genes$Gene.Symbol))
# Sets
CGC_Mutpanning_shared <- unique(CGC_genes$Gene.Symbol)[unique(CGC_genes$Gene.Symbol) %in% unique(Mutpanning_genes$Gene.Symbol)]
CGC_only <- unique(CGC_genes$Gene.Symbol)[!unique(CGC_genes$Gene.Symbol) %in% unique(Mutpanning_genes$Gene.Symbol)]
Mutpanning_only <- unique(Mutpanning_genes$Gene.Symbol)[!unique(Mutpanning_genes$Gene.Symbol) %in% unique(CGC_genes$Gene.Symbol)]
CGC_trans_fus_genes <- CGC_genes[!str_detect(CGC_genes$Mutation.Types, patter ="Mis|N|S"),]
table(CGC_trans_fus_genes$Role.in.Cancer)
table(CGC_trans_fus_genes$Mutation.Types)

dim(TCGA_AI_table_final)
table(TCGA_AI_table_final$gene_cognate)

#For the Manuscript
TCGA_AI_table_final %>%
  filter(gene == "oncogene") %>%
  pull(Gene.Symbol) %>% unique() %>% length()

TCGA_AI_table_final <- TCGA_AI_table_final %>%
  # Remove 'random' or 'essential' gene_cognate entries that are in Mutpanning_only
  filter(!(gene_cognate %in% c("random", "essential") & Gene.Symbol %in% Mutpanning_only)) %>%
  # Remove CGC trans/fusion genes based on gene_cognate and Gene.Symbol or Synonyms
  filter(!(gene_cognate %in% c("og_cognate", "oncogene", "tsg", "tsg_cognate") &
           (Gene.Symbol %in% CGC_trans_fus_genes$Gene.Symbol |
            gene_id %in% CGC_trans_fus_genes$Synonyms)))
dim(TCGA_AI_table_final)
table(TCGA_AI_table_final$gene_cognate)

# removed_genes <- unique(TCGA_AI_table_final$Gene.Symbol)[!unique(TCGA_AI_table_final$Gene.Symbol) %in% unique(df$Gene.Symbol)]

# For the Manuscript - counting genes
TCGA_AI_table_final %>%
  group_by(gene_cognate) %>%
  summarise(
    n_unique_symbols = n_distinct(Gene.Symbol)
  ) %>%
  arrange(desc(n_unique_symbols))

df <- unique(TCGA_AI_table_final[,c("chrom","pos","Gene.Symbol")])
dim(df)

# Save
path <- "/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_AI_final_tables/TCGA_AI_table_bb_models_final_2.RData"
# saveRDS(TCGA_AI_table_final, path)
TCGA_AI_table_final <- readRDS(path)
dim(TCGA_AI_table_final)

# df <- TCGA_AI_table_final %>%
#   filter(gene == "random") %>%
#   filter(SNV_varity == "stopgain")
# dim(df)
# table(df$NMD_status)

# prop.table(table(TCGA_AI_table_final$all_reg_AI))*100
# prop.table(table(abs(TCGA_AI_table_final$log_odds) > 1))*100

###################################################################################################
################################# VARIABILITY REG-AI VS CNA-AI ####################################
###################################################################################################

###################################################
################### FIGURE 1A #####################
###################################################

# mRNA-AI vs reg-AI vs CNA-AI ppt schematic

###################################################
################### FIGURE 1B #####################
###################################################

# --- Step 1: Prepare data for plotting ---

df_processed <- TCGA_AI_table_final %>%
  mutate(cancer = gsub("TCGA-", "", cancer))

# Filter for significant mRNA-AI SNVs and categorize them
plot_df <- df_processed %>%
  filter(all_mRNA_AI) %>%
  mutate(
    AI_source = case_when(
      # all_reg_AI & all_CNA_AI   ~ "Both reg-AI & CNA-AI",
      all_reg_AI & !all_CNA_AI  ~ "reg-AI",
      !all_reg_AI & all_CNA_AI  ~ "CNA-AI"
      # TRUE                    ~ "mRNA-AI only" # Significant mRNA AI, but not meeting reg or CNA criteria
    )
  ) %>%
  filter(!is.na(AI_source)) %>%
  group_by(cancer, AI_source) %>%
  summarise(n = n(), .groups = "drop") %>% # Count SNVs in each category per cancer
  group_by(cancer) %>%
  mutate(
    total_snvs_in_cancer = sum(n),
    percent = (n / total_snvs_in_cancer) * 100
  ) %>%
  ungroup()

category_levels <- c("reg-AI", "CNA-AI")
plot_df$AI_source <- factor(plot_df$AI_source, levels = category_levels)

# Order cancer types by the total percentage of reg-AI ("reg-AI only" + "Both")
cancer_order <- plot_df %>%
  filter(AI_source %in% c("reg-AI")) %>%
  group_by(cancer) %>%
  summarise(total_reg_AI_percent = sum(percent), .groups = "drop") %>%
  arrange(total_reg_AI_percent) %>% # Or arrange(desc(total_reg_AI_percent)) for descending
  pull(cancer)

plot_df$cancer <- factor(plot_df$cancer, levels = cancer_order)

# --- Step 2: Plot the data ---
category_colors <- c(
  "reg-AI"           = "#A191B6", # Greenish
  "CNA-AI"           = "#539753" # Orangish
)

# Create the plot
plot <- ggplot(plot_df, aes(x = cancer, y = percent, fill = AI_source)) +
  geom_col(position = "stack", width = 0.8) + # geom_col is equivalent to geom_bar(stat="identity")
  geom_text(
    aes(label = if_else(percent > 5, sprintf("%.0f%%", percent), "")), # Show label if > 5%
    position = position_stack(vjust = 0.5),
    size = 3, # Adjust size as needed
    color = "white", # Color for text, ensure contrast
    fontface = "bold"
  ) +
  coord_flip() + # Flip coordinates to make bars horizontal
  scale_fill_manual(
    values = category_colors,
    name = "Source of mRNA-AI", # Legend title
    labels = c( # Ensure labels match the levels if you want to customize
        "reg-AI only"           = "Regulatory AI only",
        "CNA-AI only"           = "CNA-driven AI only"
    ),
    limits = category_levels # Ensures legend order matches factor levels
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) + # Format y-axis as percentage
  labs(
    title = "Proportion of AI Sources for Significant mRNA-AI SNVs",
    x = "Cancer Type",
    y = "Percent of SNVs"
  ) +
  theme_grey(base_size = 14) + # Starting with a grey theme similar to original
  theme(
    axis.title.y = element_text(size = 15, margin = margin(r = 10)), # Cancer Type label
    axis.title.x = element_text(size = 15, margin = margin(t = 10)), # Percent of SNVs label
    axis.text.x = element_text(size = 10), 
    axis.text.y = element_text(size = 10, face = "bold"), # Cancer names
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"), # Title
    legend.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    strip.text = element_text(face = "bold") # For facets, if any were used
  )

## For the Manuscript (and plot) ##

plot_df %>%
  filter(AI_source == "reg-AI") %>%
  arrange((percent))
# plot_df %>%
#   filter(AI_source == "mRNA-AI") %>%
#   arrange((percent))
plot_df %>%
  filter(AI_source == "CNA-AI") %>%
  arrange(desc(percent))

# Calculate pan-cancer proportions
pancan_summary_df <- df_processed %>%
  filter(all_mRNA_AI) %>%
  mutate(
    AI_source = case_when( 
      all_reg_AI & !all_CNA_AI  ~ "reg-AI",
      !all_reg_AI & all_CNA_AI  ~ "CNA-AI"
    )
  ) %>%
  filter(!is.na(AI_source)) %>%
  group_by(AI_source) %>%
  summarise(
    n = n(),
    .groups = "drop" 
  ) %>%
  mutate(
    total_exclusive_snvs_pancan = sum(n),
    percent = (n / total_exclusive_snvs_pancan) * 100
  )
pancan_summary_df

# Alternative (current paper)
df <- TCGA_AI_table_final %>%
  filter(all_mRNA_AI & (all_reg_AI | all_CNA_AI))
dim(df)
# 41549 / (747634)*100
table(df$all_mRNA_AI)
prop.table(table(df$all_reg_AI,df$all_mRNA_AI))
prop.table(table(df$all_CNA_AI,df$all_mRNA_AI))

# Save
write.table(plot_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig1/Fig1B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(plot_df, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig1/Fig1B.RData")

######################################################
################### FIGURE 1 C-D #####################
######################################################

# Define colors and shapes
ai_type_colors <- c("reg-AI" = "#A191B6", "CNA-AI" = "#539753", "no-AI" = "grey80")
mrna_ai_shapes <- c("TRUE" = 17, "FALSE" = 16) # Triangle for TRUE, Circle for FALSE

# --- Plot 1: OV ---

# Prepare data for Plot 1
plot1_data <- TCGA_AI_table_final %>%
  filter(cancer == "TCGA-OV") %>%
  # Filter out NAs from columns used for aesthetics
  filter(!is.na(all_reg_AI), !is.na(all_CNA_AI), !is.na(all_mRNA_AI)) %>%
  mutate(
    # Create the 3-category AI type for color, prioritizing reg-AI
    ai_type_color = factor(case_when(
      all_reg_AI == TRUE                       ~ "reg-AI",
      all_reg_AI == FALSE & all_CNA_AI == TRUE ~ "CNA-AI",
      all_reg_AI == FALSE & all_CNA_AI == FALSE~ "no-AI",
      TRUE                                     ~ "Other"
    ), levels = c("reg-AI", "CNA-AI", "no-AI", "Other")),
    
    # Factor for shape aesthetic
    all_mRNA_AI_shape = factor(all_mRNA_AI, levels = c(TRUE, FALSE))
  ) %>%
  filter(ai_type_color != "Other")

# Create Plot 1
plot1 <- ggplot(plot1_data, aes(x = DNA_observed_ratio, y = RNA_observed_ratio)) +
  geom_point(aes(color = ai_type_color, shape = all_mRNA_AI_shape), alpha = 0.7, size = 2.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = ai_type_colors,
    name = "significant AI",
    labels = c("reg-AI" = "reg-AI", "CNA-AI" = "CNA-AI", "no-AI" = "neither")
  ) +
  scale_shape_manual(
    values = mrna_ai_shapes,
    name = "mRNA-AI",
    labels = c("TRUE" = "significant", "FALSE" = "ns")
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "OV",
    x = "DNA VAF",
    y = "RNA VAF"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5)
  )

# --- Plot 2: UCEC-MSS ---

# MSI status
sample_NMD_efficiencies_TCGA_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt"
sample_NMD_efficiencies_TCGA <- read.table(file = sample_NMD_efficiencies_TCGA_path, header = TRUE, sep = "\t")
MSI_POLE_tmp <- sample_NMD_efficiencies_TCGA[,c("sample","MSI_status","POLE")]
TCGA_AI_table_final_filt <- merge(TCGA_AI_table_final,MSI_POLE_tmp)

# Prepare data for Plot 2
plot2_data <- TCGA_AI_table_final_filt %>%
  filter(cancer == "TCGA-UCEC") %>%
  filter(POLE != "SBS10ab") %>%
  filter(MSI_status != "MSI-H") %>%
  # Filter out NAs from columns used for aesthetics
  filter(!is.na(all_reg_AI), !is.na(all_CNA_AI), !is.na(all_mRNA_AI)) %>%
  mutate(
    # Create the 3-category AI type for color, prioritizing reg-AI
    ai_type_color = factor(case_when(
      all_reg_AI == TRUE                       ~ "reg-AI",
      all_reg_AI == FALSE & all_CNA_AI == TRUE ~ "CNA-AI",
      all_reg_AI == FALSE & all_CNA_AI == FALSE~ "no-AI",
      TRUE                                     ~ "Other" # Should not be reached if NAs are filtered
    ), levels = c("reg-AI", "CNA-AI", "no-AI", "Other")),
    
    # Factor for shape aesthetic
    all_mRNA_AI_shape = factor(all_mRNA_AI, levels = c(TRUE, FALSE))
  ) %>%
  filter(ai_type_color != "Other") # Remove any "Other" cases if they somehow occur

# Create Plot 1
plot2 <- ggplot(plot2_data, aes(x = DNA_observed_ratio, y = RNA_observed_ratio)) +
  geom_point(aes(color = ai_type_color, shape = all_mRNA_AI_shape), alpha = 0.7, size = 2.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = ai_type_colors,
    name = "significant AI",
    labels = c("reg-AI" = "reg-AI", "CNA-AI" = "CNA-AI", "no-AI" = "neither")
  ) +
  scale_shape_manual(
    values = mrna_ai_shapes,
    name = "mRNA-AI",
    labels = c("TRUE" = "significant", "FALSE" = "ns")
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title = "UCEC",
    x = "DNA VAF", # Variant Allele Fraction
    y = "RNA VAF"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical", # Arrange legends vertically if they overlap
    plot.title = element_text(hjust = 0.5)
  )

# --- Combine plots side-by-side using patchwork ---
plot <- plot1 + plot2

# Save
write.table(plot1_data, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig1/Fig1C.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(plot1_data, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig1/Fig1C.RData")

write.table(plot2_data, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig1/Fig1D.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(plot2_data, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig1/Fig1D.RData")

####################################################################################################
################################# MECHANISMS OF REG-AI ANALYSIS ####################################
####################################################################################################

# THIS INCLUDES: SEI/PUFFIN-D/SPLICING

# For the manuscript

# 1) Proportions of spliceAI and Pangolin vs reg-AI

prop.table(table(TCGA_AI_table_final$DS_max_abs > 0.15))
prop.table(table(TCGA_AI_table_final$Pangolin_abs_max > 0.10))
prop.table(table(TCGA_AI_table_final$DS_max_abs > 0.15,TCGA_AI_table_final$reg_AI_type))

subset_data <- subset(TCGA_AI_table_final, reg_AI_type == "Pos reg-AI")
prop.table(table(subset_data$Pangolin_abs_max > 0.1))
prop.table(table(subset_data$DS_max_abs > 0.15))
subset_data <- subset(TCGA_AI_table_final, reg_AI_type == "Neg reg-AI")
prop.table(table(subset_data$Pangolin_abs_max > 0.1))
prop.table(table(subset_data$DS_max_abs > 0.15))
subset_data <- subset(TCGA_AI_table_final, reg_AI_type == "No reg-AI")
prop.table(table(subset_data$Pangolin_abs_max > 0.1))
prop.table(table(subset_data$DS_max_abs > 0.15))

# 2) TP53 examples
# scale of deltaScore for 3 well known examples
# From Fran's paper, recurrently syn mutations in TP53. In contrast to those in oncogenes, these mutations inactivate splice sites.
TCGA_AI_table_final_TP53 <- TCGA_AI_table_final %>%
  filter(gene=="TP53") %>%
  filter(SNV_varity == "effectively_syn")
# A) Codon 125 exon 4 --> Intron retention or cryptic splice site
df <- TCGA_AI_table_final_TP53[grep("125",TCGA_AI_table_final_TP53$AAChange),]
df <- df[grep("exon4",df$AAChange),]
df[1,]
dim(df)
table(as.character(df$alt))
df$SNV_varity_new
mean(df$DS_DG)
mean(df$DS_AG)
mean(df$DS_AL)
mean(df$DS_DL)
mean(df$larg_inc)
mean(df$larg_dec)
# table(df$log_odds < -1)
# table(df$log_odds_cna < -1)
table(df$reg_AI_type)
table(df$all_CNA_AI)

# B) 3' of exon 6 --> cryptic splice site and frameshift in mRNA
df <- TCGA_AI_table_final_TP53[grep("224",TCGA_AI_table_final_TP53$AAChange),]
df <- df[grep("exon6",df$AAChange),]
df[1,]
dim(df)
df$SNV_varity_new
mean(df$DS_DG)
mean(df$DS_AG)
mean(df$DS_AL)
mean(df$DS_DL)
mean(df$larg_inc)
mean(df$larg_dec)
# table(df$log_odds < -1)
# table(df$log_odds_cna < -1)
table(df$reg_AI_type)
table(df$all_CNA_AI)

# 3′ terminal G of exon 9 is commonly lost. 
df <- TCGA_AI_table_final_TP53[grep("331",TCGA_AI_table_final_TP53$AAChange),]
df <- df[grep("exon9",df$AAChange),]
df[1,]
dim(df)
df$SNV_type
mean(df$DS_DG)
mean(df$DS_AG)
mean(df$DS_AL)
mean(df$DS_DL)
mean(df$larg_inc)
mean(df$larg_dec)
# table(df$log_odds < -1)
# table(df$log_odds_cna < -1)
table(df$reg_AI_type)
table(df$all_CNA_AI)

TCGA_AI_table_final_TP53[1,]

# 3) Proportions of SpliceAI/Pangolin and reg-AI by SNV type

df <- TCGA_AI_table_final #%>%
  # filter(gene_cognate == "random")  #%>%
  # filter(SNV_varity != "stopgain")
  # filter(SNV_varity_new == "Synonymous")

subset_data <- subset(df, reg_AI_type == "Neg reg-AI")
dim(subset_data)
prop.table(table(subset_data$Pangolin_abs_max >= 0.1))
prop.table(table(subset_data$DS_max_abs >= 0.15))

########################################################
################### FIGURE 2A ##########################
########################################################

# Sei/Puffin/Splicing reg-AI schematic --> Done

########################################################
################### FIGURE 2B ##########################
########################################################

# Add Positive control mutations (TERT/CDC20/FOXA1)
pos_control <- read.table(paste0("/g/strcombio/fsupek_fisher/gpalou/Puffin/batches/output/artificial_sequences/pos_control_puffin_output.txt"),
                    header = TRUE, sep = "\t")
pos_control <- pos_control %>%
  # switch to rowwise for c_across()
  # filter(Gene.Symbol == "TERT") %>%
  rowwise() %>%
  mutate(
    Puffin_MaxDelta = {
      values <- c_across(c(
        PRO_CAP_50bp_TSS_sum_all_abs,
        FANTOM_CAGE_50bp_TSS_sum_all_abs,
        GRO_CAP_50bp_TSS_sum_all_abs
      ))
      idx <- which.max(abs(values))
      if (length(idx) == 0) 
        NA_real_ 
      else 
        values[idx]
    }
  ) %>%
  ungroup()

options(scipen = 999)

# A) Puffin-D 
# Q1 threshold
puffin_thres <- quantile(pos_control$Puffin_MaxDelta,seq(0,1,0.05))["25%"]
puffin_thres

# B) Sei

get_top_pct <- function(x, pct) {
  stopifnot(is.numeric(x), pct >= 0, pct <= 100)
  # 1 - pct/100 is the quantile for the bottom (1 - pct) fraction
  q <- quantile(x, probs = 1 - pct/100, na.rm = TRUE)
  unname(q)
}
# # A) Top 5%
top_thres <- 5
sei_thres <- get_top_pct(TCGA_AI_table_final$Sei_MaxDelta, top_thres)
sei_thres

df <- TCGA_AI_table_final %>%
  mutate(
    Puffin_group = case_when(
      Puffin_MaxDelta >= puffin_thres                          ~ "TN-impact",
      Puffin_MaxDelta < puffin_thres                        ~ "TN-neutral",
      TRUE                                           ~ NA_character_
    ),
    Puffin_group = factor(
      Puffin_group,
      levels = c(
        "TN-neutral",
        "TN-impact"
      )
    ),
    Sei_group = case_when(
      Sei_MaxDelta >= sei_thres                           ~ "TN-impact",
      Sei_MaxDelta < sei_thres                        ~ "TN-neutral",
      TRUE                                           ~ NA_character_
    ),
    Sei_group = factor(
      Sei_group,
      levels = c(
        "TN-neutral",
        "TN-impact"
      )
    ),
    pangolin_group = ifelse(Pangolin_abs_max >= 0.1, "High-impact",
                            ifelse(Pangolin_abs_max < 0.1, "Low-impact", NA_character_)),
    spliceAI_group = ifelse(DS_max_abs >= 0.15, "High-impact",
                            ifelse(DS_max_abs < 0.15, "Low-impact", NA_character_))
  )
table(df$Sei_group)
table(df$Puffin_group)
table(df$spliceAI_group)
table(df$pangolin_group)

plot_df <- df %>%
    filter(!is.na(SNV_varity)) %>%
    filter(SNV_varity != "stopgain") %>%   
    filter(!gene_cognate %in% c("both","both_cognate","essential","tsg","oncogene")) %>%
    mutate(gene_cognate = factor(as.character(gene_cognate)))  
table(plot_df$gene_cognate)

# — 1) Re‐build your raw data frame including “All genes” —
plot_full <- bind_rows(
  plot_df,
  plot_df %>% mutate(gene_cognate = "All genes")
) %>%
  mutate(
    gene_cognate = factor(
      gene_cognate,
      levels = c("All genes", levels(plot_df$gene_cognate))
    )
  ) #%>%
  # drop any missing group calls
  # filter(!is.na(Sei_group),
  #        !is.na(Puffin_group),
  #        !is.na(spliceAI_group))

# 2) Summarise each variable’s pct off its own non-NA denominator
summary_df <- plot_full %>%
  group_by(gene_cognate, reg_AI_type) %>%
  summarise(
    ## SpliceAI
    spliceAI_num = sum(!is.na(DS_max_abs)       & DS_max_abs       >= 0.15),
    spliceAI_den = sum(!is.na(DS_max_abs)),
    pct_spliceAI = spliceAI_num / spliceAI_den * 100,

    ## Pangolin
    pangolin_num = sum(!is.na(Pangolin_abs_max) & Pangolin_abs_max >= 0.1),
    pangolin_den = sum(!is.na(Pangolin_abs_max)),
    pct_pangolin = pangolin_num / pangolin_den * 100,

    ## Puffin
    puffin_num   = sum(!is.na(Puffin_MaxDelta)  & Puffin_MaxDelta  >= puffin_thres),
    puffin_den   = sum(!is.na(Puffin_MaxDelta)),
    pct_puffin   = puffin_num   / puffin_den   * 100,

    ## SEI
    sei_num      = sum(!is.na(Sei_MaxDelta)     & Sei_MaxDelta     >= sei_thres),
    sei_den      = sum(!is.na(Sei_MaxDelta)),
    pct_sei      = sei_num      / sei_den      * 100,

    ## “Other” just as remainder
    pct_other    = 100 - (pct_spliceAI + pct_pangolin + pct_puffin + pct_sei),
    .groups = "drop"
  ) %>%
  # pivot only the pct_ columns for plotting; the *_num and *_den stay wide
  pivot_longer(
    cols     = starts_with("pct_"),
    names_to = "category",
    values_to= "pct"
  ) %>%
  mutate(
    category = factor(
      category,
      levels = c("pct_sei",
                 "pct_puffin",
                 "pct_spliceAI",
                 "pct_pangolin",
                 "pct_other"),
      labels = c("High-Impact SEI",
                 "High-Impact Puffin",
                 "High-Impact SpliceAI",
                 "High-Impact Pangolin",
                 "Other")
    )
  )

# 1) Build a long data.frame with three “big” facets:
#    - High-Impact SEI        (one slice: Sei)
#    - High-Impact Puffin     (one slice: Puffin-D)
#    - High-Impact Splicing   (three slices: SP only, PG only, overlap)
df_plot <- bind_rows(
  # SEI facet
  summary_df %>%
    transmute(
      gene_cognate,
      reg_AI_type,
      facet_group = "High-Impact SEI",
      subcat      = "Sei",
      pct         = pct_sei
    ),

  # Puffin facet
  summary_df %>%
    transmute(
      gene_cognate,
      reg_AI_type,
      facet_group = "High-Impact Puffin",
      subcat      = "Puffin-D",
      pct         = pct_puffin
    ),

  # Splicing facet (3 slices)
  summary_df %>%
    transmute(
      gene_cognate,
      reg_AI_type,
      facet_group = "High-Impact Splicing",
      subcat      = "SpliceAI only",
      pct         = pct_sp_only
    ),
  summary_df %>%
    transmute(
      gene_cognate,
      reg_AI_type,
      facet_group = "High-Impact Splicing",
      subcat      = "Pangolin only",
      pct         = pct_pg_only
    ),
  summary_df %>%
    transmute(
      gene_cognate,
      reg_AI_type,
      facet_group = "High-Impact Splicing",
      subcat      = "SpliceAI & Pangolin",
      pct         = pct_sp_pg
    )
) %>%
  mutate(
    facet_group = factor(
      facet_group,
      levels = c("High-Impact SEI",
                 "High-Impact Puffin",
                 "High-Impact Splicing")
    ),
    subcat = factor(
      subcat,
      levels = c("Sei",
                 "Puffin-D",
                 "SpliceAI only",
                 "Pangolin only",
                 "SpliceAI & Pangolin")
    )
  )

# 2) Define your palette:
base_pal     <- brewer.pal(5, "Set2")
pal_main     <- base_pal[-4]   # drops the pink
names(pal_main) <- c("Sei", "Puffin-D", "SpliceAI only", "Pangolin only")
pal_overlap  <- base_pal[4]    # the pink
pal_full     <- c(pal_main,
                  "SpliceAI & Pangolin" = pal_overlap)

# 3) Plot
plot <- ggplot(df_plot, aes(x = reg_AI_type, y = pct, fill = subcat)) +
  geom_col() +
  facet_grid(facet_group ~ gene_cognate,
             scales = "free_y", switch = "y") +
  scale_fill_manual(values = pal_full, name = NULL) +
  coord_cartesian(ylim = c(0, 30)) +
  labs(x = NULL, y = "% of SNVs") +
  theme_classic() +
  theme(
    # no x axis text or ticks
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),

    # facet strip styling
    strip.background = element_blank(),
    strip.placement  = "outside",
    strip.text.y     = element_text(face = "bold"),

    # unified legend at bottom
    legend.position  = "bottom"
  )

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, plot, width = 150, height = 100, units = "mm") 

# Save
# write.table(plot_list, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig2/Fig2B.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(df_plot, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig2/Fig2B.RData")

# For the Manuscript (proportions %)
summary_df %>%
  filter(gene_cognate == "tsg_cognate") %>%
  filter(category == "High-Impact SEI") %>% data.frame()

df <- TCGA_AI_table_final %>%
  filter(gene_cognate == "og_cognate") %>%
  filter(SNV_varity != "stopgain")

# subset_data <- subset(df, reg_AI_type == "Pos reg-AI")
# round(prop.table(table(subset_data$DS_max_abs > 0.15)),2)
# prop.table(table(subset_data$Puffin_MaxDelta >= puffin_thres))
# subset_data <- subset(df, reg_AI_type == "Neg reg-AI")
# round(prop.table(table(subset_data$DS_max_abs > 0.15)),2)
# prop.table(table(subset_data$Puffin_MaxDelta >= puffin_thres))
# subset_data <- subset(df, reg_AI_type == "No reg-AI")
# round(prop.table(table(subset_data$DS_max_abs > 0.15)),2)

enrich_all <- function(df, category_label, gene = "tsg_cognate") {
  # map category labels to your num/den prefixes
  prefix_map <- c(
    "High-Impact SEI"       = "sei",
    "High-Impact Puffin"    = "puffin",
    "High-Impact SpliceAI"  = "spliceAI",
    "High-Impact Pangolin"  = "pangolin"
  )
  pf  <- prefix_map[category_label]
  num <- paste0(pf, "_num")
  den <- paste0(pf, "_den")
  
  # filter down to the rows you care about
  sub <- df %>%
    filter(category == category_label,
           gene_cognate == gene)
  
  # helper to pull (hi, total) for one reg_AI_type
  get_counts <- function(rt) {
    d <- sub %>% filter(reg_AI_type == rt)
    hi    <- sum(d[[num]], na.rm = TRUE)
    total <- sum(d[[den]], na.rm = TRUE)
    c(hi = hi, total = total)
  }
  
  pos  <- get_counts("Pos reg-AI")
  neg  <- get_counts("Neg reg-AI")
  none <- get_counts("No reg-AI")
  all  <- pos + neg
  
  run_tests <- function(counts1, counts2, name1, name2) {
    x <- c(counts1["hi"], counts2["hi"])
    n <- c(counts1["total"], counts2["total"])
    # chi-square test
    chi <- prop.test(x, n)
    # Fisher’s exact
    tbl <- matrix(
      c(x[1], n[1] - x[1],
        x[2], n[2] - x[2]),
      nrow = 2, byrow = TRUE,
      dimnames = list(
        Group   = c(name1, name2),
        Outcome = c("High", "NotHigh")
      )
    )
    fisher <- fisher.test(tbl)
    list(chi_sq = chi, fisher = fisher)
  }
  
  list(
    all_vs_none = run_tests(all,  none, "Reg-AI",     "No reg-AI"),
    pos_vs_none = run_tests(pos,  none, "Pos reg-AI", "No reg-AI"),
    neg_vs_none = run_tests(neg,  none, "Neg reg-AI", "No reg-AI")
  )
}

# Usage:
res <- enrich_all(summary_df, "High-Impact SEI", "tsg_cognate")

# Access results:
res$all_vs_none$chi_sq    # χ² on (Pos+Neg) vs None
res$all_vs_none$fisher    # Fisher on (Pos+Neg) vs None

res$pos_vs_none$chi_sq    # χ² on Pos vs None
res$pos_vs_none$fisher    # Fisher on Pos vs None

res$neg_vs_none$chi_sq    # χ² on Neg vs None
res$neg_vs_none$fisher    # Fisher on Neg vs None

########################################################
################### SUPP FIGURE 1A #####################
########################################################

splice_type <- "spliceAI" # spliceAI

df <- TCGA_AI_table_final %>%
  # Splicing groups
  mutate(
      pangolin_group = ifelse(Pangolin_abs_max >= 0.1, "High-impact",
                              ifelse(Pangolin_abs_max < 0.1, "Low-impact", NA_character_)),
      spliceAI_group = ifelse(DS_max_abs >= 0.15, "High-impact",
                              ifelse(DS_max_abs < 0.15, "Low-impact", NA_character_))
      )

plot_df <- df %>%
    filter(!is.na(SNV_varity)) %>%
    # filter(SNV_varity != "stopgain") %>%   
    filter(!gene_cognate %in% c("both","both_cognate","essential","tsg","oncogene")) %>%
    mutate(gene_cognate = factor(as.character(gene_cognate)))  
table(plot_df$gene_cognate)

plot_df_all <- bind_rows(
    plot_df,
    plot_df %>% mutate(gene_cognate = "All genes")
  ) %>%
  # if you want a sensible factor ordering
  mutate(
    gene_cognate = factor(
      gene_cognate,
      levels = c("All genes", levels(plot_df$gene_cognate))
    )
  )

# check
table(plot_df_all$gene_cognate)

if (splice_type == "pangolin") {
  plot_df_all <- plot_df_all %>%
    filter(!is.na(pangolin_group)) %>%
    # count & pct within each gene_cognate × reg_AI_type
    dplyr::count(gene_cognate, reg_AI_type, pangolin_group, name = "n") %>%
    # group_by(gene_cognate, pangolin_group) %>%
    group_by(gene_cognate, reg_AI_type) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()
  } else if ( splice_type == "spliceAI") {
    plot_df_all <- plot_df_all %>%
      filter(!is.na(spliceAI_group)) %>%
      # count & pct within each gene_cognate × reg_AI_type
      dplyr::count(gene_cognate, reg_AI_type, spliceAI_group, name = "n") %>%
      # group_by(gene_cognate, spliceAI_group) %>%
      group_by(gene_cognate, reg_AI_type) %>%
      mutate(pct = n / sum(n) * 100) %>%
      ungroup()
  }
colnames(plot_df_all)[3] <- "splicing_group"

plot <- ggplot(plot_df_all,
            aes(x = reg_AI_type,
                y = pct,
                fill = splicing_group)) +
                # fill = reg_AI_type)) +
  geom_col() +
  geom_text(aes(label = paste0(round(pct, 0), "%")),
            position = position_stack(vjust = 0.5),
            size = 3) +
  # two‐row facet: pos_reg_AI on top, neg_reg_AI on bottom
  facet_grid(. ~ gene_cognate, switch = "y") +
  # zoom to 0–100% without dropping any bars
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(
    x    = "",
    y    = "% of SNVs",
    fill = "Splicing impact"
  ) +
  theme_classic() +
  # scale_fill_viridis(discrete = T) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.placement  = "outside",
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold")
  )

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, plot, width = 150, height = 100, units = "mm") 

# Save
write.table(plot_df_all, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig2/Fig2B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(plot_df_all, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig2/Fig2B.RData")

##########################################################
################### SUPP. FIGURE S1A #####################
##########################################################

# TP53 plot (Alternatives using ggbio or gviz)

# 1) TP53 exon coordinates

gtf_path <- "/g/strcombio/fsupek_cancer1/gpalou/ASE_project/conversor_tables/gencode.v26.annotation.gtf"

gtf_df <- read_tsv(
  file       = gtf_path,
  comment    = "#",                     # drop any line beginning with "#"
  col_names  = c(                       # supply the 9 GTF standard columns
    "seqname", "source", "feature",
    "start",   "end",    "score",
    "strand",  "frame",  "attribute"
  ),
  col_types  = cols(
    seqname  = col_character(),
    source   = col_character(),
    feature  = col_character(),
    start    = col_integer(),
    end      = col_integer(),
    score    = col_character(),         # often “.” in GTF
    strand   = col_character(),
    frame    = col_character(),
    attribute= col_character()
  )
)

tp53_main <- gtf_df %>%
  # 1) only exons for the exact transcript
  filter(
    feature == "exon",
    str_detect(attribute, 'ENST00000269305')
  ) %>%
  # 2) derive everything in one go
  transmute(
    # strip "chr" if present, then coerce to character or numeric
    chromosome_name = str_remove(seqname, "^chr"),
    start,
    end,
    rank = row_number()              # auto 1,2,3...
  )

# 2) TCGA-AI table

TCGA_TP53_AI <- TCGA_AI_table_final %>%
  filter(gene=="TP53") %>%
  filter(SNV_varity != "stopgain")

# make sure 'pos', 'DS_max_abs' and 'Pangolin_abs_max' are numeric
TCGA_TP53_AI$pos               <- as.numeric(TCGA_TP53_AI$pos)
TCGA_TP53_AI$DS_max_abs        <- as.numeric(TCGA_TP53_AI$DS_max_abs)
TCGA_TP53_AI$Pangolin_abs_max  <- as.numeric(TCGA_TP53_AI$Pangolin_abs_max)

# 3) Create exon blocks and remove introns

gap <- 20  # how much space between exons
# Build exon_map with gaps
exon_map <- tp53_main %>%
  filter(!rank %in% c(1, 11)) %>%
  arrange(rank) %>%
  mutate(
    exon_length = end - start,
    # cum_start: first exon at 0, then sum of (length + gap) of all prior exons
    cum_start = c(0, head(cumsum(exon_length + gap), -1)),
    cum_end   = cum_start + exon_length
  ) %>%
  select(rank, start, end, exon_length, cum_start, cum_end)

# Compress exon axis
TCGA_TP53_AI <- TCGA_TP53_AI %>%
  rowwise() %>%
  mutate(
    hits      = list(which(pos >= exon_map$start & pos <= exon_map$end)),
    exon_rank = if (length(hits)==1) exon_map$rank[hits]      else NA_integer_,
    exon_start= if (length(hits)==1) exon_map$start[hits]     else NA_integer_,
    exon_end   = if (length(hits)==1) exon_map$end[hits]       else NA_integer_,
    exon_cum  = if (length(hits)==1) exon_map$cum_start[hits] else NA_integer_
  ) %>%
  ungroup() %>%
  filter(!is.na(exon_rank)) %>%
  mutate(
    # offset       = pos - exon_start,
    offset       = exon_end - pos,
    x_compressed = exon_cum + offset
  )

# 4) Known codon positions
highlight_codons <- c("125","224","331")
transcript_id   <- "ENST00000269305"

TCGA_TP53_AI <- TCGA_TP53_AI %>%
  mutate(
    annots      = str_split(AAChange, ","),
    annots_tx   = lapply(annots, function(v) v[str_detect(v, transcript_id)]),
    codon_list  = lapply(annots_tx, function(v) {
      unlist(str_extract_all(v, "(?<=:p\\.[A-Z])\\d+(?=[A-Z])"))
    }),
    highlight   = sapply(codon_list, function(x) any(x %in% highlight_codons)),
    label       = sapply(codon_list, function(x) paste(intersect(x, highlight_codons),
                                                       collapse = ","))
  ) %>%
  dplyr::select(-annots, -annots_tx)

df_hl <- filter(TCGA_TP53_AI, highlight)

codon_hl_lines <- TCGA_TP53_AI %>% 
  filter(highlight) %>% 
  distinct(codon = label, pos) %>% 
  # drop any empty labels, if they snuck in
  filter(codon != "")

# pick a color for each codon
codon_cols <- c("125" = "#D73027",  # red
                "224" = "#4575B4",  # blue
                "331" = "#91BFDB")  # light‐blue

# Remap codon positions onto the compressed exon axis
codon_hl_mapped <- codon_hl_lines %>%
  rowwise() %>%
  mutate(
    hit        = list(which(pos >= exon_map$start & pos <= exon_map$end)),
    exon_start = if (length(hit)==1) exon_map$start[hit]     else NA_integer_,
    exon_end   = if (length(hit)==1) exon_map$end[hit]       else NA_integer_,
    exon_cum   = if (length(hit)==1) exon_map$cum_start[hit] else NA_integer_,
    # flip inside each exon just like for the SNVs
    offset       = if (!is.na(exon_start)) exon_end - pos  else NA_real_,
    x_compressed = if (!is.na(offset)) exon_cum + offset else NA_real_
  ) %>%
  ungroup() %>%
  filter(!is.na(x_compressed))

# 5) Recurrent mutations
TCGA_TP53_AI <- TCGA_TP53_AI %>%
  mutate(
    codon_str = map_chr(codon_list, ~ if(length(.x) >= 1) .x[[1]] else NA_character_)
  ) %>%  # group_by(pos) %>%
  group_by(codon_str) %>%
  mutate(
    recurrent = n_distinct(sample) >= 5
  ) %>%
  mutate(recurrent = as.integer(recurrent)) %>%  # FALSE→0, TRUE→1
  ungroup()

table(TCGA_TP53_AI$recurrent)

# 5) Plot

# pivot to long so we can facet by “metric” if you like
TCGA_TP53_AI_melt <- TCGA_TP53_AI %>%
  mutate(reg_AI_type = recode_factor(reg_AI_type,
    "Pos reg-AI" = "Pos",
    "No reg-AI"     = "No",
    "Neg reg-AI" = "Neg"
  )) %>%
  pivot_longer(
    cols      = c(DS_max_abs, Pangolin_abs_max),
    names_to  = "metric",
    values_to = "delta_score"
  )

codon_cols <- c(
  "125" = "#E41A1C",
  "224" = "#377EB8",
  "331" = "#4DAF4A"
)

final_plot <- ggplot() +
  # 1) exon blocks
  geom_rect(
    data = exon_map,
    aes(xmin = cum_start, xmax = cum_end, ymin = -Inf, ymax = Inf),
    fill  = "grey80", alpha = 0.5
  ) +
  # 2) codon v‐lines underneath
  geom_vline(
    data        = codon_hl_mapped,
    aes(xintercept = x_compressed),
    colour      = codon_cols[codon_hl_mapped$codon],
    linetype    = "dashed",
    size        = 0.6,
    alpha       = 0.7,
    show.legend = FALSE
  ) +
  # 3) add codon labels at y=0.5
  geom_text(
    data        = codon_hl_mapped,
    aes(x         = x_compressed+35,
        label     = codon),
    y           = 0.37,
    colour      = codon_cols[codon_hl_mapped$codon],
    size        = 3,
    vjust       = -0.5
  ) +
  # 4) now draw the points on top
  geom_point(
    data = TCGA_TP53_AI_melt,
    aes(
      x      = x_compressed,
      y      = delta_score,
      colour = reg_AI_type,
      # colour = factor(metric),
      shape  = reg_AI_type,
      # shape = factor(metric),
      size   = recurrent
    ),
    # size  = 1.5, 
    alpha = 0.5
  ) +
  # 5) your AI colour/shape scales
  scale_colour_manual(
    name   = "sig AI",
    values = c("Pos"="#00A087FF","No"="grey","Neg"="#3C5488FF")
  ) +
  scale_shape_manual(
    name   = "sig AI",
    values = c("Pos"=19,"No"=4,"Neg"=25)
  ) +
  # define the two point sizes
  # scale_size_manual(
  #   name   = "Mutation type",
  #   values = c(
  #     "Low recurrence" = 1,
  #     "High recurrence" = 3
  #   )
  # ) +
  scale_size_continuous(
    name  = "Recurrent\nmutated codon?",
    range = c(1.5, 3),         # size 1 when 0 (singleton), size 3 when 1 (recurrent)
    breaks= c(0, 1),
    labels= c("No", "Yes")
  ) +
  # 6) x‐axis exon labels
  scale_x_continuous(
    breaks = exon_map$cum_start + exon_map$exon_length/2,
    labels = paste0("E", exon_map$rank)
  ) +
  labs(
    x     = "Exons",
    y     = "Splicing delta score",
    title = "TP53"
  ) +
  theme_classic() +
  theme(
    plot.title      = element_text(hjust=0.5,size=14),
    axis.title      = element_text(size=10),
    axis.text       = element_text(size=8),
    legend.position = "top",
    legend.title    = element_text(size=9),
    legend.text     = element_text(size=8),
    legend.box      = "horizontal"
  ) +
  guides(
    shape = guide_legend(
      override.aes = list(size = 3),   # legend point size
      nrow         = 2,
      byrow        = TRUE
    )
  )

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/test.png"
ggsave(final_figure_path, final_plot, width = 275, height = 90, units = "mm") 

plot_list <- list(plot_df = TCGA_TP53_AI_melt, exon_coords = exon_map, codon_coords = codon_hl_mapped)

# Save
# write.table(plot_list, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig1/panel_A.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(plot_list, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig1/panel_A.RData")

# For the Manuscript

recurrent_not_known <- TCGA_TP53_AI_melt %>%
    filter(exon_rank == 5) %>%
    filter(delta_score > 0.4) %>% data.frame()
recurrent_not_known[,]

df <- TCGA_TP53_AI %>%
  filter(codon_str %in% 187) %>% data.frame()
df
df[1:5,1:15]
table(df$pos,df$SNV_varity)

table(TCGA_TP53_AI$reg_AI_type, TCGA_TP53_AI$CNA_AI_type)

#######################################################
################### SUPP. FIG S1B #####################
#######################################################

# Splicing enrichment analysis

plot_df <- TCGA_AI_table_final %>%
      filter(SNV_varity != "stopgain")
plot_df$reg_AI_type <- gsub(" reg-AI", "",plot_df$reg_AI_type)

gene_cognate_types <- names(table(TCGA_AI_table_final$gene_cognate))
gene_cognate_types <- gene_cognate_types[!gene_cognate_types %in% c("both","both_cognate")]
SNV_types <- names(table(TCGA_AI_table_final$SNV_varity))
SNV_types <- c(SNV_types,"Low Impact Missense","Synonymous","All SNVs")
SNV_types <- SNV_types[SNV_types != "stopgain"]
reg_AI_cat <- c("Pos","Neg")
print(gene_cognate_types)
print(SNV_types)
print(reg_AI_cat)
variables_columns <- c("DS_AG","DS_AL","DS_DG","DS_DL","DS_max_abs","larg_inc","larg_dec","Pangolin_abs_max")

all_df_fisher_test_res <- c()

for (splicing_thres in c(0.05,0.1,0.15)) {
  for (SNV_type_char in SNV_types) {
      print(SNV_type_char)
      # SNV type
      if ( SNV_type_char %in% c("effectively_syn","nonsynonymous SNV")) {
        df_test1 <- plot_df %>%
            filter(SNV_varity == SNV_type_char)
      } else if (SNV_type_char %in% c("Low Impact Missense","Synonymous")) {
        df_test1 <- plot_df %>%
            filter(SNV_varity_new == SNV_type_char)        
      } else if (SNV_type_char == "All SNVs") {
        df_test1 <- plot_df
      }
      
      for (gene_type_char in gene_cognate_types) {
          print(gene_type_char)
          # Gene type
          df_test2 <- df_test1 %>%
              filter(gene_cognate %in% c(gene_type_char))
          
          for (reg_AI_type_char in reg_AI_cat) {
              print(reg_AI_type_char)
              # reg-AI type
              if (reg_AI_type_char == "Pos") {
                  df_test3 <- df_test2 %>%
                      filter(reg_AI_type %in% c("Pos","No"))
              } else if (reg_AI_type_char == "Neg") {
                  df_test3 <- df_test2 %>%
                      filter(reg_AI_type %in% c("Neg","No"))
              }

              # print(tail(sort(table(df_test3$Gene.Symbol)),3))
              if (nrow(df_test3) == 0) {next}
              # next

              for (col_var in variables_columns) {
                  # print(paste0("<--------",col_var,"-------->"))
                  df_fisher_test_res <- data.frame(OR = NA, p_value = NA, CI_95_low = NA, CI_95_high = NA,
                                              CI_90_low = NA, CI_90_high = NA, CI_80_low = NA, CI_80_high = NA,
                                                  SNV_type = SNV_type_char, gene_type = gene_type_char,
                                                  reg_AI_cat = reg_AI_type_char, col_var = col_var,
                                                  splicing_thres = splicing_thres)
                  # Column variable type
                  if ( col_var %in% c("DS_AG","DS_AL","DS_DG","DS_DL","DS_max_abs") ) { 
                      # Splicing
                      df_test3 <- df_test3 %>%
                          mutate(variant_class = if_else(!!sym(col_var) >= splicing_thres, "High", "Low") )
                  } else if ( col_var %in% c("larg_inc","larg_dec","Pangolin_abs_max")) {
                      # Splicing
                      df_test3 <- df_test3 %>%
                          mutate(variant_class = if_else(abs(!!sym(col_var)) >= splicing_thres, "High", "Low") )
                  }

                  df_test3$reg_AI_type <- factor(df_test3$reg_AI_type, levels = c("No",reg_AI_type_char))
                  df_test3$variant_class <- factor(df_test3$variant_class, levels = c("Low","High"))
                  # Fisher Test comparing proportions (pos/neg vs no cis-ASE) in X gene type
                  # Test is: ([X] variant class: High / Low) / (cis-ASE: pos/neg / no)
                  df_test_table <- table(df_test3$variant_class, df_test3$reg_AI_type)
                  # Use pseudocount of 1 (NO)
                  df_test_table <- (df_test_table)
                  fisher_test_res <- fisher.test(df_test_table , conf.int = TRUE, conf.level = 0.95,
                                              alternative = "two.sided")
                  df_fisher_test_res$OR <- as.numeric(fisher_test_res$estimate)
                  df_fisher_test_res$p_value <- as.numeric(fisher_test_res$p.value)
                  # CI
                  df_fisher_test_res$CI_95_low <- as.numeric(fisher_test_res$conf.int)[1]
                  df_fisher_test_res$CI_95_high <- as.numeric(fisher_test_res$conf.int)[2]
                  fisher_test_res <- fisher.test(df_test_table , conf.int = TRUE, conf.level = 0.90,alternative = "two.sided")
                  df_fisher_test_res$CI_90_low <- as.numeric(fisher_test_res$conf.int)[1]
                  df_fisher_test_res$CI_90_high <- as.numeric(fisher_test_res$conf.int)[2]
                  fisher_test_res <- fisher.test(df_test_table , conf.int = TRUE, conf.level = 0.80,alternative = "two.sided")
                  df_fisher_test_res$CI_80_low <- as.numeric(fisher_test_res$conf.int)[1]
                  df_fisher_test_res$CI_80_high <- as.numeric(fisher_test_res$conf.int)[2]
                  # Sample size
                  df_fisher_test_res$n_yes_reg_AI_High_variable <- df_test_table["High",reg_AI_type_char]
                  df_fisher_test_res$n_yes_reg_AI_Low_variable <- df_test_table["Low",reg_AI_type_char]
                  df_fisher_test_res$n_no_reg_AI_High_variable <- df_test_table["High","No"]
                  df_fisher_test_res$n_no_reg_AI_Low_variable <- df_test_table["Low","No"]
                  # Save
                  all_df_fisher_test_res <- rbind(all_df_fisher_test_res, df_fisher_test_res)
              }
          }
      }
  }
}

head(all_df_fisher_test_res)

change_names <- function(df) {
  df %>%
    # 1) rename columns
    rename(
      SNV_varity   = SNV_type,
      gene_cognate = gene_type
    ) %>%
    # 2) recode everything in one go
    mutate(
      # gene_cognate relabel & factor‐level order
      gene_cognate = recode(
        gene_cognate,
        random       = "Passenger",
        essential    = "Essential",
        og_cognate   = "Cognate OG",
        oncogene     = "Noncognate OG",
        tsg          = "Noncognate TSG",
        tsg_cognate  = "Cognate TSG",
        .default     = as.character(gene_cognate)
      ),
      gene_cognate = factor(
        gene_cognate,
        levels = c("Passenger","Essential",
                   "Noncognate OG","Cognate OG",
                   "Noncognate TSG","Cognate TSG",
                   "TP53")
      ),

      # SNV_varity relabel
      SNV_varity = recode(
        SNV_varity,
        effectively_syn       = "Eff. synonymous",
        `nonsynonymous SNV`   = "Missense",
        stopgain              = "Nonsense",
        .default              = as.character(SNV_varity)
      ),

      # col_var relabel
      col_var = recode(
        col_var,
        DS_AG            = "SP_AG",
        DS_AL            = "SP_AL",
        DS_DG            = "SP_DG",
        DS_DL            = "SP_DL",
        DS_max_abs       = "SP_MaxDelta",
        larg_inc         = "PG_Gain",
        larg_dec         = "PG_Loss",
        Pangolin_abs_max = "PG_MaxDelta",
        .default         = as.character(col_var)
      ),
      col_var = factor(
        col_var,
        levels = c("SP_MaxDelta","SP_AG","SP_AL","SP_DG","SP_DL","PG_MaxDelta","PG_Gain","PG_Loss")
      ),
      # reg_AI_cat relabel & factor-levels
      reg_AI_cat = recode(
        reg_AI_cat,
        Pos = "pos reg-AI",
        Neg = "neg reg-AI",
        .default = as.character(reg_AI_cat)
      ),
      reg_AI_cat = factor(
        reg_AI_cat,
        levels = c("pos reg-AI","neg reg-AI")
      )
    )
}

# Change name of variables
all_df_fisher_test_res_filt <- change_names(all_df_fisher_test_res)

# Remove Inf values
all_df_fisher_test_res_filt <- all_df_fisher_test_res_filt %>% 
  filter(OR != Inf & OR != -Inf & OR != 0)

# ORs FDR correction
all_df_fisher_test_res_filt <- all_df_fisher_test_res_filt %>%
  group_by(splicing_thres,reg_AI_cat,SNV_varity) %>%
  mutate(p_value_FDR_adjusted = p.adjust(p_value, method = "fdr"))

# Add labels with conditional coloring
all_df_fisher_test_res_filt <- all_df_fisher_test_res_filt %>%
    mutate(
        # < 0.15 --> "*"; < 0.05 --> "**"; < 0.01 --> "***"
        label = case_when(
            p_value_FDR_adjusted > 0.15 ~ "",
            p_value_FDR_adjusted > 0.05 ~ "*",
            p_value_FDR_adjusted > 0.01 ~ "**",
            !is.na(p_value_FDR_adjusted) ~ "***",
            TRUE ~ NA_character_
        )
    )

plot_df <- all_df_fisher_test_res_filt %>%
        filter(
          (splicing_thres == 0.1  & str_detect(col_var, "^PG")) |
          (splicing_thres == 0.15 & str_detect(col_var, "^SP"))
              ) #%>%
        # filter(!SNV_varity %in% c("Nonsense","All SNVs")) %>%
        # filter(!gene_cognate %in% c("Essential"))

plot <- ggplot(plot_df, aes(x = col_var, y = gene_cognate)) + 
    geom_tile(aes(fill = log2(OR))) +
    scale_fill_gradient2(low = "#0D5EAF", mid = "white", high = "#AF540D", midpoint = 0,
                        #  limits = c(-2.2, 3.65), breaks = seq(-2.2, 3.65, by = 1), 
                        name = expression(log[2](`OR`))) +
    geom_text(aes(label = label), nudge_y = 0, size = 3) +
    scale_color_identity() +  # Ensures colors are taken directly from 'label_color'
    facet_grid(reg_AI_cat ~ SNV_varity) +
    labs(x = "", y = "") +
    theme_classic()+
    theme(
        legend.position = "right",
        axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(size = 9),
        plot.margin = unit(c(0, 1, 0, 1), "cm")
    )

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
# ggsave(final_figure_path, plot, width = 150, height = 100, units = "mm") 

# Save
# write.table(plot_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig2/Fig2D.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# saveRDS(plot_df, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig2/Fig2D.RData")

# Save
write.table(all_df_fisher_test_res_filt, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig1/panel_B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(all_df_fisher_test_res_filt, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig1/panel_B.RData")

#######################################################
################### SUPP. FIG S1C #####################
#######################################################

# Splicing enrichment analysis --> only SpliceAI and Pangolin MaxDelta

calculate_bd_p <- function(df) {
  # must have exactly Passenger + one other
  if(nrow(df) != 2) return(NA_real_)
  other <- setdiff(df$gene_cognate, "Passenger") 
  # reorder Passenger first
  df <- df %>% arrange(match(gene_cognate, c("Passenger", other)))
  # pull out & coerce to numeric
  h_yes <- as.numeric(df$n_yes_reg_AI_High_variable[1:2])
  h_no  <- as.numeric(df$n_no_reg_AI_High_variable[1:2])
  l_yes <- as.numeric(df$n_yes_reg_AI_Low_variable[1:2])
  l_no  <- as.numeric(df$n_no_reg_AI_Low_variable[1:2])
  # bail on any missing
  if(any(is.na(c(h_yes, h_no, l_yes, l_no)))) return(NA_real_)
  # build 2×2×2 array as numeric
  vals <- c(
    # slice 1 = Passenger
    h_yes[1], h_no[1], l_yes[1], l_no[1],
    # slice 2 = other class
    h_yes[2], h_no[2], l_yes[2], l_no[2]
  )
  tbl3 <- array(
    vals,
    dim   = c(2, 2, 2),
    dimnames = list(
      SplicingImpact   = c("High","Low"),
      AlleleImbalance  = c("Yes","No"),
      GeneSet          = c("Passenger", other)
    )
  )
  # run test quietly, catch errors
  pval <- tryCatch({
    suppressWarnings(
      BreslowDayTest(tbl3, correct = FALSE)$p.value
    )
  }, error = function(e) NA_real_)
  pval
}

all_df_bd_res <- map_dfr(
  c("Cognate OG","Noncognate OG","Cognate TSG","Noncognate TSG","Essential","TP53"),
  function(gene_class){
    all_df_fisher_test_res_filt %>%
      filter(gene_cognate %in% c(gene_class, "Passenger")) %>%
      group_by(col_var, SNV_varity, reg_AI_cat, splicing_thres) %>%
      nest() %>%
      mutate(
        breslow_day_p = map_dbl(data, calculate_bd_p),
        compare_to    = gene_class
      ) %>%
      select(-data)
  }
)

df_plot <- all_df_fisher_test_res_filt %>%
  # keep exactly the two MaxDelta tests you want
  filter(
    (col_var == "PG_MaxDelta" & splicing_thres == 0.1) |
    (col_var == "SP_MaxDelta" & splicing_thres == 0.15)
  ) %>%
  # bring in your Breslow‐Day p
  left_join(
    all_df_bd_res %>%
      select(SNV_varity, reg_AI_cat, col_var, splicing_thres,
             compare_to, breslow_day_p),
    by = c(
      "SNV_varity", "reg_AI_cat", "col_var", "splicing_thres",
      "gene_cognate" = "compare_to"
    )
  ) %>%
  # for each gene+facet, compute Fisher’s combined p
  group_by(gene_cognate, SNV_varity, reg_AI_cat) %>%
  mutate(
    # only two rows per group, so this sums exactly those two p’s
    fisher_X2 = -2 * sum(log(breslow_day_p)),
    fisher_df = 2 * n(),  
    fisher_p  = pchisq(fisher_X2, df = fisher_df, lower.tail = FALSE)
  ) %>%
  ungroup() %>%
  # drop the fisher_p from the PG rows so only SP shows it
  mutate(
    fisher_p = ifelse(col_var == "SP_MaxDelta", fisher_p, NA_real_)
  )

df_plot$gene_cognate <- factor(df_plot$gene_cognate, 
      levels = levels(all_df_fisher_test_res_filt$gene_cognate))

plot_levels <- levels(df_plot$gene_cognate)

# 2) build annot_df with a dynamic y based on rank
annot_df <- df_plot %>%
  filter(gene_cognate != "Passenger") %>%
  distinct(reg_AI_cat, SNV_varity, gene_cognate, fisher_p) %>%
  mutate(
    gene_cognate = factor(gene_cognate, levels = plot_levels),
    xstart       = as.numeric(factor("Passenger", levels = plot_levels)),
    xend         = as.numeric(gene_cognate)
  ) %>%
  group_by(reg_AI_cat, SNV_varity) %>%
  arrange(xend, .by_group = TRUE) %>%
  mutate(
    rank = row_number(),
    y0   = 2 + (rank - 1) * 0.25,
    # only label the significant ones
    label = ifelse(!is.na(fisher_p) & fisher_p < 0.05,
                   paste0("p=", format(fisher_p, scientific = TRUE, digits = 2)),
                   NA_character_),
    # only keep y where you have a label
    y = ifelse(!is.na(label), y0, NA_real_)
  ) %>%
  ungroup()

# 3) re‐plot with these staggered y’s
plot <- ggplot(df_plot,
            aes(x     = gene_cognate,
                y     = log2(OR),
                color = gene_cognate,
                shape = col_var,
                group = interaction(gene_cognate, col_var))) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = log2(CI_90_low),
                    ymax = log2(CI_90_high)),
                width    = 0.3,
                position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(reg_AI_cat ~ SNV_varity, scales = "fixed") +
  coord_cartesian(ylim = c(-3.5, 5.5)) +   # extend top margin
  # coord_cartesian(ylim = c(-3, max(annot_df$y) + 0.5)) +   # extend top margin
  scale_color_manual(values = c(
    "Passenger"   = "#83c4be",
    "Cognate OG"  = "#37d576",
    "Cognate TSG" = "#7083e5",
    "TP53"        = "#7cd8ef",
    "Essential"   = "#db4538"
  )) +
  scale_shape_manual(values = c(
    "PG_MaxDelta" = 16,
    "SP_MaxDelta" = 17
  )) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    strip.text      = element_text(size = 7.5),
    plot.title      = element_text(size = 10)
  ) +
  labs(
    x     = NULL,
    y     = expression(log[2](OR)),
    shape = "MaxDelta / thres",
    title = "PG (0.1) vs SP (0.15) with staggered BD p-values"
  ) +
  geom_segment(data        = annot_df,
               aes(x    = xstart,
                   xend = xend,
                   y    = y-0.2,
                   yend = y-0.2),
               inherit.aes = FALSE) +
  geom_text(data        = annot_df,
            aes(x     = (xstart + xend)/2,
                y     = y + 0.1,
                label = label),
            inherit.aes = FALSE,
            size = 3)

# plot_list <- list(plot_df = df_plot, annot_df = annot_df)

# Save
# write.table(df_plot, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig2/Fig2E.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# saveRDS(df_plot, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig2/Fig2E.RData")
write.table(df_plot, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig1/panel_C.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(df_plot, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig1/panel_C.RData")

#######################################################
################### SUPP. FIG S1D #####################
#######################################################

plot_df <- TCGA_AI_table_final
plot_df$reg_AI_type <- gsub(" reg-AI", "",plot_df$reg_AI_type)
plot_df$Gene.Symbol <- as.character(plot_df$Gene.Symbol)

all_df_fisher_test_res_genes <- c()

genes_to_test <- plot_df %>%
  filter(gene %in% c("both","TP53","tsg","oncogene")) %>%
  pull(Gene.Symbol) %>% as.character %>% unique()

for (SNV_type_char in SNV_types) {
    print(SNV_type_char)
    # SNV type
    if ( SNV_type_char %in% c("effectively_syn","nonsynonymous SNV","stopgain")) {
        df_test1 <- plot_df %>%
            filter(SNV_varity == SNV_type_char)
    } else if (SNV_type_char %in% c("Low Impact Missense","Synonymous")) {
        df_test1 <- plot_df %>%
            filter(SNV_varity_new == SNV_type_char)        
    }

    for (reg_AI_type_char in reg_AI_cat) {
      print(reg_AI_type_char)
        n <- 0
        results <- df_test1 %>%
            group_by(Gene.Symbol) %>%
            group_modify(~ {
                df <- .x
                n <<- n + 1
                timing <- round((n/length(genes_to_test))*100,2)
                print(paste0("<<<<<<<<<------- ",timing," ------->>>>>>>>>>>"))

                # reg-AI type
                if (reg_AI_type_char == "Pos") {
                    df_test2 <- df %>%
                        filter(reg_AI_type %in% c("Pos","No"))
                } else if (reg_AI_type_char == "Neg") {
                    df_test2 <- df %>%
                        filter(reg_AI_type %in% c("Neg","No"))
                }
                # df <- df %>% 
                #     filter(cognate == TRUE)
                # print(table(df$cognate))
                
                if (nrow(df_test2) == 0) { return(tibble())}
                
                results_by_var <- map_dfr(variables_columns, function(col_var) {

                  # Column variable type
                  if ( col_var %in% c("DS_AG","DS_AL","DS_DG","DS_DL","DS_max_abs") ) { 
                      # Splicing
                      df_test2 <- df_test2 %>%
                          mutate(variant_class = if_else(!!sym(col_var) >= 0.15, "High", "Low") )
                  } else if ( col_var %in% c("larg_inc","larg_dec","Pangolin_abs_max")) {
                      # Splicing
                      df_test2 <- df_test2 %>%
                          mutate(variant_class = if_else(!!sym(col_var) >= 0.10, "High", "Low") )
                  }
                  # print(colnames(df_test2))

                  df_test_table <- table(df_test2$variant_class, df_test2$reg_AI_type) # Check dimensions only

                  df_res_table <- tibble(
                      OR = NA_real_,
                      p_value = NA_real_,
                      CI_95_low = NA_real_,
                      CI_95_high = NA_real_,
                      CI_90_low = NA_real_,
                      CI_90_high = NA_real_,
                      CI_80_low = NA_real_,
                      CI_80_high = NA_real_,
                      reg_AI_cat = reg_AI_type_char,
                      col_var = col_var,
                      n_yes_reg_AI_High_variable = NA_real_,
                      n_yes_reg_AI_Low_variable = NA_real_,
                      n_no_reg_AI_High_variable = NA_real_,
                      n_no_reg_AI_Low_variable = NA_real_,
                      message = "Insufficient data for Fisher's test"
                  )
                  # Check the dimensions of the table before performing Fisher's test
                  if (any(dim(df_test_table) < 2)) {
                      return(df_res_table)
                  } else if (any(dim(df_test_table) >= 2)) {
                      df_test2$reg_AI_type <- factor(df_test2$reg_AI_type, levels = c("No",reg_AI_type_char))
                      df_test2$variant_class <- factor(df_test2$variant_class, levels = c("Low","High"))
                      df_test_table <- table(df_test2$variant_class, df_test2$reg_AI_type)
                      fisher_test_res <- fisher.test(df_test_table , conf.int = TRUE, alternative = "two.sided")
                      # print(df_test_table)
                      # print(fisher_test_res)
                      # Save stuff
                      df_res_table["OR"] <- as.numeric(fisher_test_res$estimate)
                      df_res_table$p_value <- as.numeric(fisher_test_res$p.value)
                      # CI other levels
                      df_res_table$CI_95_low <- as.numeric(fisher_test_res$conf.int[1])
                      df_res_table$CI_95_high <- as.numeric(fisher_test_res$conf.int[2])
                      fisher_test_res <- fisher.test(df_test_table , conf.int = TRUE, conf.level = 0.90, alternative = "two.sided")
                      df_res_table$CI_90_low <- as.numeric(fisher_test_res$conf.int)[1]
                      df_res_table$CI_90_high <- as.numeric(fisher_test_res$conf.int)[2]
                      fisher_test_res <- fisher.test(df_test_table , conf.int = TRUE, conf.level = 0.80,alternative = "two.sided")
                      df_res_table$CI_80_low <- as.numeric(fisher_test_res$conf.int)[1]
                      df_res_table$CI_80_high <- as.numeric(fisher_test_res$conf.int)[2]
                      # Other
                      df_res_table$reg_AI_cat <- reg_AI_type_char
                      df_res_table$col_var <- col_var
                      df_res_table$message <- "Test completed"
                      df_res_table$n_yes_reg_AI_High_variable <- df_test_table["High",reg_AI_type_char]
                      df_res_table$n_yes_reg_AI_Low_variable <- df_test_table["Low",reg_AI_type_char]
                      df_res_table$n_no_reg_AI_High_variable <- df_test_table["High","No"]
                      df_res_table$n_no_reg_AI_Low_variable <- df_test_table["Low","No"]
                      # Return
                      return(df_res_table)
                  }
                })
                return(results_by_var)
            }) %>%
            ungroup()
        results$reg_AI_cat <- reg_AI_type_char
        results$SNV_type <- SNV_type_char
        all_df_fisher_test_res_genes <- rbind(all_df_fisher_test_res_genes,results)
    } #reg_AI_type
} #SNV type

# View the results
print(all_df_fisher_test_res_genes)

# Save
output_path <- "/g/strcombio/fsupek_cancer1/gpalou/ASE_project/Pangolin/splicing_enrichment_cancer_genes_reg_AI_2.txt"
# write.table(file = output_path,
#             all_df_fisher_test_res_genes,
#             sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

all_df_fisher_test_res_genes <- read.csv(output_path, header = TRUE, sep = "\t")
dim(all_df_fisher_test_res_genes)

colnames(all_df_fisher_test_res_genes)[1] <- "gene_symbol"

# Filters
df_clean <- all_df_fisher_test_res_genes %>%
  # 1) filter out bad ORs
  filter(
    !is.infinite(OR),      # no ±Inf
    OR > 0,                # drop zeros or negatives
    OR <= 25
  ) %>%
  # 2) bring in gene type
  # 2) bring in gene *type* under its own name
  left_join(
    TCGA_AI_table_final %>%
      dplyr::select(
        gene_symbol = Gene.Symbol,
        # gene_type   = gene      # rename here
        gene,
      ) %>%
      dplyr::distinct(),
    by = "gene_symbol"
  ) %>%
  # 3) drop unwanted 'gene' categories and SNV_type
  filter(
    !gene %in% c("random", "essential", "both"),
    SNV_type != "stopgain"
  ) %>%
  # 4) correct p‐values FDR within each group
  group_by(reg_AI_cat, SNV_type, col_var) %>%
  mutate(
    p_value_FDR = p.adjust(p_value, method = "fdr")
  ) %>%
  ungroup() %>%
  # 5) recode factors and add significance flag
  mutate(
    SNV_type = recode(
      SNV_type,
      effectively_syn       = "Eff. Syn.",
      `nonsynonymous SNV`   = "High aa-impact Miss.",
      `Low Impact Missense` = "Low aa-impact Miss.",
      .default              = SNV_type
    ),
    reg_AI_cat = recode(
      reg_AI_cat,
      Pos = "Pos reg-AI",
      Neg = "Neg reg-AI"
    ),
    gene = recode(
      gene,
      TP53 = "TSG",     # collapse TP53 → tumor suppressor
      tsg  = "TSG",
      .default = "OG"   # everything else as oncogene
    ),
    significant = case_when(
      p_value_FDR <= 0.10 ~ "<=10%",
      p_value_FDR <= 0.15 ~ "<=15%",
      p_value_FDR <= 0.25 ~ "<=25%",
      TRUE                ~ "ns"
    )
  )

# If you then want the DS_max_abs slice sorted by raw p‐value:

df_best <- df_clean %>%
  # 1) keep only the two columns of interest
  filter(col_var %in% c("Pangolin_abs_max", "DS_max_abs")) %>%
  # 2) for each gene, pick the single row with lowest FDR
  group_by(gene_symbol) %>%
  slice_min(order_by = p_value_FDR, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(p_value_FDR)

# -log10(p_value)
list_plots <- list()

# Create a named color vector for the gene categories
gene_colors <- c("OG" = "#37d576", "TSG" = "#7083e5")

for (reg_AI_type in c("Neg reg-AI","Pos reg-AI")) {

        df_tmp_filt <- df_best %>%
                filter(reg_AI_cat == reg_AI_type) %>%
                # arrange(desc(OR))
                slice_min(order_by = p_value_FDR, n = 30, with_ties = FALSE)

        df_tmp_filt$ID <- paste0(df_tmp_filt$gene_symbol,"_",df_tmp_filt$SNV_type)
        df_tmp_filt$gene_symbol <- factor(
          df_tmp_filt$gene_symbol,
          levels = unique(as.character(df_tmp_filt$gene_symbol))
        )

        plot_tmp <- ggplot(data = df_tmp_filt, mapping = aes(x = gene_symbol, y = log2(OR) )) + 
                geom_label(aes(label = "", fill = gene), color = NA, size = 2) +  # Empty labels with no text and no border
                geom_point(aes(color = significant, shape = SNV_type), size = 4, alpha = 1) + 
                geom_errorbar(aes(ymin = log2(CI_90_low), ymax = log2(CI_90_high)), position = position_dodge(width=0.75), 
                        width = 0.2, color = "#000000", size = 0.5) +
                facet_grid(. ~ reg_AI_cat) + 
                geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
                scale_fill_manual(values = c("OG" = "#37d576", "TSG" = "#7083e5"), 
                                name = "Gene category") +  # Custom legend labels               
                scale_color_manual(values = c("<=10%" = "#24c1a9", "<=15%" = "#ad2a2a", 
                                  "<=25%" = "#9a59a8","ns" = "#C0C0C0")) +  # Custom colors
                scale_shape_manual(
                    values = c(
                        "Eff. Syn." = 16,
                        "High aa-impact Mis." = 17,
                        "Low aa-impact Miss." = 18,
                        "Synonymous" = 15
                    )) +
                labs(shape = "", color = "FDR", x = "", y = expression(log[2](OR))) +
                theme_classic() + 
                theme(
                        # strip.text = element_text(size = 10),
                        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, color = gene_colors[df_tmp_filt$gene]), # Color labels by gene type
                        legend.title = element_text(size = 9),
                        legend.text = element_text(size = 8),
                        legend.key.width = unit(0.7, "cm"),
                        legend.key.size = unit(0.05, "cm"),  # Smaller legend keys
                        legend.spacing.y = unit(0.005, "cm"),
                        legend.spacing.x = unit(0.5, "cm"),
                        legend.position = "bottom") +
                guides(shape = guide_legend(nrow = 2),
                        color = guide_legend(nrow = 2),
                        fill = guide_legend(nrow = 2))  # Applies 3-row layout for all legends

        list_plots[[reg_AI_type]] <- plot_tmp
}

plot <- ggarrange(list_plots[["Neg reg-AI"]], list_plots[["Pos reg-AI"]], align = "hv", widths = c(0.6,0.4),
                ncol=2, nrow=1, common.legend = TRUE, legend="top") + theme_classic() +
                    theme(
        plot.margin = unit(c(0, 3, 0, 3), "cm")
    )

# Save
# write.table(df_best, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig2/Fig2F.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# saveRDS(df_best, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig2/Fig2F.RData")
write.table(df_best, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig1/panel_D.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(df_best, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig1/panel_D.RData")

#For the manuscript
table(all_df_fisher_test_res_genes$SNV_type)
all_df_fisher_test_res_genes %>%
  filter(col_var == "DS_max_abs") %>%
  filter(SNV_type == "Low Impact Missense") %>%
  filter(reg_AI_cat == "Neg reg-AI") %>%
  filter(gene_symbol == "CDH1") %>% data.frame()

######################################################################################
######################### PUFFIN-D (PROMOTER ACTIVITY) ###############################
######################################################################################

#######################################################
################### SUPP. FIG S2A #####################
#######################################################

# A) For the Manuscript

TCGA_AI_table_final <- TCGA_AI_table_final %>%
  filter(SNV_varity != "stopgain")

# Add Positive control mutations (TERT/CDC20/FOXA1)
pos_control <- read.table(paste0("/g/strcombio/fsupek_fisher/gpalou/Puffin/batches/output/artificial_sequences/pos_control_puffin_output.txt"),
                    header = TRUE, sep = "\t")
pos_control <- pos_control %>%
  # switch to rowwise for c_across()
  # filter(Gene.Symbol == "TERT") %>%
  rowwise() %>%
  mutate(
    Puffin_MaxDelta = {
      values <- c_across(c(
        PRO_CAP_50bp_TSS_sum_all_abs,
        FANTOM_CAGE_50bp_TSS_sum_all_abs,
        GRO_CAP_50bp_TSS_sum_all_abs
      ))
      idx <- which.max(abs(values))
      if (length(idx) == 0) 
        NA_real_ 
      else 
        values[idx]
    }
  ) %>%
  ungroup()

options(scipen = 999)

# 1) Q2 threshold
puffin_thres_Q2 <- quantile(pos_control$Puffin_MaxDelta,seq(0,1,0.05))["50%"]
puffin_thres_Q2
# % in Passenger genes
df <- TCGA_AI_table_final %>%
  filter(gene_cognate == "random")
prop.table(table(df$Puffin_MaxDelta > puffin_thres_Q2))*100
# % in cancer genes
df <- TCGA_AI_table_final %>%
  filter(gene %in% c("both","oncogene","tsg","TP53"))
prop.table(table(df$Puffin_MaxDelta > puffin_thres_Q2))*100

# 2) Q1 threshold
puffin_thres_Q1 <- quantile(pos_control$Puffin_MaxDelta,seq(0,1,0.05))["25%"]
puffin_thres_Q1
# % in Passenger genes
df <- TCGA_AI_table_final %>%
  filter(gene_cognate == "random")
prop.table(table(df$Puffin_MaxDelta > puffin_thres_Q1))*100
# % in cancer genes
df <- TCGA_AI_table_final %>%
  filter(gene %in% c("both","oncogene","tsg","TP53"))
prop.table(table(df$Puffin_MaxDelta > puffin_thres_Q1))*100

# 3) Example genes

df_tmp <- TCGA_AI_table_final %>%
  filter(gene == "tsg") %>%
  filter(Puffin_MaxDelta > puffin_thres_Q2) %>%
  # filter(reg_AI_type == "Pos reg-AI") %>%
  # filter(Gene.Symbol == "KRAS") %>% 
  data.frame()
dim(df_tmp)
sort(table(as.character(df_tmp$Gene.Symbol)))
sort(df_tmp$Puffin_MaxDelta)
table(df_tmp$reg_AI_type)
df_tmp[,c("RNA_ref","RNA_alt","AAChange","reg_AI_type","CNA_AI_type","mRNA_AI_type", "SNV_varity",
      "GRO_CAP_50bp_TSS_sum_all","PRO_CAP_50bp_TSS_sum_all","FANTOM_CAGE_50bp_TSS_sum_all")]

df_tmp <- TCGA_AI_table_final %>%
  filter(gene %in% c("tsg","oncogene")) %>%
  filter(Puffin_MaxDelta > puffin_thres_Q1)
summary(df_tmp$TSS_MUT_pos_distance)
table(df_tmp$reg_AI_type)

tbl <- table(df_tmp$reg_AI_type)
tbl
( tbl["Neg reg-AI"]+tbl["Pos reg-AI"] ) / (tbl["Neg reg-AI"]+tbl["Pos reg-AI"]+tbl["No reg-AI"])*100

# B) Puffin-D delta score vs TSS position

TCGA_AI_table_final_filt <- TCGA_AI_table_final %>%
                        filter(!gene_cognate %in% c("both","both_cognate")) %>%
                        filter(SNV_varity != "stopgain")

# Remove duplicates
duplicated_variants <- which(duplicated(TCGA_AI_table_final_filt[,c("gene_cognate","SNV_varity","Gene.Symbol","TSS_MUT_pos_distance","Puffin_MaxDelta")]))
TCGA_AI_table_final_filt <- TCGA_AI_table_final_filt[-duplicated_variants,]

# Create bins of 10,000 base pairs, with the last bin including all values >40,000
TCGA_AI_table_final_filt$TSS_MUT_pos_distance_bins <- ifelse(
  is.na(TCGA_AI_table_final_filt$TSS_MUT_pos_distance), 
  NA,  # Keep NAs as NAs
  cut(
    abs(TCGA_AI_table_final_filt$TSS_MUT_pos_distance),  # Absolute values of TSS_MUT_pos_distance
    breaks = c(0, 1000, 2000, 3000, 4000, 5000, 10000, 20000, 30000, 40000, 50000, Inf),  # Breaks at 0, 10k, 20k, 30k, 40k, and the rest
    include.lowest = TRUE,  # Include the lowest interval
    right = FALSE,  # Intervals are closed on the left
  )
)

TCGA_AI_table_final_filt$TSS_MUT_pos_distance_bins_2 <- ifelse(
  is.na(TCGA_AI_table_final_filt$TSS_MUT_pos_distance), 
  NA,  # Keep NAs as NAs
  cut(
    abs(TCGA_AI_table_final_filt$TSS_MUT_pos_distance),  # Absolute values of TSS_MUT_pos_distance
    breaks = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, Inf),  # Breaks at 0, 10k, 20k, 30k, 40k, and the rest
    include.lowest = TRUE,  # Include the lowest interval
    right = FALSE,  # Intervals are closed on the left
  )
)
table(TCGA_AI_table_final_filt$TSS_MUT_pos_distance_bins,TCGA_AI_table_final_filt$gene_cognate)
table(TCGA_AI_table_final_filt$TSS_MUT_pos_distance_bins_2,TCGA_AI_table_final_filt$gene_cognate)
# Filt
TCGA_AI_table_final_filt <- TCGA_AI_table_final_filt[!is.na(TCGA_AI_table_final_filt$TSS_MUT_pos_distance_bins),]
table(TCGA_AI_table_final_filt$TSS_MUT_pos_distance_bins,TCGA_AI_table_final_filt$gene_cognate)

# Ensure gene_cognate is a factor and set levels to control the legend order
TCGA_AI_table_final_filt <- TCGA_AI_table_final_filt %>%
  mutate(gene_cognate = recode(gene_cognate,
                               "random" = "Passenger",
                               "og_cognate" = "Cognate OG",
                               "oncogene" = "Noncognate OG",
                               "tsg_cognate" = "Cognate TSG",
                               "tsg" = "Noncognate TSG",
                               "TP53" = "TP53",
                               "pos contrl - known" = "Pos. control")) %>%
  mutate(gene_cognate = factor(gene_cognate, 
                               levels = c("Pos. control", "Passenger", 
                                          "Cognate OG", "Noncognate OG", 
                                          "Cognate TSG", "Noncognate TSG", 
                                          "TP53")))

# Calculate the top 2 and bottom 2 for each boxplot group
df_labels <- TCGA_AI_table_final_filt %>%
  group_by(gene_cognate, TSS_MUT_pos_distance_bins) %>%
  arrange(desc(Puffin_MaxDelta)) %>%
  dplyr::slice(1) %>%  # Select top 2
  ungroup()

Puffin_D_threshold <- as.numeric(puffin_thres_Q2)

plot <- ggplot(TCGA_AI_table_final_filt, mapping = aes(
                        x = factor(TSS_MUT_pos_distance_bins), y = log10(abs(Puffin_MaxDelta)),
                        fill = factor(gene_cognate) 
                        ) ) + 
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = log10(Puffin_D_threshold), ymax = Inf, 
                 fill = "lightgreen", alpha = 0.2) +
        geom_point(alpha = 0.05, position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9), size = 0.5) +  # Add jitter to avoid overlap of points
        geom_boxplot(width = 1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
        geom_hline(yintercept = log10(Puffin_D_threshold), color = "red", linetype = "dashed", size = 0.5) +
        geom_text(data = df_labels, aes(label = Gene.Symbol), check_overlap = FALSE,
          position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9),  # Align text with point dodge
          hjust = 0.5, vjust = -0.5, size = 1, color = "black") +
        scale_x_discrete(labels = c("1kb", "2kb", "3kb", "4kb", "5kb", "10kb", "20kb", "30kb", "40kb", ">50kb")) +
        scale_fill_manual(
          values = c(
            "Passenger" = "#eef350",                 # Random
            "Cognate OG" = set1_colors[5],        # OG cognate
            "Noncognate OG" = scales::muted(set1_colors[5], l = 70),        # OG non-cognate
            "Cognate TSG" = set1_colors[2],       # TSG cognate
            "Noncognate TSG" = scales::muted(set1_colors[2], l = 70),  # TSG non-cognate, make lighter
            "TP53" = set1_colors[3],              # TP53
            "Pos. control" = set1_colors[1]
          )) +
        facet_grid(. ~ ., scales = "free_y") +
        theme_bw() + 
        theme(
                legend.position = "top",
                axis.text.x = element_text(angle = 90, size = 10, hjust = 1, vjust = 0.5),
                strip.text = element_text (size = 10)
        ) +
        labs(fill = "", x = "TSS - Mutation distance (bins)", y = "log10(abs(Puffin-D score))") +
        coord_cartesian(ylim = c(-5,1))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, plot, width = 175, height = 180, units = "mm") 

write.table(TCGA_AI_table_final_filt, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig2/panel_A_B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(TCGA_AI_table_final_filt, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig2/panel_A_B.RData")

#############################################################################################
######################### NON-CODING PROMOTER MUTATIONS (WGS) ###############################
#############################################################################################

######################################################
################### SUPP. FIG S3 #####################
######################################################

# 1) Open data
# 1.1) Non-coding near TSS
input_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_cancers_WGS/TCGA_pancancer_TSS_table.txt")
TCGA_WGS_pancancer_TSS <- read.table(file = input_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dim(TCGA_WGS_pancancer_TSS)

# 1.2) SVs
# input_path <- paste0("/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_cancers_WGS/TCGA_pancancer_SVs_table.txt")
# TCGA_WGS_pancancer_SVs <- read.table(file = input_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# dim(TCGA_WGS_pancancer_SVs)

# 2) Filters
# Main isoform available
TCGA_WGS_pancancer_TSS_filt <- TCGA_WGS_pancancer_TSS #%>%
    # filter(MANE_transcript == "MANE")
# Calculate the counts of genes with and without MANE transcripts
mane_counts <- TCGA_WGS_pancancer_TSS %>%
  group_by(ensembl_gene_id) %>%
  summarize(has_mane_transcript = any(!is.na(MANE_transcript))) %>%
  summarize(
    with_mane = sum(has_mane_transcript),
    without_mane = sum(!has_mane_transcript)
  )
mane_counts
# PASS
TCGA_WGS_pancancer_TSS_filt <- TCGA_WGS_pancancer_TSS_filt %>%
    filter(FILTER == "PASS")
# Quality?
# Coverage?
TCGA_WGS_pancancer_TSS_filt <- TCGA_WGS_pancancer_TSS_filt %>%
    filter(totalCount >= 10)
dim(TCGA_WGS_pancancer_TSS_filt)
# Heterozygous??

# 3) Add TCGA-AI data
cols <- c("chrom","pos","ref","alt","gene_id","sample","SNV_varity")
TCGA_AI_table_final_filt <- TCGA_AI_table_final[,cols]
TCGA_AI_table_final_filt$TSS_distance <- NA
cols <- c("CHROM","POS","REF","ALT","ensembl_gene_id","TCGA_sample","SNV_varity","TSS_distance")
colnames(TCGA_AI_table_final_filt) <- cols
# Filter for samples with WGS
dim(TCGA_AI_table_final_filt)
TCGA_AI_table_final_filt2 <- TCGA_AI_table_final_filt %>%
  filter(TCGA_sample %in% unique(TCGA_WGS_pancancer_TSS_filt$TCGA_sample))
dim(TCGA_AI_table_final_filt2)
# Merge
colnames(TCGA_WGS_pancancer_TSS_filt)[colnames(TCGA_WGS_pancancer_TSS_filt) %in% "SNV_type"] <- "SNV_varity" 
TCGA_WGS_pancancer_TSS_filt <- TCGA_WGS_pancancer_TSS_filt[,cols]
TCGA_pancancer_full_table <- rbind(TCGA_WGS_pancancer_TSS_filt,TCGA_AI_table_final_filt2)
dim(TCGA_pancancer_full_table)
head(TCGA_pancancer_full_table)

# Add Gene.Symbol, gene_cognate, reg_AI_type
cols2 <- c("chrom","pos","ref","alt","gene_id","sample")
TCGA_pancancer_full_table <- merge(TCGA_pancancer_full_table,
          TCGA_AI_table_final[,c(cols2,"reg_AI_type","CNA_AI_type","Gene.Symbol","gene_cognate","RNA_observed_ratio","RNA_expected_ratio_bb_model1")],
            by.x = c("CHROM","POS","REF","ALT","ensembl_gene_id","TCGA_sample"),
            by.y = cols2, all.x = TRUE)
dim(TCGA_pancancer_full_table)

# Filters
TCGA_AI_and_WGS_df <- TCGA_pancancer_full_table %>%
    filter(SNV_varity != "stopgain") %>%
    filter(!gene_cognate %in% c("both","both_cognate")) %>%
    filter(is.na(TSS_distance) | TSS_distance < 1500) %>%
    mutate(mutation_type = ifelse(nchar(REF) > 1 | nchar(ALT) > 1, "Indel", "SNV"))
table(TCGA_AI_and_WGS_df$mutation_type,TCGA_AI_and_WGS_df$SNV_varity)
head(TCGA_AI_and_WGS_df)
unique(TCGA_AI_and_WGS_df$TCGA_sample) # Only 777 WGS samples

# Co-occurrence for every gene_cognate and also for SNVs-level not gene-sample-level

library(dplyr)
library(stringr)

# 1) your exonic‐coding definition
coding_terms <- c("effectively_syn", "nonsynonymous SNV", "stopgain")

# 2) flag exonic vs non‐coding on every row
TCGA_flagged <- TCGA_AI_and_WGS_df %>%
  mutate(
    is_exonic = str_detect(SNV_varity, paste(coding_terms, collapse = "|"))
  )

# 3) find which (sample,gene) pairs have >=1 non‐coding SNV
has_noncoding <- TCGA_flagged %>%
  filter(!is_exonic) %>%
  distinct(TCGA_sample, ensembl_gene_id) %>%
  mutate(has_noncoding = TRUE)

# 4) attach that flag onto every exonic SNV row
exonic_with_flag <- TCGA_flagged %>%
  filter(is_exonic) %>%
  left_join(has_noncoding,
            by = c("TCGA_sample", "ensembl_gene_id")) %>%
  mutate(has_noncoding = coalesce(has_noncoding, FALSE)) %>%
  filter(!is.na(reg_AI_type))  # drop rows w/o an AI label

# 5a) summary by the **original** gene_cognate × AI  
orig_summary <- exonic_with_flag %>%
  group_by(gene_cognate, reg_AI_type) %>%
  summarise(
    n_exonic    = n(),
    n_cooccur   = sum(has_noncoding),
    pct_cooccur = 100 * n_cooccur / n_exonic
  ) %>%
  ungroup()

# 5b) summary _only_ for the cancer‐gene set, aggregated into one row per AI  
cancer_levels <- c("tsg_cognate", "tsg", "oncogene", "og_cognate","TP53")

cancer_summary <- exonic_with_flag %>%
  filter(gene_cognate %in% cancer_levels) %>%
  group_by(reg_AI_type) %>%
  summarise(
    n_exonic    = n(),
    n_cooccur   = sum(has_noncoding),
    pct_cooccur = 100 * n_cooccur / n_exonic
  ) %>%
  ungroup() %>%
  # tag it so it lines up with your orig_summary
  mutate(gene_cognate = "cancer_genes") %>%
  select(gene_cognate, everything())

# 6) stack them
final_summary <- bind_rows(orig_summary, cancer_summary) %>%
  arrange(gene_cognate, reg_AI_type)

data.frame(final_summary)

plot_df <- final_summary %>%
  mutate(
    non_cooccur = 100 - pct_cooccur
  ) %>%
  select(gene_cognate, reg_AI_type, pct_cooccur, non_cooccur) %>%
  pivot_longer(
    cols      = c(pct_cooccur, non_cooccur),
    names_to  = "outcome",
    values_to = "percent"
  ) %>%
  mutate(
    outcome = recode(outcome,
                     pct_cooccur = "Co-occur",
                     non_cooccur = "No co-occur"),
    # ensure the stacking order: "Co-occur" first → bottom
    outcome = factor(outcome, levels = c("No co-occur","Co-occur"))
  )

# draw
plot <- ggplot(plot_df, aes(x = reg_AI_type, y = percent, fill = outcome,  order = outcome  )) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = percent_format(scale = 1)) +   # show “%” scale
  facet_wrap(~ gene_cognate) +
  labs(
    x    = "reg-AI type",
    y    = "Percent of exonic SNVs",
    fill = ""
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text  = element_text(size = 10)
  )

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, plot, width = 150, height = 150, units = "mm") 

# Plot
df_plot <- final_summary %>%
    filter(gene_cognate %in% c("cancer_genes","random")) %>%
    group_by(gene_cognate) %>%
    rowwise() %>%
    mutate(
      ci = list(prop.test(
        x = n_cooccur,
        n = n_exonic,
        conf.level = 0.95
      )$conf.int * 100)   # scale to percentage
    ) %>%
    unnest_wider(ci, names_sep = "_") %>%
    rename(lower_ci = ci_1, upper_ci = ci_2)

  # Step 1: Calculate the differences, adding gene_cancer to keep track of categories
difference_df <- df_plot %>%
  group_by(reg_AI_type) %>%
  summarize(
    diff = abs(diff(pct_cooccur)),  # Absolute difference between the two categories
    start_prop = min(pct_cooccur),  # Start of the arrow at the lower bar
    end_prop = max(pct_cooccur),    # End of the arrow at the upper bar
    gene_cognate = gene_cognate[which.min(pct_cooccur)]  # Identify which group is at start (lowest value)
  ) %>%
  ungroup()

# 1) Re‐level & re‐label the reg_AI_type factor in both data frames:
new_levels <- c("Neg reg-AI", "No reg-AI", "Pos reg-AI")
new_labels <- c("Negative",   "No",      "Positive")

df_plot$reg_AI_type <- factor(df_plot$reg_AI_type,
                              levels = new_levels,
                              labels = new_labels)

difference_df$reg_AI_type <- factor(difference_df$reg_AI_type,
                                    levels = new_levels,
                                    labels = new_labels)

# Step 2: Update the plot to add error bars
plot <- ggplot(df_plot,
               aes(x    = reg_AI_type,
                   y    = pct_cooccur,
                   fill = gene_cognate)) +
  # bars + error bars
  geom_bar(stat     = "identity",
           position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin     = lower_ci,
                    ymax     = upper_ci),
                position = position_dodge(width = 0.9),
                width    = 0.25) +
  # arrows + Δ labels (as before)
  geom_segment(data = difference_df,
               aes(x     = reg_AI_type,
                   xend  = reg_AI_type,
                   y     = start_prop,
                   yend  = end_prop),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "black",
               size  = 0.5) +
  geom_text(data = difference_df,
            aes(x     = reg_AI_type,
                y     = end_prop + 5,
                label = paste0("Δ", round(diff, 1), "%")),
            size   = 3,
            vjust  = 0,
            hjust  = 0.5) +
  # labels
  labs(x    = "reg-AI",
       y    = "Proportion (%) of SNVs",
       fill = "") +
  # 3) Custom fill labels
  scale_fill_brewer(palette = "Set3",
                    labels  = c("Cancer genes", "Passenger genes")) +
  # theme + y‐axis limits
  theme_classic() +
  theme(legend.position = "right",
        axis.text.x    = element_text(hjust = 0.5)) +
  ylim(0, 100)

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/test.png"
ggsave(final_figure_path, plot, width = 150, height = 100, units = "mm") 

res_list = list(plot_df = df_plot, annot_df = difference_df)
# Save
# write.table(res_list, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig3/panel_A.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(res_list, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig3/panel_A.RData")

# 5. now extract the two examples for oncogene & Pos reg-AI that do co-occur
# oncogene_pos_examples <- exonic_with_flag %>%
#   filter(
#     gene_cognate == "oncogene",
#     reg_AI_type == "Neg reg-AI",
#     has_noncoding
#   ) #%>%
#   # slice_head(n = 5)

# # 6. view them
# oncogene_pos_examples

### Cancer genes OR

# 1. Subset to cancer_genes & the two reg-AI types
cs <- subset(final_summary,
             gene_cognate == "cancer_genes" &
             reg_AI_type %in% c("Pos reg-AI", "No reg-AI"))

# 2. Compute co-occur and no-co-occur counts
counts <- transform(cs,no_cooccur = n_exonic - n_cooccur)

# 3. Build contingency matrix:
mat <- with(counts,
            matrix(c(n_cooccur, no_cooccur),
                   nrow = 2, byrow = FALSE,
                   dimnames = list(
                     reg_AI_type = reg_AI_type,
                     outcome     = c("Cooccur", "NoCooccur")
                   )
            ))
print(mat)

mat <- mat[c(2,1),c(2,1)]
# 4. Fisher’s exact test (two-sided)
f2 <- fisher.test(mat, alternative = "two.sided")
print(f2)

# 5. Fisher’s exact test (one-sided enrichment in Pos reg-AI)
f1 <- fisher.test(mat, alternative = "greater")
print(f1)

### Cancer genes vs passenger (pos reg-AI) OR

# 1. Subset to cancer_genes & the two reg-AI types
cs <- subset(final_summary,
             gene_cognate %in% c("cancer_genes","random") &
             reg_AI_type %in% c("Pos reg-AI"))

# 2. Compute co-occur and no-co-occur counts
counts <- transform(cs,no_cooccur = n_exonic - n_cooccur)

# 3. Build contingency matrix:
mat <- with(counts,
            matrix(c(n_cooccur, no_cooccur),
                   nrow = 2, byrow = FALSE,
                   dimnames = list(
                     gene_cognate = gene_cognate,
                     outcome     = c("Cooccur", "NoCooccur")
                   )
            ))
print(mat)

mat <- mat[c(2,1),c(2,1)]
# 4. Fisher’s exact test (two-sided)
f2 <- fisher.test(mat, alternative = "two.sided")
print(f2)

# 5. Fisher’s exact test (one-sided enrichment in Pos reg-AI)
f1 <- fisher.test(mat, alternative = "greater")
print(f1)

###################################################################################################
################################# SELECTION IN MRNA AND REG-AI ####################################
###################################################################################################

###################################################
################### FIGURE 3A #####################
###################################################

# Schematic (REMOVED, no schematic at the end)

#####################################################################
################### FIGURE 3B and SUPP. FIG S4A #####################
#####################################################################

TCGA_AI_table_final_filt <- TCGA_AI_table_final %>%
  mutate(SNV_type = as.character(SNV_varity)) %>%
  mutate(SNV_type = ifelse(SNV_type == "stopgain",
                           paste("stopgain:", NMD_status),
                           SNV_type)) %>%
  mutate(SNV_type = factor(SNV_type, levels = c("effectively_syn", 
                                                "nonsynonymous SNV", 
                                                "stopgain: NMD-triggering", 
                                                "stopgain: NMD-evading"))) %>%
  filter(!gene_cognate %in% c("both_cognate","both","essential"))

table(TCGA_AI_table_final_filt$SNV_type)

# 1. pivot longer so that both AI vars live in one column
df_long <- TCGA_AI_table_final_filt %>%
  pivot_longer(
    cols      = c(reg_AI_type, mRNA_AI_type),
    names_to  = "AI_dimension",
    values_to = "AI_level"
  ) %>%
  # 2. clean up names for plotting
  mutate(
    AI_dimension = recode(AI_dimension,
                          reg_AI_type   = "reg-AI",
                          mRNA_AI_type  = "mRNA-AI"),
    AI_level     = recode(AI_level,
                          `Pos reg-AI` = "Pos",
                          `No reg-AI`  = "No",
                          `Neg reg-AI` = "Neg",
                          # if your mRNA levels are named differently, map them here:
                          `Pos mRNA-AI`   = "Pos",
                          `No mRNA-AI` = "No",
                          `Neg mRNA-AI`    = "Neg"
    )
  ) %>%
  # 3. compute counts & proportions within each AI_dimension × AI_level × gene_cognate
  group_by(gene_cognate, SNV_type, AI_dimension, AI_level, gene) %>%
  summarise(Freq = n(), .groups = "drop") %>%
  group_by(gene_cognate, AI_dimension, AI_level) %>%
  mutate(Prop = Freq / sum(Freq))

generate_plot_swapAI <- function(df, output_file) {
  # ensure AI_dimension is a true factor in the right order:
  df <- df %>%
    mutate(
      gene_cognate  = factor(gene_cognate,
                             levels = c("random", "oncogene", "og_cognate",
                                        "tsg", "tsg_cognate", "TP53")),
      AI_dimension = factor(AI_dimension, levels = c("reg-AI", "mRNA-AI")),
      AI_level     = factor(AI_level,     levels = c("Pos", "No", "Neg"))
    )

  p <- ggplot(df,
              aes(
                x    = AI_dimension,
                y    = Prop,
                fill = factor(SNV_type,
                              levels = c("nonsynonymous SNV",
                                         "stopgain: NMD-triggering",
                                         "stopgain: NMD-evading",
                                         "effectively_syn"))
              )) +
    geom_col() +
    geom_text(
      aes(label = ifelse(Prop < 0.05,
                         "", 
                         scales::percent(Prop, accuracy = 1))
      ),
      position = position_stack(vjust = 0.5),
      size = 3,
      na.rm = TRUE
    ) +
    scale_fill_manual(
      limits = c("effectively_syn",
                 "stopgain: NMD-triggering",
                 "stopgain: NMD-evading",
                 "nonsynonymous SNV"),
      labels = c(
        "nonsynonymous SNV"        = "High aa-impact missense",
        "effectively_syn"          = "ES",
        "stopgain: NMD-triggering" = "NMD-Tr nonsense",
        "stopgain: NMD-evading"    = "NMD-Ev nonsense"
      ),
      values = c("#23cee9", "#bc502c", "#f0edd8", "#ffdd00")
    ) +
    # rename the two x-ticks if you like:
    scale_x_discrete(labels = c("reg-AI" = "reg-AI",
                                "mRNA-AI" = "mRNA-AI")) +

    # now facets: rows = your Pos/No/Neg, cols = each gene
    facet_grid(
      rows = vars(AI_level),
      cols = vars(gene_cognate),
      scales = "free_y",
      switch = "y",
      labeller = labeller(
        AI_level = c(Pos = "Pos", No = "No", Neg = "Neg"),
        gene_cognate     = c(
          random   = "Passenger",
          oncogene = "OG",
          og_cognate = "Cog. OG",
          tsg      = "TSG",
          tsg_cognate = "Cog. TSG",
          TP53     = "TP53"
        )
      )
    ) +
    theme_classic() +
    theme(
      strip.placement = "outside",   # move the row labels to the very left
      strip.text      = element_text(size = 9),
      axis.title.x    = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.text.x     = element_text(size = 9),
      legend.position = "top",
      legend.title    = element_blank()
    ) +
    labs(y = NULL)

  ggsave(output_file, p, width = 200, height = 120, units = "mm")
}

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
generate_plot_swapAI(df_long, final_figure_path)

# Save
write.table(df_long, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig3/Fig3B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(df_long, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig3/Fig3B.RData")

write.table(df_long, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig4/panel_A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(df_long, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig4/panel_A.RData")

###################################################
################### FIGURE 3C #####################
###################################################

df_tmp <- TCGA_AI_table_final %>%
  filter(SNV_varity != "stopgain")

# Sei_threshold <- get_top_pct(TCGA_AI_table_final$Sei_MaxDelta, 5)
Sei_threshold <- 0.3178995 # pre-calculated top 10%
# percentile 25% of TERT pre-calculated
# Puffin_threshold <- quantile(pos_control$Puffin_MaxDelta,seq(0,1,0.05))["25%"]
Puffin_threshold <- 0.06665037

# 1️⃣ Define for each AI‐type the levels and pairwise comparisons
ai_configs <- list(
  mRNA = list(
    col    = "mRNA_AI_type",
    levels = c("No mRNA-AI", "Pos mRNA-AI", "Neg mRNA-AI"),
    comps  = list(
      Pos_vs_No = c("No mRNA-AI", "Pos mRNA-AI"),
      Neg_vs_No = c("No mRNA-AI", "Neg mRNA-AI")
    )
  ),
  reg = list(
    col    = "reg_AI_type",
    levels = c("No reg-AI", "Pos reg-AI", "Neg reg-AI"),
    comps  = list(
      Pos_vs_No = c("No reg-AI", "Pos reg-AI"),
      Neg_vs_No = c("No reg-AI", "Neg reg-AI")
    )
  )
)

# 2️⃣ Generalized function to run all Fisher tests for one AI‐type
run_OR_tests <- function(data, ai_type, cancer_type, exclude_splice = FALSE) {
  cfg   <- ai_configs[[ai_type]]
  ai_col <- sym(cfg$col)

  # optionally drop splice‐disruptors
  if (exclude_splice) {
    data <- data %>%
      filter(
        (DS_max_abs <= 0.1 | Pangolin_abs_max <= 0.1) |
          is.na(DS_max_abs) | is.na(Pangolin_abs_max)
      ) %>%
      filter(Puffin_MaxDelta <= Puffin_threshold | is.na(Puffin_MaxDelta)) %>%
      filter(Sei_MaxDelta <= Sei_threshold | is.na(Sei_MaxDelta))
  }
  if (cancer_type == "pancancer") {
    data <- data
  } else {
    data <- data %>%
    filter(cancer == cancer_type)
  }

  genes <- unique(as.character(data$gene_cognate))

  map_dfr(genes, function(g) {
    map2_dfr(cfg$comps, names(cfg$comps), function(levels, comp_name) {
      # 1) Subset and re‐factor the AI column
      sub <- data %>%
        filter(gene_cognate == g, !!ai_col %in% levels) %>%
        mutate(!!ai_col := factor(!!ai_col, levels = levels))
      
      # 2) Build a 2×2 contingency table:
      tbl2 <- table(sub$SNV_varity == "nonsynonymous SNV", sub[[cfg$col]])
      
      # 3) Check that the table is at least 2×2:
      if (nrow(tbl2) < 2 || ncol(tbl2) < 2) {
        # If we don't have 2 rows AND 2 columns, return NAs for this combination
        return(tibble(
          gene_cognate    = g,
          comparison      = comp_name,
          OR              = NA_real_,
          CI_lower        = NA_real_,
          CI_upper        = NA_real_,
          p_two_sided     = NA_real_,
          p_upper_tail    = NA_real_,
          p_lower_tail    = NA_real_,
          AI_type         = ai_type,
          excluded_splice = exclude_splice
        ))
      }
      
      # 4) Otherwise, run all three Fisher tests
      ft_ts <- fisher.test(tbl2, alternative = "two.sided", conf.int    = TRUE)
      ft_gt <- fisher.test(tbl2, alternative = "greater")
      ft_lt <- fisher.test(tbl2, alternative = "less")
      
      # 5) Return a tibble with all the requested fields
      tibble(
        gene_cognate    = g,
        comparison      = comp_name,
        OR              = unname(ft_ts$estimate),
        CI_lower        = ft_ts$conf.int[1],
        CI_upper        = ft_ts$conf.int[2],
        p_two_sided     = ft_ts$p.value,
        p_upper_tail    = ft_gt$p.value,
        p_lower_tail    = ft_lt$p.value,
        AI_type         = ai_type,
        excluded_splice = exclude_splice
      )
    })
  })
}

# 3️⃣ Run it for both AI‐types and both settings
all_AI_res <- bind_rows(
  run_OR_tests(df_tmp, "mRNA", cancer_type = "pancancer", exclude_splice = FALSE),
  run_OR_tests(df_tmp, "reg",  cancer_type = "pancancer", exclude_splice = FALSE)
)

all_AI_res_no_splice <- bind_rows(
  run_OR_tests(df_tmp, "mRNA", cancer_type = "pancancer", exclude_splice = TRUE),
  run_OR_tests(df_tmp, "reg",  cancer_type = "pancancer", exclude_splice = TRUE)
)

all_AI_final <- rbind(all_AI_res,all_AI_res_no_splice)

# 4️⃣ For the Manuscript
all_AI_res %>%
  filter(gene_cognate == "random", !excluded_splice) %>%
  data.frame()

all_AI_res_no_splice %>%
  filter(gene_cognate == "random", excluded_splice) %>%
  data.frame()

# A) Plot

# Prepare data
plot_df <- all_AI_final %>%
  # keep only Pos_vs_No, drop unwanted categories
  filter(
    comparison == "Pos_vs_No",
    !gene_cognate %in% c("both", "both_cognate", "essential")
  ) %>%
  mutate(
    splice_group = ifelse(excluded_splice,
                          "W/o SNVs impacting splicing or transcription",
                          "All SNVs"),
    # significance categories for legend
    sig_cat = case_when(
      p_two_sided < 0.001 ~ "*** p < 0.001",
      p_two_sided < 0.01  ~ "** p < 0.01",
      p_two_sided < 0.05  ~ "* p < 0.05",
      TRUE                ~ NA_character_
    ),
    # recode gene names
    gene_cognate = fct_recode(gene_cognate,
      Passenger  = "random",
      `Noncog. OG`         = "oncogene",
      `Cog. OG`  = "og_cognate",
      `Noncog. TSG`        = "tsg",
      `Cog. TSG` = "tsg_cognate",
      TP53       = "TP53"
    ),
    # enforce the exact order you want
    gene_cognate = fct_relevel(
      gene_cognate,
      "Passenger",
      "Noncog. OG",
      "Cog. OG",
      "Noncog. TSG",
      "Cog. TSG",
      "TP53"
    ),
  sig_cat = factor(
      sig_cat,
      levels = c("* p < 0.05", "** p < 0.01", "*** p < 0.001")
    )
  )

final_plot <- ggplot(plot_df,
                     aes(
                       x     = gene_cognate,
                       y     = log2(OR),
                       fill  = splice_group,
                       group = splice_group  # so dodge knows there are two sub‐points
                     )) +
  # horizontal line at log2(OR)=0
  geom_hline(yintercept = 0,
             linetype    = "dashed",
             color       = "black") +

  # 1) point + CI ribbon, dodged by splice_group
  geom_pointrange(
    aes(
      ymin = log2(CI_lower),
      ymax = log2(CI_upper)
    ),
    position = position_dodge2(
      width    = 0.6,
      preserve = "single"
    ),
    shape  = 21,    # filled circle
    color  = "black",
    size   = 0.8
  ) +

  # 2) put the “star / cross / asterisk” shape on top of each point when it’s significant
  geom_point(
    aes(shape = sig_cat),
    position = position_dodge2(
      width    = 0.6,
      preserve = "single"
    ),
    size = 3,
    color = "black",
    na.rm = TRUE
  ) +

  # facet by AI_type
  facet_wrap(~ AI_type,
             labeller = labeller(
               AI_type = c(mRNA = "mRNA-AI", reg = "reg-AI")
             )) +

  theme_bw() +
  theme(
    legend.position   = "top",
    legend.box        = "vertical",
    axis.text.x       = element_text(angle = 45, hjust = 1),
    legend.spacing.y  = unit(0.1, "cm"),
    legend.key.height = unit(0.4, "cm"),
    strip.background  = element_rect(fill = "white")
  ) +

  # viridis fill scale for the two splice_group categories
  scale_fill_viridis_d(
    option    = "viridis",
    begin     = 0.3,
    end       = 0.7,
    direction = 1
  ) +

  # Manually assign shapes to each sig_cat level
  scale_shape_manual(
    values = c(
      `* p < 0.05`   = 10,  # open circle‐plus
      `** p < 0.01`  = 4,   # “X”
      `*** p < 0.001`= 8    # “asterisk”
    ),
    na.translate = FALSE    # drop NA so that rows with sig_cat=NA simply draw no point
  ) +

  # axis/legend labels
  labs(
    x     = NULL,
    fill  = NULL,
    shape = NULL,
    y     = "log₂(OR)"
  ) +

  # Force the fill‐legend keys to show black‐bordered circles:
  guides(
    fill = guide_legend(
      override.aes = list(
        shape  = 21,
        colour = "black"
      )
    ),
    shape = guide_legend(
      override.aes = list(
        fill = NA  # leave the star/X/asterisk unfilled
      )
    )
  )

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, final_plot, width = 175, height = 250, units = "mm") 

# Save
write.table(plot_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig3/Fig3C.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(plot_df, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig3/Fig3C.RData")

# B) By cancer type
cancers <- unique(as.character(df_tmp$cancer))

all_AI_final_cancers <- c()

for (cancer_type in cancers) {
  print(cancer_type)
  all_AI_res <- bind_rows(
    run_OR_tests(df_tmp, "mRNA", cancer_type = cancer_type, exclude_splice = FALSE),
    run_OR_tests(df_tmp, "reg",  cancer_type = cancer_type, exclude_splice = FALSE)
  )

  all_AI_res_no_splice <- bind_rows(
    run_OR_tests(df_tmp, "mRNA", cancer_type = cancer_type, exclude_splice = TRUE),
    run_OR_tests(df_tmp, "reg",  cancer_type = cancer_type, exclude_splice = TRUE)
  )

  all_AI_final <- rbind(all_AI_res,all_AI_res_no_splice)
  all_AI_final$cancer_type <- cancer_type
  all_AI_final_cancers <- rbind(all_AI_final_cancers,all_AI_final)

}

# Prepare data
plot_df <- all_AI_final_cancers %>%
  # keep only Pos_vs_No, drop unwanted categories
  filter(
    OR != Inf & OR != -Inf & !is.na(OR) & OR != 0,
    comparison == "Pos_vs_No",
    !gene_cognate %in% c("both", "both_cognate", "essential")
  ) %>%
  mutate(
    splice_group = ifelse(excluded_splice,
                          "No splice-disrupting SNVs",
                          "All SNVs"),
    # significance categories for legend
    sig_cat = case_when(
      p_two_sided < 0.001 ~ "*** p < 0.001",
      p_two_sided < 0.01  ~ "** p < 0.01",
      p_two_sided < 0.05  ~ "* p < 0.05",
      TRUE                ~ NA_character_
    ),
    # recode gene names
    gene_cognate = fct_recode(gene_cognate,
      Passenger  = "random",
      `Noncog. OG`         = "oncogene",
      `Cog. OG`  = "og_cognate",
      `Noncog. TSG`        = "tsg",
      `Cog. TSG` = "tsg_cognate",
      TP53       = "TP53"
    ),
    # enforce the exact order you want
    gene_cognate = fct_relevel(
      gene_cognate,
      "Passenger",
      "Noncog. OG",
      "Cog. OG",
      "Noncog. TSG",
      "Cog. TSG",
      "TP53"
    ),
  sig_cat = factor(
      sig_cat,
      levels = c("* p < 0.05", "** p < 0.01", "*** p < 0.001")
    )
  )

# Get the list of unique gene_cognate levels:
all_genes <- unique(plot_df$gene_cognate)

# Output directory for per‐gene plots
out_dir <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/per_gene/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Loop over each gene_cognate:
for (g in all_genes) {
  
  # 1) Subset to this one gene_cognate
  df_gene <- plot_df %>% filter(gene_cognate == g)
  
  # 2) Compute average OR by cancer_type (averaging over splice_group)
  #    We assume OR is numeric and non‐missing for each combination.
  avg_order <- df_gene %>%
    group_by(cancer_type) %>%
    summarize(avg_OR = mean(OR, na.rm = TRUE)) %>%
    arrange(desc(avg_OR)) %>%
    pull(cancer_type)
  
  # 3) Re‐factor cancer_type in df_gene so that levels follow avg_order
  df_gene <- df_gene %>%
    mutate(
      cancer_type = factor(cancer_type, levels = avg_order)
    )
  
  # 4) Build the plot for this single gene_cognate
  p <- ggplot(df_gene,
              aes(
                x     = cancer_type,
                y     = log2(OR),
                fill  = splice_group,
                group = splice_group
              )) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    
    # 4a) point + CI, dodged by splice_group
    geom_pointrange(
      aes(
        ymin = log2(CI_lower),
        ymax = log2(CI_upper)
      ),
      position = position_dodge2(width = 0.6, preserve = "single"),
      shape  = 21,
      color  = "black",
      size   = 0.8
    ) +
    
    # 4b) significance stars on top
    geom_point(
      aes(shape = sig_cat),
      position = position_dodge2(width = 0.6, preserve = "single"),
      size = 3,
      color = "black",
      na.rm = TRUE
    ) +
    
    # 4c) facet by AI_type
    facet_wrap(~ AI_type,
               labeller = labeller(AI_type = c(
                 mRNA = "mRNA-AI",
                 reg  = "reg-AI"
               )),
               scales = "free_x",      # allow different x‐scales if desired
               nrow   = 1               # one row of facets
    ) +
    
    theme_bw() +
    theme(
      legend.position   = "top",
      legend.box        = "vertical",
      axis.text.x       = element_text(angle = 45, hjust = 1),
      legend.spacing.y  = unit(0.1, "cm"),
      legend.key.height = unit(0.4, "cm"),
      strip.background  = element_rect(fill = "white"),
      plot.title        = element_text(hjust = 0.5)
    ) +
    
    scale_fill_viridis_d(
      option    = "viridis",
      begin     = 0.3,
      end       = 0.7,
      direction = 1
    ) +
    
    scale_shape_manual(
      values       = c(
        `* p < 0.05`   = 10,  # open circle‐plus
        `** p < 0.01`  = 4,   # “X”
        `*** p < 0.001`= 8    # “asterisk”
      ),
      na.translate = FALSE
    ) +
    
    labs(
      title = paste("gene_cognate =", g),
      x     = "Cancer type (ordered by avg OR)",
      y     = "log₂(OR)",
      fill  = NULL,
      shape = NULL
    ) +
    
    guides(
      fill = guide_legend(
        override.aes = list(
          shape  = 21,
          colour = "black"
        )
      ),
      shape = guide_legend(
        override.aes = list(fill = NA)
      )
    )
  
  # 5) Print to console (optional, if you want to see it interactively)
  print(p)
  
  # 6) Save to file (PNG, 300 dpi)
  out_file <- file.path(out_dir, paste0("OR_by_cancer_", g, ".png"))
  ggsave(
    filename = out_file,
    plot     = p,
    width    = 300,         # mm
    height   = 200,         # mm (adjust as needed)
    units    = "mm",
    dpi      = 300
  )
}

# C) Nonsense mutations ORs

# All nonsense - reg-AI
df_tmp <- TCGA_AI_table_final %>%
      mutate(SNV_type = if_else(SNV_varity == "stopgain","nonsense","other")) %>%
      filter(reg_AI_type %in% c("Neg reg-AI","No reg-AI")) %>%
      mutate(reg_AI_type = as.character(reg_AI_type)) %>%
      mutate(SNV_type = factor(SNV_type, levels = c("other","nonsense"))) %>%
      mutate(reg_AI_type = factor(reg_AI_type, levels = c("No reg-AI","Neg reg-AI")))
tbl <- table(df_tmp$SNV_type,df_tmp$reg_AI_type)
fisher.test(tbl)

# All nonsense - mRNA-AI
df_tmp <- TCGA_AI_table_final %>%
      mutate(SNV_type = if_else(SNV_varity == "stopgain","nonsense","other")) %>%
      filter(mRNA_AI_type %in% c("Neg mRNA-AI","No mRNA-AI")) %>%
      mutate(mRNA_AI_type = as.character(mRNA_AI_type)) %>%
      mutate(SNV_type = factor(SNV_type, levels = c("other","nonsense"))) %>%
      mutate(mRNA_AI_type = factor(mRNA_AI_type, levels = c("No mRNA-AI","Neg mRNA-AI")))
tbl <- table(df_tmp$SNV_type,df_tmp$mRNA_AI_type)
fisher.test(tbl)

# NMD-triggering nonsense - reg-AI - for each gene category
df_tmp <- TCGA_AI_table_final %>%
      filter(gene_cognate == "random") %>%
      mutate(SNV_type = if_else(SNV_varity == "stopgain","nonsense","other")) %>%
      filter( (NMD_status == "NMD-triggering" & SNV_varity == "stopgain") |
            SNV_varity != "stopgain") %>%
      filter(reg_AI_type %in% c("Neg reg-AI","No reg-AI")) %>%
      mutate(reg_AI_type = as.character(reg_AI_type)) %>%
      mutate(SNV_type = factor(SNV_type, levels = c("other","nonsense"))) %>%
      mutate(reg_AI_type = factor(reg_AI_type, levels = c("No reg-AI","Neg reg-AI")))
tbl <- table(df_tmp$SNV_type,df_tmp$reg_AI_type)
tbl
fisher.test(tbl)

# Cognate TSGs NMD-evading
df_tmp <- TCGA_AI_table_final %>%
      filter(gene_cognate == "tsg_cognate") %>%
      mutate(SNV_type = if_else(SNV_varity == "stopgain","nonsense","other")) %>%
      filter( (NMD_status == "NMD-evading" & SNV_varity == "stopgain") |
            SNV_varity != "stopgain") %>%
      filter(reg_AI_type %in% c("Pos reg-AI","No reg-AI")) %>%
      mutate(reg_AI_type = as.character(reg_AI_type)) %>%
      mutate(SNV_type = factor(SNV_type, levels = c("other","nonsense"))) %>%
      mutate(reg_AI_type = factor(reg_AI_type, levels = c("No reg-AI","Pos reg-AI")))
tbl <- table(df_tmp$SNV_type,df_tmp$reg_AI_type)
fisher.test(tbl)

#######################################################
################### SUPP. FIG S4B #####################
#######################################################

# ── 1. data ────────────────────────────────────────────────────────────────
df_long <- TCGA_AI_table_final %>%
  filter(gene_cognate %in% c("tsg_cognate", "TP53", "random"),
        SNV_type == "nonsynonymous SNV") %>%
  mutate(
    delta_RNA_regAI = RNA_observed_ratio - RNA_expected_ratio_bb_model1,
    delta_RNA_mRNAI = RNA_observed_ratio - RNA_expected_ratio_null_model,
    group = ifelse(gene_cognate == "random",
                   "Passengers", "cognate TSGs")
  ) %>%
  pivot_longer(
    cols      = c(delta_RNA_regAI, delta_RNA_mRNAI),
    names_to  = "AI_type",
    values_to = "delta_RNA"
  ) %>%
  mutate(
    AI_type = recode(AI_type,
                     delta_RNA_regAI = "reg-AI",
                     delta_RNA_mRNAI = "mRNA-AI")
  )

# colour palette we will use in several places
col_vals <- c("Passengers"   = "grey60",
              "cognate TSGs" = "#1B9E77")

# ── 2. plot ────────────────────────────────────────────────────────────────
plot <- ggplot(
      df_long,
      aes(VARITY_ER, delta_RNA,
          colour = group, shape = group,
          size   = group,  alpha = group)
    ) +

  # points ---------------------------------------------------------------
  geom_point() +

  # coloured trend lines -------------------------------------------------
  geom_smooth(method = "lm", se = FALSE) +

  # thin black outline for the grey Passenger trend ---------------------
  geom_smooth(
    data = df_long %>% filter(group == "Passengers"),
    aes(VARITY_ER, delta_RNA),
    colour = "black", linetype = "dashed",
    size = 0.4, se = FALSE, inherit.aes = FALSE
  ) +

  # Pearson labels: one for each group, coloured -------------------------
  # Passenger correlation label (grey)
  stat_cor(
    data     = df_long %>% filter(group == "Passengers"),
    mapping  = aes(x = VARITY_ER, y = delta_RNA, colour = "black"),
    method   = "pearson",
    label.x.npc = "left",  label.y.npc = 0.95,
    p.accuracy  = 0.001, r.accuracy = 0.01,
    size     = 3.8, inherit.aes = FALSE
  ) +

  # Cognate-TSG correlation label (green)
  stat_cor(
    data     = df_long %>% filter(group == "cognate TSGs"),
    mapping  = aes(x = VARITY_ER, y = delta_RNA, colour = group),
    method   = "pearson",
    label.x.npc = "left",  label.y.npc = 0.85,
    p.accuracy  = 0.001, r.accuracy = 0.01,
    size     = 3.8, inherit.aes = FALSE
  ) +

  # scales ----------------------------------------------------------------
  scale_colour_manual(values = col_vals) +
  scale_shape_manual (values = c("Passengers"   = 4,  # cross
                                 "cognate TSGs" = 16)) +
  scale_size_manual  (values = c("Passengers"   = 1.0,
                                 "cognate TSGs" = 2.0)) +
  scale_alpha_manual (values = c("Passengers"   = 0.25,
                                 "cognate TSGs" = 0.6)) +

  # facet, labels, theme --------------------------------------------------
  facet_wrap(~AI_type, nrow = 1, scales = "free_y") +
  labs(x = "VARITY score",
       y = "Allelic imbalance",
       colour = "") +
  theme_classic() +
  theme(legend.position = "top") +
  guides(             # keep only colour legend
    shape = "none",
    size  = "none",
    alpha = "none"
  )

# 3. export ------------------------------------------------------------------
ggsave("/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png",
       plot, width = 150, height = 75, units = "mm")

write.table(df_long, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig4/panel_B.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(df_long, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig4/panel_B.RData")

# For the Manuscript
# df <- input_figureD[input_figureD$group == "cognate TSGs",]
# cor_res <- cor.test(df$VARITY_ER,df$log_odds_total, method = "pearson")
# cor_res$p.value
# cor.test(df_tmp$delta_RNA_regAI,df_tmp$VARITY_ER)
# cor.test(df_tmp$delta_DNA_CNA_AI,df_tmp$VARITY_ER)

#####################################################################
################### FIGURE 3D and SUPP. FIG S4C #####################
#####################################################################

# mRNA-AI and reg-AI AI-dNdES gene-level BB-models for selection

# AI dN/dES test, two gene sets:

# •	(A) 356 genes --> (CGC genes) + (Mutpanning in Passenger but not found in CGC). Removed translocation/etc from CGC.
# o	10% FDR --> 8 genes reg-AI, 18 genes (1 Passenger) mRNA-AI
# o	25% FDR --> 22 reg-AI (1 random), 35 (2 Passenger) mRNA-AI
# o	Inflation (lambda) --> reg-AI (1.41), mRNA-AI (1.62)

# •	(B) 793 genes --> (CGC Genes) + (Mutpanning in Passenger but not found in CGC)
# o	10% FDR --> 4 genes reg-AI, 20 genes (2 Passenger) mRNA-AI
# o	25% FDR --> 22 reg-AI (5 random), 38 (6 Passenger) mRNA-AI
# o	Inflation (lambda) --> reg-AI (1.08), mRNA-AI (1.25)

# I am currently using gene set (A)

# conda activate R_figures
source("/home/gpalou/projects/ASE/Figures/ggplot_themes.R")
# For genome-wide
TCGA_AI_table_final <- readRDS("/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_AI_final_tables/TCGA_AI_table_bb_models_final.RData")
# For cancer-genes selected
TCGA_AI_table_final <- readRDS("/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_AI_final_tables/TCGA_AI_table_bb_models_final_2.RData")
dim(TCGA_AI_table_final)

# A) Read data

# TCGA AI table
TCGA_AI_table_final <- TCGA_AI_table_final %>%
  filter(SNV_varity != "stopgain")
dim(TCGA_AI_table_final)
# CGC genes & Mutpanning_genes

# Select cancer genes
# putative_cancer column genes are: All cancer genes + some random that are known to be cognate from Mutpanning (MutPanningGeneTumorPairs.csv)
TCGA_AI_table_final$putative_cancer <- FALSE
# All CGC genes
TCGA_AI_table_final[TCGA_AI_table_final$gene %in% c("TP53", "tsg", "oncogene", "both"),]$putative_cancer <- TRUE 
# Add genes that are in Mutpanning (regardless of cognate)
Mutpanning_all_genes <- unique(as.character(Mutpanning_genes$Gene.Symbol))
TCGA_AI_table_final$putative_cancer[TCGA_AI_table_final$Gene.Symbol %in% 
          TCGA_AI_table_final$Gene.Symbol[TCGA_AI_table_final$gene == "random" & 
          TCGA_AI_table_final$Gene.Symbol %in% Mutpanning_all_genes]] <- TRUE 
# Remove translocation genes from CGC
# CGC_genes_to_remove <- unique(as.character(CGC_genes_table$Gene.Symbol[CGC_genes_table$Mutation.Types == "T"]))
# df_merged[df_merged$Gene.Symbol %in% c(CGC_genes_to_remove),"putative_cancer"] <- FALSE
# Remove any gene that appear only once in Mutpanning
# Mutpanning_genes_to_remove <- unique(as.character(names(table(cognate_table$Gene.Symbol)[table(cognate_table$Gene.Symbol) <= 1])))
# df_merged[df_merged$Gene.Symbol %in% c(Mutpanning_genes_to_remove),"putative_cancer"] <- FALSE
table(TCGA_AI_table_final$putative_cancer)

# Filters
test_type <- "genome_wide" #cancer_genes // genome_wide
if (test_type == "genome_wide") {
  # # If genome_wide --> Remove low gene exp genes:
  gene_exp_median <- TCGA_AI_table_final %>%
    group_by(Gene.Symbol) %>%
    summarise(gene_exp_median = median(gene_exp))
  genes_remove <- as.character(gene_exp_median[which(gene_exp_median$gene_exp_median < 1),] %>% pull(Gene.Symbol))
  print(dim(TCGA_AI_table_final))
  TCGA_AI_table_final <- TCGA_AI_table_final %>%
    filter(!Gene.Symbol %in% genes_remove)
} else if (test_type == "cancer_genes") {
  filter1 <- TCGA_AI_table_final$putative_cancer
  TCGA_AI_table_final <- TCGA_AI_table_final[which(filter1),]
}
print(dim(TCGA_AI_table_final))

# Test by gene

# List of RNA models
RNA_model_formulas <- list(
  RNA_bb_model1 = as.formula(
    "cbind(RNA_alt, RNA_ref) ~ MUT_CNA_fraction + sqrt(gene_exp) + factor(SNV_varity)"
  ),
  RNA_null_model = as.formula(
    "cbind(RNA_alt, RNA_ref) ~ purity + sqrt(gene_exp) + purity:sqrt(gene_exp) + factor(SNV_varity)"
  )
)

# Combine all models into one list
all_models <- RNA_model_formulas
# all_models <- c(RNA_model_formulas, DNA_model_formulas)
list_genes <- unique(as.character(TCGA_AI_table_final$Gene.Symbol))
list_genes
length(list_genes)

#For the Manuscript
TCGA_AI_table_final %>%
  filter(gene == "oncogene") %>%
  pull(Gene.Symbol) %>% unique() %>% length()

# Initialize result dataframe
AI_dNdES_test_res <- data.frame()

# Helper function: safely get count or return 0
get_count <- function(tbl, type) {
  if (type %in% names(tbl)) return(tbl[type])
  else return(0)
}

remove_subclonal <- "no" #yes // no

# Loop through each gene and each model
for (gene_char in list_genes) {
  print ( (which(list_genes %in% gene_char) / length(list_genes)) * 100 )

  df_gene <- TCGA_AI_table_final %>% 
    filter(Gene.Symbol == gene_char)
  # Remove splicing mutations
  # df_gene <- df_gene %>%
  #   filter(Pangolin_abs_max <= 0.15)
  # Remove subclonal variants ?
  if (remove_subclonal == "yes") {
    df_gene <- df_gene %>%
      filter((DNA_alt/(DNA_alt+DNA_ref))/purity >= 0.2)
  }

  if (nrow(df_gene) < 10) next
  # Count SNV types
  snv_counts <- table(df_gene$SNV_varity)
  # print(snv_counts)
  syn_count <- get_count(snv_counts, "effectively_syn")
  non_syn_count <- get_count(snv_counts, "nonsynonymous SNV")

  # Require at least 3 of each type
  if (syn_count < 3 || non_syn_count < 3) next

  for (model_name in names(all_models)) {
    model_formula <- all_models[[model_name]]
    model_type <- ifelse(grepl("^RNA_", model_name), "RNA", "DNA")

    try({
      model <- VGAM::vglm(
        formula = model_formula,
        family = VGAM::betabinomial,
        na.action = na.exclude,
        data = df_gene
      )

      coef_table <- summary(model)@coef3
      row_id <- "factor(SNV_varity)nonsynonymous SNV"

      if (row_id %in% rownames(coef_table)) {
        est <- coef_table[row_id, "Estimate"]
        z <- coef_table[row_id, "z value"]
        p <- 2 * pnorm(-abs(z))
        p_high <- pnorm(z, lower.tail = FALSE)
        p_low <- pnorm(z, lower.tail = TRUE)

        AI_dNdES_test_res <- rbind(AI_dNdES_test_res, data.frame(
          Gene = gene_char,
          Gene_ID = unique(as.character(df_gene$gene)),
          Model = model_name,
          Model_Type = model_type,
          Estimate = est,
          Z_value = z,
          P_value = p,
          Higher_one_sided_p = p_high,
          Lower_one_sided_p = p_low
        ))
      }
    }, silent = TRUE)
  }
}

# View top hits
head(AI_dNdES_test_res)
dim(AI_dNdES_test_res)
table(AI_dNdES_test_res$Model)

path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/bb_models/AI_dNdES_test"
if (test_type == "genome_wide") {
  # output_path <- paste0(path,"/genome_wide/AI_dNdES_test_res_remove_low_gene_exp.txt")
  output_path <- paste0(path,"/genome_wide/AI_dNdES_test_res.txt")
} else if (test_type == "cancer_genes" & remove_subclonal == "yes") {
  output_path <- paste0(path,"/AI_dNdES_test_res_geneset_A_subclonal_removed.txt")
} else if (test_type == "cancer_genes" & remove_subclonal == "no") {
  output_path <- paste0(path,"/AI_dNdES_test_res_geneset_A.txt")
}
print(output_path)

write.table(AI_dNdES_test_res, file = output_path, 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

AI_dNdES_test_res <- read.table(output_path, header = TRUE, sep = "\t")
# For the Manuscript (number of genes before filterings)
length(unique(AI_dNdES_test_res$Gene))
AI_dNdES_test_res %>%
  filter(Model == "RNA_null_model") %>%
  pull(Gene_ID) %>% table()

AI_dNdES_test_res <- AI_dNdES_test_res %>%
  filter(Z_value <= 15 & Z_value >= -15)

AI_dNdES_test_res <- AI_dNdES_test_res %>%
  group_by(Model) %>%
  mutate(
    FDR_two_sided = p.adjust(P_value, method = "fdr"),
    FDR_higher = p.adjust(Higher_one_sided_p, method = "fdr"),
    FDR_lower = p.adjust(Lower_one_sided_p, method = "fdr")
  ) %>%
  ungroup()

# For the Manuscript
# Number of genes after filterings
length(unique(AI_dNdES_test_res$Gene))
AI_dNdES_test_res %>%
  filter(Model == "RNA_null_model") %>%
  pull(Gene_ID) %>% table()
# gene categories
AI_dNdES_test_res %>%
  filter(Model == "RNA_null_model") %>%
  filter(FDR_higher <= 0.10) %>%
  pull(Gene_ID) %>% table()
AI_dNdES_test_res %>%
  filter(Model == "RNA_bb_model1") %>%
  filter(FDR_higher <= 0.10) %>%
  pull(Gene_ID) %>% table()
# genes
genes <- AI_dNdES_test_res %>%
  filter(Model == "RNA_null_model") %>%
  filter(FDR_higher <= 0.25) %>%
  # filter(Gene_ID == "oncogene") %>%
  pull(Gene)
genes

# Fisher test of enrichment of cancer genes over random genes (genome_wide scan only)

df2 <- AI_dNdES_test_res %>%
  filter(Model == "RNA_bb_model1") %>%
  mutate(group = case_when(
    Gene_ID %in% c("TP53", "tsg", "oncogene") ~ "driver",
    Gene_ID == "random"                      ~ "random",
    TRUE                                     ~ NA_character_
  )) %>%
  mutate(group = factor(group, levels = c("random", "driver"))) %>%
  filter(!is.na(group))

ct <- table(
  Group        = df2$group,
  Significant  = df2$FDR_higher <= 0.1
)
print(ct)

ft <- fisher.test(ct)
ft

ft$estimate     # the odds ratio
ft$p.value      # the p-value

# B) ## VOLCANO --> AI-dNdES test

# FDR threshold to use
FDR_10 <- 0.10

# Column names
p_col   <- "Higher_one_sided_p"
fdr_col <- "FDR_higher"

# Prepare the data
df <- AI_dNdES_test_res %>%
  # filter(Model == "RNA_bb_model1") %>%
  # 1) rename models for the legend
  mutate(
    Model = recode(Model,
                   RNA_bb_model1    = "reg-AI",
                   RNA_null_model   = "mRNA-AI"),
    # 2) force TP53 into the TSG group
    Gene_ID = case_when(
      Gene == "TP53"  ~ "tsg",
      TRUE             ~ Gene_ID
    ),
    # 3) relabel Gene_ID categories
    Gene_ID = recode(Gene_ID,
                          oncogene   = "OG",
                          tsg         = "TSG",
                          both        = "Both",
                          random      = "Mutpanning")
  ) %>%
  arrange(.data[[fdr_col]]) %>%
  # 4) flag only those under 10% FDR
  mutate(
    signif10 = .data[[fdr_col]] < FDR_10,
    to_label = ifelse(signif10, Gene, NA_character_)
  )

# Compute per‐model 10% cutoff in –log10 space
thresholds <- df %>%
  group_by(Model) %>%
  summarise(
    max_p10 = max(.data[[p_col]][.data[[fdr_col]] < FDR_10], na.rm = TRUE)
  ) %>%
  transmute(Model, y10 = -log10(max_p10 + 1e-10))

# Build the combined volcano
plot <- ggplot(df,
            aes(
              x = Estimate,
              y = -log10(.data[[p_col]] + 1e-10),
              colour = Model,
              shape  = Gene_ID,
              alpha   = signif10
            )) +
  geom_point(size = 3) +
  scale_alpha_manual(
    values = c(`TRUE` = 1, `FALSE` = 0.5),  # non-significant (FALSE) get 25% opacity
    guide  = FALSE
  ) +
  # single FDR=10% line per model
  geom_hline(
    data = thresholds,
    aes(yintercept = y10, colour = Model),
    linetype = "dashed"
  ) +
  # label only FDR<10% genes
  geom_text_repel(
    aes(label = to_label),
    max.overlaps = 1000,
    show.legend  = FALSE
  ) +
  geom_vline(xintercept = 0, colour = "black") +
  # pastel colors for the two models
  scale_colour_manual(
    values = c("reg-AI" = "#A191B6", "mRNA-AI" = "#b35a4f")
  ) +
  # custom shapes
  scale_shape_manual(
    values = c("OG"         = 16,   # circle
               "TSG"        = 15,   # square
               "Mutpanning" = 8,    # asterisk
               "Both"       = 17)   # triangle
  ) +
  # axis labels & limits
  labs(
    x = "AI-dNdES effect size",
    y = expression(log[10](p-value))
  ) +
  coord_cartesian(xlim = c(-0.5, 1.3)) +
  # clean up legends & theme
  theme_light() +
  theme(
    plot.title     = element_blank(),
    legend.title   = element_blank(),
    legend.position= "top",
    axis.title     = element_text(size = 14),
    axis.text      = element_text(size = 12),
    legend.text    = element_text(size = 12)
  )

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, plot, width = 150, height = 150, units = "mm") 

# Save
write.table(df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig3/Fig3D.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(df, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig3/Fig3D.RData")

write.table(df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig4/panel_C.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(df, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig4/panel_C.RData")

# For the Manuscript

gene_char <- "IDH1"
AI_dNdES_test_res %>%
  filter(Model == "RNA_null_model") %>%
  filter(FDR_higher < 0.10) %>%
  filter(Gene != "TP53") %>% 
  arrange(Estimate) %>% data.frame()

df <- TCGA_AI_table_final %>%
  filter(Gene.Symbol == gene_char) %>%
  filter(reg_AI_type != "Neg reg-AI") %>%
  mutate(SNV_type = as.character(SNV_type)) %>%
  mutate(reg_AI_type = as.character(reg_AI_type))
table(df$SNV_type,df$reg_AI_type)
prop.table(table(df$SNV_type,df$reg_AI_type))

######################################################
################# SUPP. TABLE. S2 ####################
######################################################

# Genome-wide Passenger genes at FDR < 25%

SuppTableS2 <- AI_dNdES_test_res %>%
  # 1) rename Gene_ID → “gene type”
  rename(`gene type` = Gene_ID) %>%
  # 2) recode Model and gene type values
  mutate(
    Model = recode(Model,
                   RNA_bb_model1 = "reg-AI",
                   RNA_null_model = "mRNA-AI"),
    `gene type` = recode(`gene type`,
                         random = "Passenger",
                         tsg    = "TSG",
                         both   = "both TSG/OG")
  ) %>%
  # 3) drop Model_Type
  select(-Model_Type) %>%
  filter(FDR_higher < 0.10) %>%
  filter(`gene type` == "Passenger") %>%
  arrange((FDR_higher))
# check
head(SuppTableS2)

write.table(SuppTableS2, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppTables/SuppTableS2.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

######################################################
################# SUPP. TABLE. S3 ####################
######################################################

# Significant cancer genes at FDR < 25%

SuppTableS3 <- AI_dNdES_test_res %>%
  # 1) rename Gene_ID → “gene type”
  rename(`gene type` = Gene_ID) %>%
  # 2) recode Model and gene type values
  mutate(
    Model = recode(Model,
                   RNA_bb_model1 = "reg-AI",
                   RNA_null_model = "mRNA-AI"),
    `gene type` = recode(`gene type`,
                         random = "Passenger",
                         tsg    = "TSG",
                         "TP53" = "TSG",
                         both   = "both TSG/OG")
  ) %>%
  # 3) drop Model_Type
  select(-Model_Type) %>%
  filter(FDR_higher < 0.25) %>%
  arrange((FDR_higher))
# check
head(SuppTableS3)

write.table(SuppTableS3, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppTables/SuppTableS3.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#######################################################
################### SUPP. FIG S4D #####################
#######################################################

# A) QQ-plot from Volcano plots

gg_qqplot <- function(ps, lambda_value, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # annotate("text", x = 0.5, y = 6, label = paste("λ ==", lambda_value), parse = TRUE) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

plot_df <- AI_dNdES_test_res %>%
        filter(Model == "RNA_null_model")

# Plot
chisq <- qchisq(1-plot_df$Higher_one_sided_p, 1)
lambda <- round(median(chisq,na.rm=TRUE)/qchisq(0.5,1),2)
lambda

plot <- gg_qqplot(plot_df$Higher_one_sided_p, lambda_value = lambda) +
  ggplot_theme_bw() +
  labs(title = paste0("mRNA-AI")) +
  theme(
    axis.ticks = element_line(size = 0.5),
    panel.grid = element_blank()
    # panel.grid = element_line(size = 0.5, color = "grey80")
  )

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
# ggsave(final_figure_path, plot, width = 150, height = 150, units = "mm") 

write.table(AI_dNdES_test_res, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig4/panel_D.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(AI_dNdES_test_res, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig4/panel_D.RData")

######################################################
################# SUPP. FIG S4E ######################
######################################################

# B) ## VOLCANO --> AI-dNdES test

# FDR threshold to use
FDR_10 <- 0.10

# Column names
p_col   <- "Higher_one_sided_p"
fdr_col <- "FDR_higher"

# Prepare the data
df <- AI_dNdES_test_res %>%
  # filter(Model == "RNA_bb_model1") %>%
  # 1) rename models for the legend
  mutate(
    Model = recode(Model,
                   RNA_bb_model1    = "reg-AI",
                   RNA_null_model   = "mRNA-AI"),
    # 2) force TP53 into the TSG group
    Gene_ID = case_when(
      Gene == "TP53"  ~ "tsg",
      TRUE             ~ Gene_ID
    ),
    # 3) relabel Gene_ID categories
    Gene_ID = recode(Gene_ID,
                          oncogene   = "OG",
                          tsg         = "TSG",
                          both        = "Both",
                          random      = "Mutpanning")
  ) %>%
  arrange(.data[[fdr_col]]) %>%
  # 4) flag only those under 10% FDR
  mutate(
    signif10 = .data[[fdr_col]] < FDR_10,
    to_label = ifelse(signif10, Gene, NA_character_)
  )

# Compute per‐model 10% cutoff in –log10 space
thresholds <- df %>%
  group_by(Model) %>%
  summarise(
    max_p10 = max(.data[[p_col]][.data[[fdr_col]] < FDR_10], na.rm = TRUE)
  ) %>%
  transmute(Model, y10 = -log10(max_p10 + 1e-10))

# Build the combined volcano
plot <- ggplot(df,
            aes(
              x = Estimate,
              y = -log10(.data[[p_col]] + 1e-10),
              colour = Model,
              shape  = Gene_ID,
              alpha   = signif10
            )) +
  geom_point(size = 3) +
  scale_alpha_manual(
    values = c(`TRUE` = 1, `FALSE` = 0.5),  # non-significant (FALSE) get 25% opacity
    guide  = FALSE
  ) +
  # single FDR=10% line per model
  geom_hline(
    data = thresholds,
    aes(yintercept = y10, colour = Model),
    linetype = "dashed"
  ) +
  # label only FDR<10% genes
  geom_text_repel(
    aes(label = to_label),
    max.overlaps = 1000,
    show.legend  = FALSE
  ) +
  geom_vline(xintercept = 0, colour = "black") +
  # pastel colors for the two models
  scale_colour_manual(
    values = c("reg-AI" = "#A191B6", "mRNA-AI" = "#b35a4f")
  ) +
  # custom shapes
  scale_shape_manual(
    values = c("OG"         = 16,   # circle
               "TSG"        = 15,   # square
               "Mutpanning" = 8,    # asterisk
               "Both"       = 17)   # triangle
  ) +
  # axis labels & limits
  labs(
    x = "AI-dNdES effect size",
    y = expression(log[10](p-value))
  ) +
  coord_cartesian(xlim = c(-0.5, 1.3)) +
  # clean up legends & theme
  theme_light() +
  theme(
    plot.title     = element_blank(),
    legend.title   = element_blank(),
    legend.position= "top",
    axis.title     = element_text(size = 14),
    axis.text      = element_text(size = 12),
    legend.text    = element_text(size = 12)
  )

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, plot, width = 150, height = 150, units = "mm") 

# AI-dNdEs test after removing subclonal variants

write.table(df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig4/panel_E.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(df, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig4/panel_E.RData")

##################### ON-GOING ##################

# A) ORs fisher test By cancer type

cancers <- unique(as.character(TCGA_AI_table_final$cancer))
# Helper function: safely get count or return 0
get_count <- function(tbl, type) {
  if (type %in% names(tbl)) return(tbl[type])
  else return(0)
}

AI_dNdES_test_res_all <- c()

for (cancer_type in c("pancancer",cancers)) {
  print(cancer_type)
  # Initialize result dataframe
  AI_dNdES_test_res <- data.frame()
  if (cancer_type == "pancancer") {
    df_cancer <- TCGA_AI_table_final
  } else {
    df_cancer <- TCGA_AI_table_final %>%
          filter(cancer == cancer_type)
  }
  # Loop through each gene and each model
  for (gene_char in list_genes) {
    # print ( (which(list_genes %in% gene_char) / length(list_genes)) * 100 )

    df_gene <- df_cancer %>% 
      filter(Gene.Symbol == gene_char)
    # Remove splicing mutations
    # df_gene <- df_gene %>%
    #   filter(Pangolin_abs_max <= 0.15)

    if (nrow(df_gene) < 10) next
    # Count SNV types
    snv_counts <- table(df_gene$SNV_varity)
    # print(snv_counts)
    syn_count <- get_count(snv_counts, "effectively_syn")
    non_syn_count <- get_count(snv_counts, "nonsynonymous SNV")

    # Require at least 3 of each type
    if (syn_count < 3 || non_syn_count < 3) next

    for (model_name in names(all_models)) {
      model_formula <- all_models[[model_name]]
      model_type <- ifelse(grepl("^RNA_", model_name), "RNA", "DNA")

      try({
        model <- VGAM::vglm(
          formula = model_formula,
          family = VGAM::betabinomial,
          na.action = na.exclude,
          data = df_gene
        )

        coef_table <- summary(model)@coef3
        row_id <- "factor(SNV_varity)nonsynonymous SNV"

        if (row_id %in% rownames(coef_table)) {
          est <- coef_table[row_id, "Estimate"]
          z <- coef_table[row_id, "z value"]
          p <- 2 * pnorm(-abs(z))
          p_high <- pnorm(z, lower.tail = FALSE)
          p_low <- pnorm(z, lower.tail = TRUE)

          AI_dNdES_test_res <- rbind(AI_dNdES_test_res, data.frame(
            cancer_type = cancer_type,
            Gene = gene_char,
            Gene_ID = unique(df_gene$gene),
            Model = model_name,
            Model_Type = model_type,
            Estimate = est,
            Z_value = z,
            P_value = p,
            Higher_one_sided_p = p_high,
            Lower_one_sided_p = p_low
          ))
        }
      }, silent = TRUE)
    }
  }

  AI_dNdES_test_res_all <- rbind(AI_dNdES_test_res_all,AI_dNdES_test_res)
}

# View top hits
head(AI_dNdES_test_res_all)
dim(AI_dNdES_test_res_all)
table(AI_dNdES_test_res_all$Model)

output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/bb_models/AI_dNdES_test/_AI_dNdES_test_res_by_cancer.txt")
write.table(AI_dNdES_test_res_all, file = output_path, 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

AI_dNdES_test_res_all <- read.table(output_path, header = TRUE, sep = "\t")
AI_dNdES_test_res_all <- AI_dNdES_test_res_all %>%
  filter(Z_value <= 15 & Z_value >= -15)

AI_dNdES_test_res_all <- AI_dNdES_test_res_all %>%
  group_by(Model,cancer_type) %>%
  mutate(
    FDR_two_sided = p.adjust(P_value, method = "fdr"),
    FDR_higher = p.adjust(Higher_one_sided_p, method = "fdr"),
    FDR_lower = p.adjust(Lower_one_sided_p, method = "fdr")
  ) %>%
  ungroup()

# B) ## VOLCANO --> AI-dNdES test

# FDR threshold to use
FDR_10 <- 0.10

# Column names
p_col   <- "Higher_one_sided_p"
fdr_col <- "FDR_higher"

# Prepare the data
df <- AI_dNdES_test_res_all %>%
  # 1) rename models for the legend
  mutate(
    Model = recode(Model,
                   RNA_bb_model1    = "reg-AI",
                   RNA_null_model   = "mRNA-AI"),
    # 2) force TP53 into the TSG group
    Gene_ID = case_when(
      Gene == "TP53"  ~ "tsg",
      TRUE             ~ Gene_ID
    ),
    # 3) relabel Gene_ID categories
    Gene_ID = recode(Gene_ID,
                          oncogene   = "OG",
                          tsg         = "TSG",
                          both        = "Both",
                          random      = "Mutpanning")
  ) %>%
  arrange(.data[[fdr_col]]) %>%
  # 4) flag only those under 10% FDR
  mutate(
    signif10 = .data[[fdr_col]] < FDR_10,
    to_label = ifelse(signif10, Gene, NA_character_)
  )

# Compute per‐model 10% cutoff in –log10 space
thresholds <- df %>%
  group_by(Model) %>%
  summarise(
    max_p10 = max(.data[[p_col]][.data[[fdr_col]] < FDR_10], na.rm = TRUE)
  ) %>%
  transmute(Model, y10 = -log10(max_p10 + 1e-10))

# Build the combined volcano
final_plot <- ggplot(df,
            aes(
              x = Estimate,
              y = -log10(.data[[p_col]] + 1e-10),
              colour = Model,
              shape  = Gene_ID,
              alpha   = signif10
            )) +
  geom_point(size = 3) +
  facet_wrap(cancer_type ~ ., scales = "free_y") +
  scale_alpha_manual(
    values = c(`TRUE` = 1, `FALSE` = 0.5),  # non-significant (FALSE) get 25% opacity
    guide  = FALSE
  ) +
  # single FDR=10% line per model
  geom_hline(
    data = thresholds,
    aes(yintercept = y10, colour = Model),
    linetype = "dashed"
  ) +
  # label only FDR<10% genes
  geom_text_repel(
    aes(label = to_label),
    max.overlaps = 1000,
    show.legend  = FALSE
  ) +
  geom_vline(xintercept = 0, colour = "black") +
  # pastel colors for the two models
  scale_colour_manual(
    values = c("reg-AI" = "#A191B6", "mRNA-AI" = "#b35a4f")
  ) +
  # custom shapes
  scale_shape_manual(
    values = c("OG"         = 16,   # circle
               "TSG"        = 15,   # square
               "Mutpanning" = 8,    # asterisk
               "Both"       = 17)   # triangle
  ) +
  # axis labels & limits
  labs(
    x = "AI-dNdES effect size",
    y = expression(log[10](p-value))
  ) +
  # coord_cartesian(xlim = c(-0.5, 1.3)) +
  # clean up legends & theme
  theme_light() +
  theme(
    plot.title     = element_blank(),
    legend.title   = element_blank(),
    legend.position= "top",
    axis.title     = element_text(size = 14),
    axis.text      = element_text(size = 12),
    legend.text    = element_text(size = 12)
  )

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, final_plot, width = 475, height = 450, units = "mm") 

########################################

# BY CANCER TYPE: mRNA-AI and reg-AI AI-dNdES gene-level BB-models for selection

# 1) Start from your “volcano” data frame (called df in your code)
#    which already has Model, cancer_type, Gene, FDR_higher, etc.
#
#    We want to:
#    a) Keep only the columns we need: Model, cancer_type, Gene, FDR_higher
#    b) Flag significance at 10% FDR
#    c) Gather into “long” form (if necessary)
#    d) Filter to genes that are ever-significant in any cancer_type/Model
#    e) Compute a plotting value, e.g. -log10(FDR_higher)
#    f) Re‐factor `Gene` so that your rows appear in the order you like
#    g) Re‐factor `cancer_type` so that columns are in a sensible order

FDR_10 <- 0.10

heatmap_df <- AI_dNdES_test_res_all %>%
  # (a) Keep only relevant columns
  select(
    Model,
    cancer_type,
    Gene,
    FDR_higher
  ) %>%
  # (b) Create a logical “is_signif” at 10% FDR
  mutate(
    is_signif = (FDR_higher < FDR_10)
  ) %>%
  # (c) (Optional) If your data were not already one‐row per (Model,cancer_type,Gene),
  #     you’d pivot or group/summarize here. But in your code, it already is one‐row each.
  #
  # (d) Find all genes that are significant in at least one cancer_type (across both Models)
  group_by(Gene) %>%
  mutate(
    ever_signif = any(is_signif, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(ever_signif) %>% 
  # (e) Compute a numeric “fill” value. We’ll do –log10(FDR_higher) so that lower FDR = bigger fill.
  #     Add a tiny epsilon so that log10(0) is avoided if FDR_higher=0.
  mutate(
    log10FDR = -log10(FDR_higher + 1e-10)
  ) %>%
  # (f) Re-factor Gene so that it appears in the order of “most cancer types significant → least”
  #     We can order by how many times it is significant (across all cancer types and both Models),
  #     or by the minimum FDR across all. Here we order by number of significant hits:
  group_by(Gene) %>%
  mutate(
    n_signif_hits = sum(is_signif, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # Then use fct_reorder: genes with more significant hits go at top
  mutate(
    Gene = fct_reorder(Gene, n_signif_hits, .desc = TRUE)
  ) %>%
  # (g) Re-factor cancer_type in some desired order. For example,
  #     you could keep the natural factor levels already on cancer_type,
  #     or reorder by how many signif hits in that cancer_type, etc.
  #     Let’s reorder cancer_type by the total number of significant genes (across all Genes & Models).
  group_by(cancer_type) %>%
  mutate(
    n_sig_in_cancer = sum(is_signif, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    cancer_type = fct_reorder(cancer_type, n_sig_in_cancer, .desc = TRUE)
  ) %>%
  # Keep only the columns needed for plotting
  select(
    Model,
    cancer_type,
    Gene,
    is_signif,
    log10FDR
  )

# # 2. Build the tile‐plot with ggplot2
# heatmap_plot <- ggplot(heatmap_df,
#                        aes(
#                          x    = cancer_type,
#                          y    = Gene,
#                          fill =  ifelse(is_signif, log10FDR, NA_real_)
#                        )) +
#   # A tile for every combination. Even those with is_signif=FALSE will get a light fill (low log10FDR).
#   geom_tile(data = filter(heatmap_df, Model == "RNA_bb_model1"),
#           colour = "black", size = 0.2) +
#   scale_fill_viridis_c(
#     option    = "inferno",
#     direction = -1,
#     na.value  = "grey80", 
#     name      = expression(-log[10](FDR))
#   )  +
#   geom_point(
#     data = filter(heatmap_df, Model == "RNA_null_model" & is_signif),
#     aes(
#       x = cancer_type,
#       y = Gene
#     ),
#     shape    = 4,       # “×” mark
#     color    = "black",
#     size     = 3
#     # position = position_dodge2(width = 0.8, preserve = "single")
#   ) +
#   labs(
#     x = "",
#     y = "",
#     title = "AI-dNdES test by cancer type",
#     subtitle = ""
#   ) +
#   theme_classic(base_size = 13) +
#   theme(
#     axis.text.x        = element_text(angle = 45, hjust = 1),
#     panel.grid.major   = element_blank(),
#     panel.grid.minor   = element_blank(),
#     strip.text         = element_text(face = "bold", size = 14),
#     legend.position    = "right"
#   )

# # 3. (Optional) Save it to disk
# ggsave(
#   "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/genes_by_cancer_heatmap.png",
#   heatmap_plot,
#   width = 200, 
#   height = 300, 
#   units = "mm",
#   dpi = 300
# )

combined_df <- heatmap_df %>%
  # 1) Remove “TCGA-” prefix from the cancer_type factor for cleaner x‐axis labels
  mutate(
    cancer_type = str_remove(as.character(cancer_type), "^TCGA-")
  ) %>%
  # 2) Group by Gene + cancer_type, and figure out whether each Model was significant
  group_by(Gene, cancer_type) %>%
  summarize(
    reg_sig  = any(Model == "RNA_bb_model1" & is_signif, na.rm = TRUE),
    mrna_sig = any(Model == "RNA_null_model" & is_signif, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  # 3) Create a “status” factor with exactly these four categories (plus possible NA)
  mutate(
    status = case_when(
      reg_sig  & !mrna_sig  ~ "reg-AI only",
      !reg_sig & mrna_sig   ~ "mRNA-AI only",
      reg_sig  & mrna_sig   ~ "Both",
      !reg_sig & !mrna_sig  ~ "None",
      TRUE                  ~ NA_character_
    ),
    status = factor(
      status,
      levels = c("reg-AI only", "mRNA-AI only", "Both", "None")
    )
  )

# 1) Compute, for each cancer_type, the number of genes with non‐NA status:
cancer_order <- combined_df %>%
  group_by(cancer_type) %>%
  summarize(
    n_non_na = sum(!is.na(status)),
    .groups = "drop"
  ) %>%
  arrange(desc(n_non_na)) %>%
  pull(cancer_type)

# 1) Compute, for each Gene, how many cancer types are ‘significant’ (status != "None" & !is.na):
gene_order <- combined_df %>%
  group_by(Gene) %>%
  summarize(
    n_cancers_signif = sum(!is.na(status)),
    .groups = "drop"
  ) %>%
  arrange((n_cancers_signif)) %>%
  pull(Gene)

# 2) Re‐factor cancer_type in combined_df so that levels follow that descending count:
combined_df <- combined_df %>%
  mutate(
    cancer_type = factor(cancer_type, levels = cancer_order),
    Gene = factor(Gene, levels = as.character(gene_order))
  )

# Inspect the first few rows
print(combined_df[1:20, ])

# 1) Define a named vector of fill colors for each status
my_colors <- c(
  "reg-AI only"   = "#A191B6",  # a red‐purple
  "mRNA-AI only"  = "#bd5d4a",  # a blue
  "Both"          = "#99D594",  # a green
  "None"          = "#CCCCCC"   # a light gray
)

# 2) Build the tile plot
heatmap_plot <- ggplot(combined_df,
                       aes(
                         x    = cancer_type,
                         y    = Gene,
                         fill = status
                       )) +
  geom_tile(color = "black", size = 0.2) +
  scale_fill_manual(
    values   = my_colors,
    na.value = "white",        # any missing status becomes white
    drop     = FALSE,          # keep all four legend keys, even if some have zero counts
    name     = ""              # blank legend title
  ) +
  labs(
    x     = "",
    y     = "",
    title = "AI-dNdES test (FDR < 10%)",
    subtitle = paste0(
      "Gray = ns  |  White = missing"
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "right"
  )

# 3) (Optional) Save to disk
ggsave(
  "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/genes_by_cancer_status.png",
  heatmap_plot,
  width  = 200,
  height = 300,
  units  = "mm",
  dpi    = 300
)

###################################################
################### FIGURE 3E #####################
###################################################

# 1) Define the genes in the exact order you want
# genes <- c("TP53","IDH1","STK11","BAP1","KEAP1")
# genes <- c("FAT1","CUL3","MET","NF1","PPP2R1A")
# genes <- c("CIC","SMAD3","CDKN2A","RAC1","TGFBR2")
# genes <- c("KRAS","TSC1","EIF1AX","DAXX")
# genes <- c("TP53","CDKN2A","KRAS","KEAP1","TGFBR2","STK11","BAP1")
# genes <- c("IDH1","FAT1","CUL3","MET","NF1","PPP2R1A","BAP1","CIC","SMAD3","RAC1","TSC1","DAXX")

genes <- c("TP53","KRAS","CDKN2A")

# CDKN2A --> reg-AI
# KEAP1, STK11, BAP1, KRAS, TGFBR2 --> CNA-AI

# transcript_id   <- "ENST00000269305" # TP53
# transcript_id <- "ENST00000345146" #IDH1
# transcript_id <- "ENST00000460680" #BAP1
# transcript_id <- "ENST00000171111" #KEAP1
# transcript_id <- "ENST00000326873" #STK11

# 2) Subset & extract codons for all five genes, preserving factor order
df_plot2 <- TCGA_AI_table_final %>%
  filter(Gene.Symbol %in% genes, SNV_varity != "stopgain") %>%
  # make sure Gene.Symbol is a factor in the order you want
  mutate(Gene.Symbol = factor(Gene.Symbol, levels = genes)) %>%
  pivot_longer(
    cols      = c(reg_AI_type, mRNA_AI_type, CNA_AI_type),
    names_to  = "AI_type_category",
    values_to = "AI_type"
  ) %>%
  mutate(
    AI_status = word(AI_type, 1),
    AI_source = case_when(
      AI_type_category == "reg_AI_type"  ~ "reg-AI",
      AI_type_category == "mRNA_AI_type" ~ "mRNA-AI",
      AI_type_category == "CNA_AI_type"  ~ "CNA-AI"
    )
  ) %>% 
  mutate(
    SNV_varity = recode(
      SNV_varity,
      "nonsynonymous SNV" = "High aa-impact missense",
      "effectively_syn"   = "ES"
    )
  )

# Create the 6 bins based on your criteria
df_plot2 <- df_plot2 %>%
  mutate(
    WT_bin = case_when(
      WT_CNA2 == 0 ~ "del",
      WT_CNA2 == 1 & MUT_CNA2 == 1 ~ "neutral", 
      WT_CNA2 > 1 ~ "amp"
    ),
    MUT_bin = case_when(
      MUT_CNA2 == 0 ~ "del",
      MUT_CNA2 == 1 ~ "neutral",
      MUT_CNA2 > 1 ~ "amp"
    ),
    # New MUT_CNA_fraction bins
    CNFmut_bin = case_when(
      MUT_CNA_fraction < 0.4 ~ "<0.4",
      MUT_CNA_fraction >= 0.4 & MUT_CNA_fraction <= 0.6 ~ "0.4-0.6",
      MUT_CNA_fraction > 0.6 ~ ">0.6"
    ),
    CNFmut_bin = factor(CNFmut_bin, levels = c("<0.4", "0.4-0.6", ">0.6")),    # Combine them
    CNA_bin = paste0("WT ",WT_bin,"/ MUT ",MUT_bin),
    # Create ordered factors
    WT_bin = factor(WT_bin, levels = c("del", "neutral", "amp")),
    MUT_bin = factor(MUT_bin, levels = c("del", "neutral", "amp")),
    # CNA_bin = factor(CNA_bin, levels = c(
    #   "WT del / MUT del", "WT del / MUT neut", "WT del / MUT amp",
    #   "WT neut / MUT del", "WT neut / MUT neut", "WT neut / MUT amp", 
    #   "WT amp / MUT del", "WT amp / MUT neut", "WT amp / MUT amp"
    # ))
  ) %>%
  filter(!str_detect(CNA_bin, "NA"))
table(df_plot2$CNA_bin)
table(df_plot2$WT_bin)
table(df_plot2$MUT_bin)
table(df_plot2$CNFmut_bin)

# df_plot2 <- df_plot2 %>%
#   filter(MUT_bin == "neutral")

# 3) Plot with a Gene.Symbol × AI_source grid
final_plot <- ggplot(df_plot2, aes(
    x     = WT_bin,
    y     = RNA_ASE,
    # shape = AI_status
  )) +
  # Single boxplot per bin (no color grouping)
  geom_boxplot() +
  # Colored points overlaid on boxplots
  # Colored and shaped points overlaid on boxplots
  geom_point(aes(color = AI_status, shape = AI_status), 
             size = 1.5, alpha = 0.75, position = position_jitter(width = 0.2)) +
  scale_shape_manual(
    values = c("Neg" = 25, "No" = 4, "Pos" = 19)
  ) +
  labs(x = expression(ACN[wt]), y = "RNA-VAF", shape = "", color = "") +
  theme_bw() +
  # **THIS** is the key: one row per gene, one column per AI source
  facet_grid(
    rows = vars(Gene.Symbol),
    cols = vars(AI_source),
    scales = "free_x",
    space  = "free_x"
  ) +
  # cooler colours for Pos / No / Neg
  scale_color_manual(
    name   = "sig AI",
    values = c(
      "Pos" = "#4DBBD5FF",  # teal‐blue
      "No"  = "#00A087FF",  # aqua
      "Neg" = "#3C5488FF"   # deep indigo
    )
  ) +
  # manual fill colours for SNV_varity
  scale_fill_manual(
    name   = "",
    values = c(
      "High aa-impact missense" = "#e3d579",
      "ES"                       = "#23cee9"
    )
  ) +
  # tell ggside to collapse all y-side panels into one on the right:
  ggside(
    collapse = "y",    # collapse all y-side panels into one
    y.pos    = "right" # place that single panel at the right
  ) +
  guides(shape = "none") +
  # add the y-axis density **only** in the "reg-AI" column
  geom_ysidedensity(data        = df_plot2 %>% filter(AI_source == "reg-AI")
                    ,mapping = aes(fill = SNV_varity, shape = NULL, color = NULL), 
                    alpha = 0.7, 
                    size = 0.3, 
                    position = "identity") +
  theme(
    ggside.panel.scale = 0.3,
    # shrink y-side panel to 15% of main panel width:
    axis.text.x = element_text(size = 8),
    ggside.panel.scale.y = 0.15,
    legend.position    = "top",
    legend.box         = "horizontal",
    strip.background   = element_rect(fill = "white")
  )

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, final_plot, width = 165, height = 125, units = "mm") 

# # Save
write.table(df_plot2, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig3/Fig3E.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(df_plot2, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig3/Fig3E.RData")

# df <- TCGA_AI_table_final %>%
#   filter(Gene.Symbol %in% "IDH1", SNV_varity != "stopgain")
# prop.table(table(df$all_CNA_AI))*100
# prop.table(table(df$all_reg_AI))*100
# table(df$reg_AI_type, df$SNV_varity)

######################################################
################### SUPP. FIG S5 #####################
######################################################

##### Reviewer #2 Question 3: TP53 and TSGs LOH #####

# path <- "/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_AI_final_tables/TCGA_AI_table_bb_models_final_2.RData"
# TCGA_AI_table_final <- readRDS(path)
dim(TCGA_AI_table_final)

setDT(TCGA_AI_table_final)

TCGA_AI_table_final[, `:=`(
  # assign major/minor based on whether mutant fraction > 0.5
  MUT_CNA2 = ifelse(
    MUT_CNA_fraction > 0.5,
    # estimate_MUT_CNA >= 2,
    Major_Copy_Number,
    Minor_Copy_Number
  ),
  WT_CNA2  = ifelse(
    MUT_CNA_fraction > 0.5,
    # estimate_MUT_CNA >= 2,
    Minor_Copy_Number,
    Major_Copy_Number
  )
)]

# cor.test(TCGA_AI_table_final$estimate_MUT_CNA,TCGA_AI_table_final$WT_CNA2)
# cor.test(TCGA_AI_table_final$estimate_MUT_CNA,TCGA_AI_table_final$MUT_CNA2)
# cor.test(TCGA_AI_table_final$WT_CNA,TCGA_AI_table_final$WT_CNA2)
# cor.test(TCGA_AI_table_final$MUT_CNA,TCGA_AI_table_final$MUT_CNA2)

# Choose genes to plot:

# Significant genes
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/bb_models/AI_dNdES_test/AI_dNdES_test_res_geneset_A.txt")
AI_dNdES_test_res <- read.table(output_path, header = TRUE, sep = "\t")

AI_dNdES_test_res <- AI_dNdES_test_res %>%
  filter(Z_value <= 15 & Z_value >= -15)

AI_dNdES_test_res <- AI_dNdES_test_res %>%
  group_by(Model) %>%
  mutate(
    FDR_two_sided = p.adjust(P_value, method = "fdr"),
    FDR_higher = p.adjust(Higher_one_sided_p, method = "fdr"),
    FDR_lower = p.adjust(Lower_one_sided_p, method = "fdr")
  ) %>%
  ungroup()

genes <- AI_dNdES_test_res %>%
  filter(Model == "RNA_null_model") %>%
  filter(FDR_higher < 0.1) %>%
  select(Gene,Gene_ID) %>%
  filter(Gene_ID %in% c("TP53","tsg")) %>%
  pull(Gene) %>% as.character()

# genes <- c("TP53","KRAS","CDKN2A")

# 2) Subset & extract codons for all five genes, preserving factor order
df_plot2 <- TCGA_AI_table_final %>%
  filter(Gene.Symbol %in% genes, SNV_varity != "stopgain") %>%
  # make sure Gene.Symbol is a factor in the order you want
  mutate(Gene.Symbol = factor(Gene.Symbol, levels = genes)) %>%
  pivot_longer(
    cols      = c(reg_AI_type, mRNA_AI_type, CNA_AI_type),
    names_to  = "AI_type_category",
    values_to = "AI_type"
  ) %>%
  mutate(
    AI_status = word(AI_type, 1),
    AI_source = case_when(
      AI_type_category == "reg_AI_type"  ~ "reg-AI",
      AI_type_category == "mRNA_AI_type" ~ "mRNA-AI",
      AI_type_category == "CNA_AI_type"  ~ "CNA-AI"
    )
  ) %>% 
  mutate(
    SNV_varity = recode(
      SNV_varity,
      "nonsynonymous SNV" = "High aa-impact missense",
      "effectively_syn"   = "ES"
    )
  )
dim(df_plot2)
data.frame(df_plot2[1:5,c(1:10,395:400)])

# Create the 6 bins based on your criteria
df_plot2 <- df_plot2 %>%
  mutate(
    WT_bin = case_when(
      WT_CNA2 == 0 ~ "del",
      WT_CNA2 == 1 & MUT_CNA2 == 1 ~ "neutral", 
      WT_CNA2 > 1 ~ "amp"
    ),
    MUT_bin = case_when(
      MUT_CNA2 == 0 ~ "del",
      MUT_CNA2 == 1 ~ "neutral",
      MUT_CNA2 > 1 ~ "amp"
    ),
    # New MUT_CNA_fraction bins
    CNFmut_bin = case_when(
      MUT_CNA_fraction < 0.4 ~ "<0.4",
      MUT_CNA_fraction >= 0.4 & MUT_CNA_fraction <= 0.6 ~ "0.4-0.6",
      MUT_CNA_fraction > 0.6 ~ ">0.6"
    ),
    CNFmut_bin = factor(CNFmut_bin, levels = c("<0.4", "0.4-0.6", ">0.6")),    # Combine them
    CNA_bin = paste0("WT ",WT_bin,"/ MUT ",MUT_bin),
    # Create ordered factors
    WT_bin = factor(WT_bin, levels = c("del", "neutral", "amp")),
    MUT_bin = factor(MUT_bin, levels = c("del", "neutral", "amp")),
    # CNA_bin = factor(CNA_bin, levels = c(
    #   "WT del / MUT del", "WT del / MUT neut", "WT del / MUT amp",
    #   "WT neut / MUT del", "WT neut / MUT neut", "WT neut / MUT amp", 
    #   "WT amp / MUT del", "WT amp / MUT neut", "WT amp / MUT amp"
    # ))
  ) %>%
  filter(!str_detect(CNA_bin, "NA"))
table(df_plot2$CNA_bin)
table(df_plot2$WT_bin)
table(df_plot2$MUT_bin)
table(df_plot2$CNFmut_bin)

# 3) Plot with a Gene.Symbol × AI_source grid
final_plot <- ggplot(df_plot2, aes(
    x     = WT_bin,
    y     = RNA_ASE,
    # shape = AI_status
  )) +
  # Single boxplot per bin (no color grouping)
  geom_boxplot() +
  # Colored points overlaid on boxplots
  # Colored and shaped points overlaid on boxplots
  geom_point(aes(color = AI_status, shape = AI_status), 
             size = 1.5, alpha = 0.75, position = position_jitter(width = 0.2)) +
  scale_shape_manual(
    values = c("Neg" = 25, "No" = 4, "Pos" = 19)
  ) +
  labs(x = expression(CNF[mut]), y = "RNA-VAF", shape = "", color = "") +
  theme_bw() +
  # **THIS** is the key: one row per gene, one column per AI source
  facet_grid(
    rows = vars(Gene.Symbol),
    cols = vars(AI_source),
    scales = "free_x",
    space  = "free_x"
  ) +
  # cooler colours for Pos / No / Neg
  scale_color_manual(
    name   = "sig AI",
    values = c(
      "Pos" = "#4DBBD5FF",  # teal‐blue
      "No"  = "#00A087FF",  # aqua
      "Neg" = "#3C5488FF"   # deep indigo
    )
  ) +
  # manual fill colours for SNV_varity
  scale_fill_manual(
    name   = "SNV varity",
    values = c(
      "High aa-impact missense" = "#e3d579",
      "ES"                       = "#23cee9"
    )
  ) +
  # tell ggside to collapse all y-side panels into one on the right:
  ggside(
    collapse = "y",    # collapse all y-side panels into one
    y.pos    = "right" # place that single panel at the right
  ) +
  guides(shape = "none") +
  # add the y-axis density **only** in the "reg-AI" column
  geom_ysidedensity(data        = df_plot2 %>% filter(AI_source == "reg-AI")
                    ,mapping = aes(fill = SNV_varity, shape = NULL, color = NULL), 
                    alpha = 0.7, 
                    size = 0.3, 
                    position = "identity") +
  theme(
    ggside.panel.scale = 0.3,
    # shrink y-side panel to 15% of main panel width:
    ggside.panel.scale.y = 0.15,
    legend.position    = "top",
    legend.box         = "horizontal",
    strip.background   = element_rect(fill = "white")
  )

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, final_plot, width = 200, height = 350, units = "mm") 

write.table(df_plot2, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig5/panel_A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(df_plot2, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig5/panel_A.RData")

# For the Manuscript
df_plot2 %>%
  # filter(Gene.Symbol == "TP53") %>%
  filter(WT_bin == "del") %>%
  group_by(SNV_varity) %>%
  summarize(median(RNA_ASE))

df_plot2 %>%
  # filter(Gene.Symbol == "TP53") %>%
  filter(WT_CNA2  >= 1) %>%
  group_by(SNV_varity) %>%
  summarize(median(RNA_ASE))

# From within LOH cases, % of positive reg-AI vs CNA-AI and mRNA-AI?
df <- df_plot2 %>%
  # filter(Gene.Symbol == "TP53") %>%
  # filter(MUT_CNA2 == 1) %>%
  filter(AI_source == "CNA-AI") %>%
  # filter(AI_type != "No reg-AI") %>%
  filter(WT_bin == "del") %>%
  mutate(AI_type == as.character(AI_type))
df$AI_type <- as.character(df$AI_type)
round(prop.table(table(df$AI_type,df$WT_bin)),4)*100

# From outside LOH cases, % of positive reg-AI vs CNA-AI and mRNA-AI?
df <- df_plot2 %>%
  # filter(Gene.Symbol == "TP53") %>%
  # filter(MUT_CNA2 == 1) %>%
  filter(AI_source == "mRNA-AI") %>%
  # filter(AI_type != "No reg-AI") %>%
  filter(WT_bin != "del") %>%
  mutate(AI_type == as.character(AI_type))
df$AI_type <- as.character(df$AI_type)

round(prop.table(table(df$AI_type,df$WT_bin)[,2] + table(df$AI_type,df$WT_bin)[,3]),4)*100

######################################################
################### SUPP. FIG S6 #####################
######################################################

##### Reviewer #2 Question 6: KRAS results #####

KRAS_AI <- TCGA_AI_table_final %>%
  filter(Gene.Symbol == "KRAS") %>%
  filter(SNV_varity != "stopgain")
dim(KRAS_AI)
sort(table(KRAS_AI$cancer))
table(KRAS_AI$SNV_varity)
table(KRAS_AI$MUT_CNA,KRAS_AI$WT_CNA)
  
# 1) Is there any sig AI in passengers (ES) vs missense?
kras_AI <- KRAS_AI %>%
  filter(SNV_varity == "nonsynonymous SNV") %>% # nonsynonymous SNV // effectively_syn
  select(all_mRNA_AI, all_reg_AI, all_CNA_AI)
# Calculate proportions for each AI type
prop.table(table(kras_AI$all_mRNA_AI)) * 100
prop.table(table(kras_AI$all_reg_AI)) * 100
prop.table(table(kras_AI$all_CNA_AI)) * 100

# 2) CN gain vs CN neutral AI-dNdES test to show effect of CN
KRAS_AI <- KRAS_AI %>%
  # filter(MUT_CNA2 == 1 & WT_CNA2 == 1)
  filter(MUT_CNA2 > 1 | WT_CNA2 > 1)
dim(KRAS_AI)
table(KRAS_AI$SNV_varity)

a <- sort(table(as.character(KRAS_AI$AAChange)))
a

RNA_bb_model1 <- as.formula(
  "cbind(RNA_alt, RNA_ref) ~ MUT_CNA_fraction + sqrt(gene_exp) + factor(SNV_varity)"
)

model <- VGAM::vglm(
  formula = RNA_bb_model1,
  family = VGAM::betabinomial,
  na.action = na.exclude,
  data = KRAS_AI
)

coef_table <- summary(model)@coef3
row_id <- "factor(SNV_varity)nonsynonymous SNV"

if (row_id %in% rownames(coef_table)) {
  est <- coef_table[row_id, "Estimate"]
  z <- coef_table[row_id, "z value"]
  p <- 2 * pnorm(-abs(z))
  p_high <- pnorm(z, lower.tail = FALSE)
  p_low <- pnorm(z, lower.tail = TRUE)

  AI_dNdES_test_res <- data.frame(
    Gene = "KRAS",
    Estimate = est,
    Z_value = z,
    P_value = p,
    Higher_one_sided_p = p_high,
    Lower_one_sided_p = p_low
  )
}
AI_dNdES_test_res

# 3) Plot for selected oncogenes: CN MUT del / neutral (MUT/WT) / CN MUT gain

# Significant genes
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/bb_models/AI_dNdES_test/AI_dNdES_test_res_geneset_A.txt")
AI_dNdES_test_res <- read.table(output_path, header = TRUE, sep = "\t")

AI_dNdES_test_res <- AI_dNdES_test_res %>%
  filter(Z_value <= 15 & Z_value >= -15)

AI_dNdES_test_res <- AI_dNdES_test_res %>%
  group_by(Model) %>%
  mutate(
    FDR_two_sided = p.adjust(P_value, method = "fdr"),
    FDR_higher = p.adjust(Higher_one_sided_p, method = "fdr"),
    FDR_lower = p.adjust(Lower_one_sided_p, method = "fdr")
  ) %>%
  ungroup()

genes <- AI_dNdES_test_res %>%
  filter(Model == "RNA_null_model") %>% # RNA_null_model // RNA_bb_model1
  filter(FDR_higher < 0.25) %>%
  select(Gene,Gene_ID) %>%
  filter(Gene_ID %in% c("oncogene")) %>%
  pull(Gene) %>% as.character()
genes
# genes <- c("TP53","KRAS","CDKN2A")

# 2) Subset & extract codons for all five genes, preserving factor order
df_plot2 <- TCGA_AI_table_final %>%
  filter(Gene.Symbol %in% genes, SNV_varity != "stopgain") %>%
  # make sure Gene.Symbol is a factor in the order you want
  mutate(Gene.Symbol = factor(Gene.Symbol, levels = genes)) %>%
  pivot_longer(
    cols      = c(reg_AI_type, mRNA_AI_type, CNA_AI_type),
    names_to  = "AI_type_category",
    values_to = "AI_type"
  ) %>%
  mutate(
    AI_status = word(AI_type, 1),
    AI_source = case_when(
      AI_type_category == "reg_AI_type"  ~ "reg-AI",
      AI_type_category == "mRNA_AI_type" ~ "mRNA-AI",
      AI_type_category == "CNA_AI_type"  ~ "CNA-AI"
    )
  ) %>% 
  mutate(
    SNV_varity = recode(
      SNV_varity,
      "nonsynonymous SNV" = "High aa-impact missense",
      "effectively_syn"   = "ES"
    )
  )
dim(df_plot2)
data.frame(df_plot2[1:5,c(1:10,395:400)])

# Create the 6 bins based on your criteria
df_plot2 <- df_plot2 %>%
  mutate(
    WT_bin = case_when(
      WT_CNA2 == 0 ~ "del",
      WT_CNA2 == 1 & MUT_CNA2 == 1 ~ "neutral", 
      WT_CNA2 > 1 ~ "amp"
    ),
    MUT_bin = case_when(
      MUT_CNA2 == 0 ~ "del",
      WT_CNA2 == 1 & MUT_CNA2 == 1 ~ "neutral",
      MUT_CNA2 > 1 ~ "amp"
    ),
    # New MUT_CNA_fraction bins
    CNFmut_bin = case_when(
      MUT_CNA_fraction < 0.4 ~ "<0.4",
      MUT_CNA_fraction >= 0.4 & MUT_CNA_fraction <= 0.6 ~ "0.4-0.6",
      MUT_CNA_fraction > 0.6 ~ ">0.6"
    ),
    CNFmut_bin = factor(CNFmut_bin, levels = c("<0.4", "0.4-0.6", ">0.6")),    # Combine them
    CNA_bin = paste0("WT ",WT_bin,"/ MUT ",MUT_bin),
    # Create ordered factors
    WT_bin = factor(WT_bin, levels = c("del", "neutral", "amp")),
    MUT_bin = factor(MUT_bin, levels = c("del", "neutral", "amp")),
    # CNA_bin = factor(CNA_bin, levels = c(
    #   "WT del / MUT del", "WT del / MUT neut", "WT del / MUT amp",
    #   "WT neut / MUT del", "WT neut / MUT neut", "WT neut / MUT amp", 
    #   "WT amp / MUT del", "WT amp / MUT neut", "WT amp / MUT amp"
    # ))
  ) %>%
  filter(!str_detect(CNA_bin, "NA"))
table(df_plot2$CNA_bin)
table(df_plot2$WT_bin)
table(df_plot2$MUT_bin)
table(df_plot2$CNFmut_bin)

# 3) Plot with a Gene.Symbol × AI_source grid
final_plot <- ggplot(df_plot2, aes(
    x     = MUT_bin,
    y     = log2(RNA_ASE/DNA_ASE),
    # shape = AI_status
  )) +
  # Single boxplot per bin (no color grouping)
  geom_boxplot(outlier.shape = NA) +
  # Colored points overlaid on boxplots
  # Colored and shaped points overlaid on boxplots
  geom_point(aes(color = AI_status, shape = AI_status, 
             size = AI_status, alpha = AI_status), position = position_jitter(width = 0.2)) +
  scale_shape_manual(
    values = c("Neg" = 25, "No" = 4, "Pos" = 19)
  ) +
  scale_size_manual(
    values = c(
      "Pos" = 1.5,
      "No"  = 0.75,           # override “No” to 0.5
      "Neg" = 1.5
    ),
    guide = "none"           # hide redundant legend
  ) +
  scale_alpha_manual(
    values = c(
      "Pos" = 0.75,
      "No"  = 0.25,          # override “No” to 0.25
      "Neg" = 0.75
    ),
    guide = "none"
  ) +
  # labs(x = expression(CNF[mut]), y = "RNA-VAF", shape = "", color = "") +
  labs(x = "MUT bins", y = "RNA-VAF", shape = "", color = "") +
  theme_bw() +
  # **THIS** is the key: one row per gene, one column per AI source
  facet_grid(
    rows = vars(Gene.Symbol),
    cols = vars(AI_source),
    scales = "free_x",
    space  = "free_x"
  ) +
  # cooler colours for Pos / No / Neg
  scale_color_manual(
    name   = "sig AI",
    values = c(
      "Pos" = "#4DBBD5FF",  # teal‐blue
      "No"  = "grey80",  # aqua
      "Neg" = "#3C5488FF"   # deep indigo
    )
  ) +
  # manual fill colours for SNV_varity
  scale_fill_manual(
    name   = "SNV varity",
    values = c(
      "High aa-impact missense" = "#e3d579",
      "ES"                       = "#23cee9"
    )
  ) +
  # coord_cartesian(ylim = c(0,2))+
  # tell ggside to collapse all y-side panels into one on the right:
  ggside(
    collapse = "y",    # collapse all y-side panels into one
    y.pos    = "right" # place that single panel at the right
  ) +
  guides(shape = "none") +
  # add the y-axis density **only** in the "reg-AI" column
  geom_ysidedensity(data        = df_plot2 %>% filter(AI_source == "reg-AI")
                    ,mapping = aes(fill = SNV_varity, shape = NULL, color = NULL), 
                    alpha = 0.7, 
                    size = 0.3, 
                    position = "identity") +
  theme(
    ggside.panel.scale = 0.3,
    # shrink y-side panel to 15% of main panel width:
    ggside.panel.scale.y = 0.15,
    legend.position    = "top",
    legend.box         = "horizontal",
    strip.background   = element_rect(fill = "white")
  )

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, final_plot, width = 200, height = 350, units = "mm") 

write.table(df_plot2, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig6/panel_A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(df_plot2, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig6/panel_A.RData")

# For the Manuscript
df_plot2 %>%
  # filter(Gene.Symbol == "KRAS") %>%
  filter(MUT_bin == "amp") %>%
  group_by(SNV_varity) %>%
  summarize(median(RNA_ASE))

df_plot2 %>%
  # filter(Gene.Symbol == "KRAS") %>%
  filter(MUT_bin  == "neutral") %>%
  group_by(SNV_varity) %>%
  summarize(median(RNA_ASE))

# From within CNamp cases, % of positive reg-AI vs CNA-AI and mRNA-AI?
df <- df_plot2 %>%
  # filter(Gene.Symbol == "TP53") %>%
  filter(AI_source == "CNA-AI") %>%
  # filter(AI_type != "No reg-AI") %>%
  filter(MUT_bin == "amp") %>%
  mutate(AI_type == as.character(AI_type))
df$AI_type <- as.character(df$AI_type)
round(prop.table(table(df$AI_type,df$MUT_bin)),4)*100

# From outside LOH cases, % of positive reg-AI vs CNA-AI and mRNA-AI?
df <- df_plot2 %>%
  # filter(Gene.Symbol == "TP53") %>%
  # filter(MUT_CNA2 == 1) %>%
  filter(AI_source == "mRNA-AI") %>%
  # filter(AI_type != "No reg-AI") %>%
  filter(WT_bin != "del") %>%
  mutate(AI_type == as.character(AI_type))
df$AI_type <- as.character(df$AI_type)

round(prop.table(table(df$AI_type,df$WT_bin)[,2] + table(df$AI_type,df$WT_bin)[,3]),4)*100

###################################################################################################
##################################### SURVIVAL CURVES #############################################
###################################################################################################

#Do not activate any conda
library(ggplot2)
library(gridExtra)
library(stringr)
library(ggsignif)
library(tidyr)
library(dplyr)
library(viridis)
library(ggrepel)
library(survival)
library(survminer)
library(forestmodel) # This is the package I don't have in R_figures conda env
library(readxl)
library(purrr)
library(glue)

ggsave_workaround <- function(g){survminer:::.build_ggsurvplot(x = g,
                                                               surv.plot.height = NULL,
                                                               risk.table.height = NULL,
                                                               ncensor.plot.height = NULL)}   

# Significant genes
output_path <- paste0("/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/bb_models/AI_dNdES_test/AI_dNdES_test_res_geneset_A.txt")
AI_dNdES_test_res <- read.table(output_path, header = TRUE, sep = "\t")

AI_dNdES_test_res <- AI_dNdES_test_res %>%
  filter(Z_value <= 15 & Z_value >= -15)

AI_dNdES_test_res <- AI_dNdES_test_res %>%
  group_by(Model) %>%
  mutate(
    FDR_two_sided = p.adjust(P_value, method = "fdr"),
    FDR_higher = p.adjust(Higher_one_sided_p, method = "fdr"),
    FDR_lower = p.adjust(Lower_one_sided_p, method = "fdr")
  ) %>%
  ungroup()

# TCGA AI table
TCGA_AI_table_final <- TCGA_AI_table_final %>%
  filter(SNV_varity != "stopgain")

# TCGA metadata
TCGA_metadata <- read.table("/g/strcombio/fsupek_cancer1/gpalou/TCGA_metadata/TCGA_clinical_all.tsv", 
                header = TRUE, sep = "\t", stringsAsFactors = FALSE)
TCGA_metadata <- TCGA_metadata[,colnames(TCGA_metadata) %in% c("submitter_id", "gender", "race", "age_at_diagnosis", "vital_status", "days_to_death", "days_to_last_follow_up", "tumor_stage")]
colnames(TCGA_metadata) <- c("sample", "gender", "race", "tumor_stage","age_at_diagnosis", "vital_status", "days_to_death", "days_to_last_follow_up")
project_names <- read.table("/g/strcombio/fsupek_cancer1/gpalou/TCGA_RNAseq_quantification/TCGA_projects_names.txt")

#####################################################################
################### FIGURE 4A and SUPP. FIG S7A #####################
#####################################################################

# A) Survival curves

build_surv_df <- function(
  AI_dNdES_test_res,
  TCGA_AI_table_final,
  TCGA_metadata,
  gene_selected = NULL,
  cancer_type = NULL,
  AI_type,        
  sel_FDR,        
  gene_set,       
  var_type,       
  pos_threshold,  
  neg_threshold,  
  n_pos_bins     
) {
  # A) Select genes by FDR and model
  model_name <- if (AI_type == "reg-AI") "RNA_bb_model1" else "RNA_null_model"
  if (is.null(gene_selected)) {
  sel_genes <- AI_dNdES_test_res %>%
    filter(Model == model_name, FDR_higher < sel_FDR) %>%
    pull(Gene) %>% as.character()
  } else {
    sel_genes <- gene_selected
  }

  if (gene_set == "noTP53") {
    sel_genes <- setdiff(sel_genes, "TP53")
  } else if (gene_set == "onlyTP53") {
    sel_genes <- intersect(sel_genes, "TP53")
  }

  TCGA_AI_table_final_filt <- TCGA_AI_table_final %>%
    filter(
      SNV_varity == "nonsynonymous SNV",
      Gene.Symbol %in% sel_genes,
      # keep all cancers when pancancer, otherwise only the chosen one
      cancer_type == "pancancer" | cancer == cancer_type
    )

  # B) Build AI per sample exactly as in your original code
  AI_affected <-  TCGA_AI_table_final_filt %>%
    group_by(sample) %>%
    mutate(
      RNA_bb_model1_AI_diff   = RNA_observed_ratio - RNA_expected_ratio_bb_model1,
      RNA_null_model_AI_diff  = RNA_observed_ratio - RNA_expected_ratio_null_model,
      DNA_null_model_AI_diff  = DNA_observed_ratio - DNA_expected_ratio_null_model,
      DNA_bb_model1_AI_diff   = DNA_observed_ratio - DNA_expected_ratio_bb_model1,
      DNA_AI_diff             = DNA_null_model_AI_diff - DNA_bb_model1_AI_diff
    ) %>%
    summarise(
      AI = if (var_type == "qual") {
        if (AI_type == "reg-AI") {
          # any positive reg-AI in that sample
          as.numeric(any(reg_AI_type == "Pos reg-AI"))
        } else if (AI_type == "mRNA-AI") {
          as.numeric(any(mRNA_AI_type == "Pos mRNA-AI"))
        }
      } else {
        if (AI_type == "reg-AI") {
          mean(RNA_bb_model1_AI_diff, na.rm = TRUE)
        } else if (AI_type == "mRNA-AI") {
          mean(RNA_null_model_AI_diff, na.rm = TRUE)
        } else if (AI_type == "CNA-AI") {
          mean(DNA_AI_diff, na.rm = TRUE)
        }
      },
      .groups = "drop"
    )

  # C) Merge to get survival times
  surv_df <- TCGA_metadata %>%
    inner_join(AI_affected, by = "sample") %>%
    mutate(
      days_to_death          = as.numeric(days_to_death),
      days_to_last_follow_up = as.numeric(days_to_last_follow_up),
      time   = if_else(vital_status == "dead", days_to_death, days_to_last_follow_up) / 365.25,
      status = as.integer(vital_status == "dead")
    ) %>%
    filter(!is.na(time))

  if (var_type == "quant") {
    # 1) Get the global min
    min_AI <- min(surv_df$AI, na.rm = TRUE)

    # 2) Compute quantiles on the positives
    pos_vals   <- surv_df$AI[surv_df$AI > pos_threshold]
    pos_breaks <- unname(quantile(
      pos_vals,
      probs = seq(0, 1, length.out = n_pos_bins + 1),
      na.rm = TRUE
    ))

    # 3) Glue together your breaks
    # breaks_new <- c(min_AI, neg_threshold, pos_threshold, pos_breaks[-1])
    breaks_new <- c(neg_threshold, pos_threshold, pos_breaks[-1])

    # 4) Build the exact same labels
    neg_ref_labels <- c(
      sprintf("[%s,%s)", round(min_AI, 2), neg_threshold),
      sprintf("[%s,%s]",    neg_threshold,      pos_threshold)
    )
    lower <- round(head(pos_breaks, -1), 2)
    upper <- round(tail(pos_breaks, -1), 2)
    pos_labels <- paste0("(", lower, ",", upper, "]")

    all_labels <- unname(c(neg_ref_labels, pos_labels))
    all_labels <- all_labels[-1]

    # 5) Cut (right=FALSE, include.lowest=TRUE) and refactor
    surv_df$AI <- cut(
      surv_df$AI,
      breaks         = breaks_new,
      labels         = all_labels,
      include.lowest = TRUE,
      right          = FALSE
    )
    surv_df$AI <- factor(
      surv_df$AI,
      levels = c(neg_ref_labels[2], setdiff(all_labels, neg_ref_labels[2]))
    )
  }

  message("Bin counts:\n")
  print(table(surv_df$AI))
  return(surv_df)
}

surv_df <- build_surv_df(
  AI_dNdES_test_res      = AI_dNdES_test_res,
  TCGA_AI_table_final    = TCGA_AI_table_final,
  TCGA_metadata          = TCGA_metadata,
  cancer_type            = "pancancer",
  gene_selected          = NULL,
  AI_type                = "mRNA-AI",
  sel_FDR                = 0.25,
  gene_set               = "noTP53",
  var_type               = "quant",
  pos_threshold          = 0.1,
  neg_threshold          = -0.1,
  n_pos_bins             = 3
)

# Fit and Plot
km_fit <- survfit(Surv(time, status) ~ AI, data = surv_df)

plot <- ggsurvplot(
  km_fit, data = surv_df,
  risk.table       = TRUE,
  pval             = TRUE,
  conf.int         = FALSE,
  xlab             = "Time (years)",
  ylab             = "Survival probability",
  # legend.labs      = c("Not affected","Affected"),
  legend.title     = "AI",
  xlim = c(0,15),
  surv.median.line = "v",
  ggtheme          = theme_bw(base_size = 12)
  )

plot_to_save <- ggsave_workaround(plot)

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, plot_to_save, width = 175, height = 180, units = "mm") 

# Save
write.table(surv_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig4/Fig4A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(surv_df, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig4/Fig4A.RData")

surv_df <- build_surv_df(
  AI_dNdES_test_res      = AI_dNdES_test_res,
  TCGA_AI_table_final    = TCGA_AI_table_final,
  TCGA_metadata          = TCGA_metadata,
  cancer_type            = "pancancer",
  gene_selected          = NULL,
  AI_type                = "reg-AI",
  sel_FDR                = 0.25,
  gene_set               = "noTP53",
  var_type               = "quant",
  pos_threshold          = 0.1,
  neg_threshold          = -0.1,
  n_pos_bins             = 3
)

write.table(surv_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig7/panel_A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(surv_df, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig7/panel_A.RData")

# Manuscript --> P value (log-rank)
surv_pvalue(km_fit,data = surv_df)
km_fit

# Median survival time
summary(km_fit)$table[1,"median"]
mean(summary(km_fit)$table[-c(1:2),"median"],na.rm=TRUE)

#####################################################################
################### FIGURE 4B and SUPP. FIG S7B #####################
#####################################################################

# B) Cox Regression

# Covariates
covariates <- read.table("/g/strcombio/fsupek_home/gpalou/analysis_results/NMD_project/cancers/pancancer/NMD_efficiencies_TCGA.txt", header = TRUE, sep = "\t")
covariates <- covariates[,colnames(covariates) %in% c("sample", "endogenous_purity", "CNV_burden", "endogenous_LF", "cancer_subtype","MSI_status")]
data_path <- "/g/strcombio/fsupek_cancer1/gpalou/ASE_project"
heterogeneity <- read_excel(paste0(data_path,"/conversor_tables/Erik_van_Dijk_et_al_2021_CNH.xlsx"))
heterogeneity <- as.data.frame(heterogeneity)[,c(1,2,13)]
heterogeneity$sample <- substr(heterogeneity$Samplename, start = 1, stop = 12)
heterogeneity$CNH <- as.numeric(heterogeneity$CNH)
CNA_signatures <- read.table("/g/strcombio/fsupek_cancer1/gpalou/TCGA_CNV_PCA/pancancer_sparse_PCA_ind_3e-04_robust_no_num_PCs_100.txt", header = TRUE, sep = "\t")
cols <- colnames(CNA_signatures)[which( colSums(CNA_signatures) != 0 )]
CNA_signatures <- CNA_signatures[,cols]
CNA_signatures$sample <- gsub("\\.", "-", rownames(CNA_signatures))

HR_cox_model <- function( sel_FDR, gene_set, AI_type, var_type, cancer_type, gene_selected = NULL,
                          AI_table, TCGA_metadata, covariate_df, heterogeneity, CNA_signatures
                        ) {
  # A) select genes by FDR_higher < sel_FDR

  if (is.null(gene_selected)) {
    sel_genes <- AI_dNdES_test_res %>%
      filter(
        Model == if_else(AI_type == "reg-AI", "RNA_bb_model1","RNA_null_model"),
        FDR_higher < sel_FDR
      ) %>%
      pull(Gene) %>% as.character()
  } else {
    sel_genes <- gene_selected
  }
 
  # apply gene_set
  if (gene_set=="noTP53")  sel_genes <- setdiff(sel_genes, "TP53")
  if (gene_set=="onlyTP53") sel_genes <- intersect(sel_genes, "TP53")
  
  # B) subset AI data & compute AI_diff
  df <- AI_table %>%
    # filter(WT_CNA2 =0) %>% # LOH-only
    filter(SNV_varity=="nonsynonymous SNV",
           Gene.Symbol %in% sel_genes) %>%
    mutate(
      RNA_bb_model1_AI_diff = RNA_observed_ratio - RNA_expected_ratio_bb_model1,
      RNA_null_model_AI_diff = RNA_observed_ratio - RNA_expected_ratio_null_model,
      DNA_null_model_AI_diff = DNA_observed_ratio - DNA_expected_ratio_null_model,
      DNA_bb_model1_AI_diff = DNA_observed_ratio - DNA_expected_ratio_bb_model1,
      DNA_AI_diff = DNA_null_model_AI_diff - DNA_bb_model1_AI_diff
    )
  if(nrow(df)==0) return(NULL)
  # C) Cancer type
  if (cancer_type == "pancancer") {
      df <- df
      TCGA_metadata_filt <- TCGA_metadata
  } else {
      df <- df %>%
              filter(cancer %in% cancer_type)
      TCGA_metadata_filt <- TCGA_metadata %>%
              filter(sample %in% as.character(df$sample))
  }
  
  # C) make per-sample ASE calls
  AI_by_sample <- df %>%
    group_by(sample,cancer) %>%
    summarize(
      AI = if (var_type=="qual") {
        if (AI_type=="reg-AI") {
          (reg_AI_type == "Pos reg-AI") #& (mRNA_AI_type == "Pos mRNA-AI")
        } else if (AI_type=="mRNA-AI") {
          mRNA_AI_type == "Pos mRNA-AI"
        } else if (AI_type == "CNA-AI") {
          (CNA_AI_type == "Pos CNA-AI") #& (mRNA_AI_type == "Pos mRNA-AI")
        }
      } else {
        if (AI_type=="reg-AI") {
          mean(RNA_bb_model1_AI_diff)
        } else if (AI_type=="mRNA-AI") {
          mean(RNA_null_model_AI_diff)
        } else if (AI_type=="CNA-AI") {
          mean(DNA_null_model_AI_diff)
        } 
      },
      .groups="drop"
    )

  # D) merge clinical + covariates

  surv_df <- TCGA_metadata_filt %>%
    inner_join(AI_by_sample, by="sample") %>%
    left_join(covariate_df,   by="sample") %>%
    left_join(heterogeneity,  by="sample") %>%
    left_join(CNA_signatures, by="sample") %>%
    mutate(
      days_to_death          = as.numeric(na_if(days_to_death,          "--")),
      days_to_last_follow_up = as.numeric(na_if(days_to_last_follow_up, "--")),
      age_at_diagnosis       = as.numeric(na_if(age_at_diagnosis,       "--")),
      time   = if_else(vital_status=="dead", days_to_death, days_to_last_follow_up) / 365.25,
      status = as.integer(vital_status=="dead")
    ) %>%
    filter(!is.na(time))

  if(nrow(surv_df)<10) return(NULL)
  if (var_type == "qual") { 
    if ( sum(surv_df$AI,na.rm=TRUE)<5) {return(NULL)}
  }

  # NEW VARIABLE CREATION:
  if (cancer_type == "pancancer" & length(gene_selected) > 1) {
    surv_df$cancer <- relevel(factor(surv_df$cancer), ref = "TCGA-BRCA")
  } else if ((cancer_type == "pancancer" & length(gene_selected) == 1) ) {
    ref_cancer <- names(rev(sort(table(as.character(surv_df$cancer))))[1])
    surv_df$cancer <- relevel(factor(surv_df$cancer), ref = ref_cancer)
  } else {
    surv_df$cancer <- factor(as.character(surv_df$cancer))
  }
  
  surv_df <- surv_df %>%
    mutate(
      age_group = factor(ntile(age_at_diagnosis/365.25, 5)),
      CNA_burden = factor(ntile(CNV_burden, 4)),
      CNhet     = factor(ntile(CNH, 4)),
      Stage     = factor( na_if(as.character(Stage), "NA") )
    ) %>%
    rename(
      Age    = age_group,
      Purity = endogenous_purity,
      Cancer = cancer,
      Gender = gender
    )
  covs <- c("AI", "Age","Gender","Purity", "CNA_burden", "CNhet", "Stage", "Cancer")
  # Remove samples with any NA in any covariate
  if (cancer_type %in% c("TCGA-LGG","TCGA-GBM","TCGA-PRAD","TCGA-SARC")) {
    surv_df <- surv_df %>%
      drop_na(any_of(covs_mod <- setdiff(covs, "Stage")))
  } else if (cancer_type == "TCGA-SKCM") {
    surv_df <- surv_df %>%
      drop_na(any_of(covs_mod <- setdiff(covs, "CNA_burden")))     
  } else if (cancer_type %in% c("TCGA-LAML", "TCGA-THYM", "TCGA-DLBC")) {
    surv_df <- surv_df %>%
      drop_na(any_of(covs_mod <- setdiff(covs, "Purity"))) 
  } else {
    surv_df <- surv_df %>%
      drop_na(any_of(covs))
  }

  # if (nrow(surv_df) == 0) { return(NULL) }

  print(
    surv_df %>%
      group_by(Age) %>%
      summarise(
        age_range = range(age_at_diagnosis / 365.25, na.rm = TRUE)
      )
  )

  if (var_type == "quant") {

    # 1) First extract the full range of your data:
    min_AI <- min(surv_df$AI, na.rm=TRUE)
    # 2) Decide how many *equal‐count* bins you want for the positives:
    n_pos_bins <- 3  
    # 3) Compute the quantile “cut points” on AI > 0.1:
    pos_vals    <- surv_df$AI[ surv_df$AI >  0.1 ]
    if(length(pos_vals) == 0){return(NULL)}
    pos_breaks  <- quantile(pos_vals,
                            probs = seq(0, 1, length.out = n_pos_bins + 1),
                            na.rm = TRUE)
    # pos_breaks is length 6: [0% ,20%,40%,60%,80%,100%] of the positive values
    # 4) Build your final breaks vector by gluing together:
    #    - the minimum (to catch everything below –0.1)
    #    - the two fixed cutpoints (–0.1, +0.1)
    #    - and the 20%,40%,…,100% quantiles of the positives
    breaks_new  <- c(-0.1, 0.1, pos_breaks[-1])  # drop the 0% (pos_breaks[1]) since 0.1 is lower
    # breaks_new  <- c(min_AI,-0.1, 0.1, pos_breaks[-1])  # drop the 0% (pos_breaks[1]) since 0.1 is lower
    # 5) Make human‐readable labels for each bin:
    neg_ref_labels <- c(
      sprintf("[%s,%s)", round(min_AI, 2), -0.1),
      sprintf("[%s,%s]",    -0.1,      0.1)
    )      #  for the positive slices, construct e.g. "(0.12,0.23]"
    pos_labels <- 
      paste0("(", 
            round(pos_breaks[-length(pos_breaks)], 2),  # lower edges: 20%,40%…,80%
            ",", 
            round(pos_breaks[-1], 2),                  # upper edges: 40%,60%…,100%
            "]")
    all_labels <- c(neg_ref_labels, pos_labels)
    all_labels <- all_labels[-1]
    # 6) Finally cut the data:
    surv_df$AI <- cut(
      surv_df$AI,
      breaks       = breaks_new,
      labels       = all_labels,
      include.lowest = TRUE,
      right          = FALSE
    )
    # Reorder levels so that "[-0.1,0.1)" is the reference level
    surv_df$AI <- factor(
      surv_df$AI,
      levels = c("[-0.1,0.1]", setdiff(all_labels, "[-0.1,0.1]"))
    )
  }
  print(table(surv_df$AI))
  
  # E) define your Cox covariates (now using age_group, cancer, CNA_quartiles, CNH_quartiles, Stage)
    # --- 0) quick exit for small cohorts (non‐pan‐cancer only) ----
      # note: cov “AI” replaces your old ASE variable
  if (nrow(surv_df) == 0) { return(NULL) }
  if (var_type == "qual") {
    if (sum(surv_df$AI, na.rm = TRUE) < 5) { return(NULL) }
  }
  
  # --- 1) compute your missingness / uniqueness flags ----
  stageNA <- (sum(is.na(surv_df$Stage) | surv_df$Stage == "NA") /
              nrow(surv_df)) > 0.9
  isGenderUnique1 <- length(table( surv_df[!is.na(surv_df$AI), "Gender"] )) == 1
  isGenderUnique2 <- length(table( surv_df[!is.na(surv_df$AI) & !is.na(surv_df$CNA_burden), "Gender"] ) ) == 1
  # these all have fixed gender
  grep_match <- grep("BRCA|UCEC|CESC|UCS|TGCT|OV|PRAD", cancer_type)
  # --- 2) start from your full pan‐cancer covariate list ----
  if (cancer_type == "pancancer") {
    covs_mod <- covs
  } else {
    covs_mod <- setdiff(covs, "Cancer")
  }
  
  # --- 3) apply each drop rule in turn ----
  # a) drop Stage if >90% missing
  if (stageNA) { covs_mod <- setdiff(covs_mod, "Stage") }
  # b) drop Gender if unique or sex‐specific cancer
  if (length(grep_match) == 1 || isGenderUnique1 || isGenderUnique2) {
    covs_mod <- setdiff(covs_mod, "Gender")
  }
  # --- 4) add PCs for multi‐gene pan‐cancer ----
  # if (cancer_type == "pancancer" && length(gene_type) > 1) {
  #   covs_mod <- c(covs_mod, dim_vars)
  # }
  # --- 5) build the formula & fit ----
  cox_fml <- as.formula(
    paste0("Surv(time, status) ~ ", paste(covs_mod, collapse = " + "))
  )
  # cox_fml <- as.formula(paste("Surv(time,status) ~", paste(covs, collapse=" + ")))
  # coxph(Surv(time, status) ~ AI + Age + Purity + CNA_burden + CNhet + Stage, data=surv_df, na.action=na.exclude)

  # F) fit the model
  # if (cancer_type != "pancancer") {
  surv_df <- surv_df %>% 
    group_by(Cancer) %>% 
    filter(n() >= 10, sum(status) >= 5) %>% 
    ungroup()
  # }

  surv_obj <- with(surv_df, Surv(time, status))
  cmod <- tryCatch(coxph(cox_fml, data=surv_df, na.action=na.exclude), error=function(e) NULL)

  if(is.null(cmod)) return(NULL)
  # G) extract AI hazard ratio
  s   <- summary(cmod)
  # print(s)
  idx <- grep("^AI", rownames(s$coefficients))
  if(length(idx)==0) return(NULL)
  df_res <- tibble(
    term    = rownames(s$coefficients)[idx],
    HR      = s$coefficients[idx,"exp(coef)"],
    CI_low  = s$conf.int[idx,"lower .95"],
    CI_high = s$conf.int[idx,"upper .95"],
    p_value = s$coefficients[idx,"Pr(>|z|)"],
    z_score = s$coefficients[idx,"z"]
  )

  format <- forest_model_format_options(colour = "firebrick", shape = 15, text_size = 5, point_size = 5)
  # Save
  return( list(coxmodel = cmod, format = format, AI_affected = AI_by_sample, surv_obj = surv_obj, df_res = df_res) )

}

HR_cox_res_list <- HR_cox_model(sel_FDR = 0.25, gene_set = "noTP53", 
            AI_type = "mRNA-AI", var_type = "quant", cancer_type = "pancancer",
            AI_table = TCGA_AI_table_final, 
            TCGA_metadata = TCGA_metadata, 
            covariate_df = covariates, 
            heterogeneity = heterogeneity, 
            CNA_signatures = CNA_signatures) 

plot <- forest_model(HR_cox_res_list$coxmodel, covariates = c("AI", "Gender", "Age", "Stage", "CNhet"), 
        format_options = HR_cox_res_list$format)

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, plot, width = 175, height = 180, units = "mm") 

# Save
# write.table(surv_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig4/Fig4B.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(HR_cox_res_list, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig4/Fig4B.RData")

HR_cox_res_list <- HR_cox_model(sel_FDR = 0.25, gene_set = "noTP53", 
            AI_type = "reg-AI", var_type = "quant", cancer_type = "pancancer",
            AI_table = TCGA_AI_table_final, 
            TCGA_metadata = TCGA_metadata, 
            covariate_df = covariates, 
            heterogeneity = heterogeneity, 
            CNA_signatures = CNA_signatures) 

plot <- forest_model(HR_cox_res_list$coxmodel, covariates = c("AI", "Gender", "Age", "Stage", "CNhet"), 
        format_options = HR_cox_res_list$format)

# write.table(HR_cox_res_list, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig7/panel_B.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(HR_cox_res_list, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig7/panel_B.RData")

#####################################################################
################### FIGURE 4C and SUPP. FIG S7C #####################
#####################################################################

# A) KM curve for TP53

surv_df <- build_surv_df(
  AI_dNdES_test_res      = AI_dNdES_test_res,
  TCGA_AI_table_final    = TCGA_AI_table_final,
  TCGA_metadata          = TCGA_metadata,
  AI_type                = "mRNA-AI",
  cancer_type            = "pancancer",
  gene_selected         = "TP53",
  sel_FDR                = 0.25,
  gene_set               = "all",
  var_type               = "quant",
  pos_threshold          = 0.1,
  neg_threshold          = -0.1,
  n_pos_bins             = 3
)

# Fit and Plot
km_fit <- survfit(Surv(time, status) ~ AI, data = surv_df)

plot <- ggsurvplot(
  km_fit, data = surv_df,
  risk.table       = TRUE,
  pval             = TRUE,
  conf.int         = FALSE,
  xlab             = "Time (years)",
  ylab             = "Survival probability",
  # legend.labs      = c("Not affected","Affected"),
  legend.title     = "AI",
  xlim = c(0,15),
  surv.median.line = "v",
  ggtheme          = theme_bw(base_size = 12)
  )

plot_to_save <- ggsave_workaround(plot)

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, plot_to_save, width = 175, height = 180, units = "mm") 

# Save
write.table(surv_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig4/Fig4C.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(surv_df, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig4/Fig4C.RData")

# Manuscript --> P value (log-rank)
surv_pvalue(km_fit,data = surv_df)
km_fit

# Median survival time
summary(km_fit)$table[1,"median"]
mean(summary(km_fit)$table[-c(1:2),"median"],na.rm=TRUE)

# reg-AI

surv_df <- build_surv_df(
  AI_dNdES_test_res      = AI_dNdES_test_res,
  TCGA_AI_table_final    = TCGA_AI_table_final,
  TCGA_metadata          = TCGA_metadata,
  AI_type                = "reg-AI",
  cancer_type            = "pancancer",
  gene_selected         = "TP53",
  sel_FDR                = 0.25,
  gene_set               = "all",
  var_type               = "quant",
  pos_threshold          = 0.1,
  neg_threshold          = -0.1,
  n_pos_bins             = 3
)

write.table(surv_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig7/panel_C.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(surv_df, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig7/panel_C.RData")

# B) # COX REGRESSION BY GENE TYPE 

# 1) mRNA-AI

sel_genes <- AI_dNdES_test_res %>%
  filter(
    Model == "RNA_null_model", # RNA_bb_model1 RNA_null_model
    FDR_higher < 0.25
  ) %>%
  pull(Gene) %>% as.character()

by_gene_all_cox_models <- list()
by_gene_all_cox_res <- c()

for (sel_gene in sel_genes ) {
  print(paste0("-------------",sel_gene,"--------------"))

  # prepare a place to store the warning
  warn_obj <- NULL

  # run your model, catching any warning in warn_obj
  HR_cox_model_res   <- withCallingHandlers(
    HR_cox_model(sel_FDR    = 0.25,
                gene_set   = "all",
                AI_type    = "mRNA-AI",
                var_type   = "quant",
                gene_selected = sel_gene,
                cancer_type= "pancancer",
                AI_table   = TCGA_AI_table_final,
                TCGA_metadata = TCGA_metadata,
                covariate_df = covariates,
                heterogeneity = heterogeneity,
                CNA_signatures = CNA_signatures),
    warning = function(w) {
      warn_obj <<- w            # stash the warning condition
      invokeRestart("muffleWarning")  # prevent it from printing
    }
  )

  if ( !is.null(HR_cox_model_res$df_res) ) {
    # now `res` is your list and `warn_obj` is the warning condition
    if (is.null(warn_obj)) {
      warn_msg <- "none"
    } else {
      warn_msg <- conditionMessage(warn_obj)
      print(warn_msg)
    }
    HR_cox_model_res$df_res$warn_msg <- warn_msg
    HR_cox_model_res$df_res$sel_gene <- sel_gene
    by_gene_all_cox_models[[sel_gene]] <- HR_cox_model_res
    by_gene_all_cox_res <- rbind(by_gene_all_cox_res,HR_cox_model_res$df_res)
  }
}

by_gene_all_cox_res_filt <- by_gene_all_cox_res %>%
  arrange(desc(HR)) %>%
  filter(HR < 20) %>%
  filter(warn_msg == "none") %>%
  filter(CI_high < 20 & CI_low < 20) %>%
  filter(!is.infinite(CI_low) & !is.infinite(CI_high)) %>%
  mutate(p_value_FDR_adjust = p.adjust(p_value, method = "fdr")) %>%
  arrange(sel_gene,term)
data.frame(by_gene_all_cox_res_filt)

# 2) reg-AI

sel_genes <- AI_dNdES_test_res %>%
  filter(
    Model == "RNA_bb_model1",
    FDR_higher < 0.25
  ) %>%
  pull(Gene) %>% as.character()

by_gene_all_cox_models <- list()
by_gene_all_cox_res <- c()

for (sel_gene in sel_genes ) {
  print(paste0("-------------",sel_gene,"--------------"))

  # prepare a place to store the warning
  warn_obj <- NULL

  # run your model, catching any warning in warn_obj
  HR_cox_model_res   <- withCallingHandlers(
    HR_cox_model(sel_FDR    = 0.25,
                gene_set   = "all",
                AI_type    = "reg-AI",
                var_type   = "quant",
                gene_selected = sel_gene,
                cancer_type= "pancancer",
                AI_table   = TCGA_AI_table_final,
                TCGA_metadata = TCGA_metadata,
                covariate_df = covariates,
                heterogeneity = heterogeneity,
                CNA_signatures = CNA_signatures),
    warning = function(w) {
      warn_obj <<- w            # stash the warning condition
      invokeRestart("muffleWarning")  # prevent it from printing
    }
  )

  if ( !is.null(HR_cox_model_res$df_res) ) {
    # now `res` is your list and `warn_obj` is the warning condition
    if (is.null(warn_obj)) {
      warn_msg <- "none"
    } else {
      warn_msg <- conditionMessage(warn_obj)
      print(warn_msg)
    }
    HR_cox_model_res$df_res$warn_msg <- warn_msg
    HR_cox_model_res$df_res$sel_gene <- sel_gene
    by_gene_all_cox_models[[sel_gene]] <- HR_cox_model_res
    by_gene_all_cox_res <- rbind(by_gene_all_cox_res,HR_cox_model_res$df_res)
  }
}

by_gene_all_cox_res_filt <- by_gene_all_cox_res %>%
  arrange(desc(HR)) %>%
  filter(HR < 20) %>%
  filter(warn_msg == "none") %>%
  filter(CI_high < 20 & CI_low < 20) %>%
  filter(!is.infinite(CI_low) & !is.infinite(CI_high)) %>%
  mutate(p_value_FDR_adjust = p.adjust(p_value, method = "fdr")) %>%
  arrange(sel_gene,term)
data.frame(by_gene_all_cox_res_filt)

###################################################
################### FIGURE 4D #####################
###################################################

# A) KM curve for HNSC

# surv_df <- build_surv_df(
#   AI_dNdES_test_res      = AI_dNdES_test_res,
#   TCGA_AI_table_final    = TCGA_AI_table_final,
#   TCGA_metadata          = TCGA_metadata,
#   AI_type                = "mRNA-AI",
#   gene_selected         = NULL,
#   cancer_type            = "TCGA-HNSC", #LGG, LUAD, PAAD, STAD, HNSC, UCEC, COAD
#   sel_FDR                = 0.25,
#   gene_set               = "all",
#   var_type               = "qual",
#   pos_threshold          = 0.1,
#   neg_threshold          = -0.1,
#   n_pos_bins             = 3
# )

# # Fit and Plot
# km_fit <- survfit(Surv(time, status) ~ AI, data = surv_df)

# plot <- ggsurvplot(
#   km_fit, data = surv_df,
#   risk.table       = TRUE,
#   pval             = TRUE,
#   conf.int         = FALSE,
#   xlab             = "Time (years)",
#   ylab             = "Survival probability",
#   # legend.labs      = c("Not affected","Affected"),
#   legend.title     = "AI",
#   xlim = c(0,15),
#   surv.median.line = "v",
#   ggtheme          = theme_bw(base_size = 12)
#   )

# plot_to_save <- ggsave_workaround(plot)

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
# ggsave(final_figure_path, plot_to_save, width = 175, height = 180, units = "mm") 

# # Save
# # write.table(surv_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig4/Fig4C.txt", 
# #                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# # saveRDS(surv_df, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig4/Fig4C.RData")

# # Manuscript --> P value (log-rank)
# surv_pvalue(km_fit,data = surv_df)
# km_fit

# # Median survival time
# summary(km_fit)$table[1,"median"]
# mean(summary(km_fit)$table[-c(1:2),"median"],na.rm=TRUE)

# B) COX REGRESSION BY CANCER TYPE

# 1) mRNA-AI

by_cancer_all_cox_models <- list()
by_cancer_all_cox_res <- c()

for (cancer_type in unique(as.character(TCGA_AI_table_final$cancer)) ) {
  print(paste0("-------------",cancer_type,"--------------"))

  # prepare a place to store the warning
  warn_obj <- NULL

  # run your model, catching any warning in warn_obj
  HR_cox_model_res   <- withCallingHandlers(
    HR_cox_model(sel_FDR    = 0.25,
                gene_set   = "all",
                AI_type    = "mRNA-AI",
                var_type   = "qual",
                cancer_type= cancer_type,
                AI_table   = TCGA_AI_table_final,
                TCGA_metadata = TCGA_metadata,
                covariate_df = covariates,
                heterogeneity = heterogeneity,
                CNA_signatures = CNA_signatures),
    warning = function(w) {
      warn_obj <<- w            # stash the warning condition
      invokeRestart("muffleWarning")  # prevent it from printing
    }
  )

  if ( !is.null(HR_cox_model_res$df_res) ) {
    # now `res` is your list and `warn_obj` is the warning condition
    
    if (!is.null(warn_obj)) {
      warn_msg <- conditionMessage(warn_obj)
    } else {
      warn_msg <- "success"
    }
    HR_cox_model_res$df_res$warn_msg <- warn_msg
    HR_cox_model_res$df_res$cancer_type <- cancer_type
    by_cancer_all_cox_models[[cancer_type]] <- HR_cox_model_res
    by_cancer_all_cox_res <- rbind(by_cancer_all_cox_res,HR_cox_model_res$df_res)
  }
}

by_cancer_all_cox_res_filt <- by_cancer_all_cox_res %>%
  filter(
    across(c(HR, CI_low, CI_high), ~ . < 10),
    !is.infinite(CI_low),
    !is.infinite(CI_high),
    warn_msg != "Ran out of iterations and did not converge"
    # !cancer_type %in% c("TCGA-SARC","TCGA-UCS") # cancer types with < 50 samples
  ) %>%
  mutate( p_value_FDR_adjust = p.adjust(p_value, method = "fdr")) %>%
  arrange(HR)
data.frame(by_cancer_all_cox_res_filt)

by_cancer_all_cox_models_filt <- by_cancer_all_cox_models[names(by_cancer_all_cox_models) %in% by_cancer_all_cox_res_filt$cancer_type]
coxmodels <- lapply(by_cancer_all_cox_models_filt, function(x){x[["coxmodel"]]})

# Save
# write.table(coxmodels, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig4/Fig4D.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(coxmodels, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/Fig4/Fig4D.RData")

# 2) reg-AI

by_cancer_all_cox_models <- list()
by_cancer_all_cox_res <- c()

for (cancer_type in unique(as.character(TCGA_AI_table_final$cancer)) ) {
  print(paste0("-------------",cancer_type,"--------------"))

  # prepare a place to store the warning
  warn_obj <- NULL

  # run your model, catching any warning in warn_obj
  HR_cox_model_res   <- withCallingHandlers(
    HR_cox_model(sel_FDR    = 0.25,
                gene_set   = "all",
                AI_type    = "reg-AI",
                var_type   = "qual",
                cancer_type= cancer_type,
                AI_table   = TCGA_AI_table_final,
                TCGA_metadata = TCGA_metadata,
                covariate_df = covariates,
                heterogeneity = heterogeneity,
                CNA_signatures = CNA_signatures),
    warning = function(w) {
      warn_obj <<- w            # stash the warning condition
      invokeRestart("muffleWarning")  # prevent it from printing
    }
  )

  if ( !is.null(HR_cox_model_res$df_res) ) {
    # now `res` is your list and `warn_obj` is the warning condition
    
    if (!is.null(warn_obj)) {
      warn_msg <- conditionMessage(warn_obj)
    } else {
      warn_msg <- "success"
    }
    HR_cox_model_res$df_res$warn_msg <- warn_msg
    HR_cox_model_res$df_res$cancer_type <- cancer_type
    by_cancer_all_cox_models[[cancer_type]] <- HR_cox_model_res
    by_cancer_all_cox_res <- rbind(by_cancer_all_cox_res,HR_cox_model_res$df_res)
  }
}

by_cancer_all_cox_res_filt <- by_cancer_all_cox_res %>%
  filter(
    across(c(HR, CI_low, CI_high), ~ . < 10),
    !is.infinite(CI_low),
    !is.infinite(CI_high),
    warn_msg != "Ran out of iterations and did not converge"
  ) %>%
  mutate( p_value_FDR_adjust = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_value,cancer_type)
data.frame(by_cancer_all_cox_res_filt)

#################################################################################################
##################################### ASCA (ATAC-seq) ANALYSIS ##################################
#################################################################################################

#######################################################
################### SUPP. FIG S8A #####################
#######################################################

# ASCA data
ASCA_table_merged <- readRDS("/g/strcombio/fsupek_cancer1/gpalou/ASE_project/TCGA_ASCA_table.RData")

# Merge
cols <- c("chrom","pos","ref","alt","sample","Gene.Symbol","RNA_ref_ASCA","RNA_alt_ASCA")
TCGA_AI_table_final <- merge(TCGA_AI_table_final,ASCA_table_merged[,cols], all.x = TRUE)

# Beta-binomial model for reg-AI (ASCA)
RNA_ASCA_bb_model1 <- VGAM::vglm(
  cbind(RNA_alt_ASCA, RNA_ref_ASCA) ~ MUT_CNA_fraction + sqrt(gene_exp),
  family = VGAM::betabinomial,
  na.action = na.exclude,
  data = TCGA_AI_table_final
)

# Get expected proportions from the model
expected_prop <- fitted(RNA_ASCA_bb_model1)[, 1]
TCGA_AI_table_final$RNA_ASCA_expected_ratio <- expected_prop
# Observed proportions
TCGA_AI_table_final$RNA_VAF_ASCA <- TCGA_AI_table_final$RNA_alt_ASCA  / (TCGA_AI_table_final$RNA_alt_ASCA + 
        TCGA_AI_table_final$RNA_ref_ASCA)
TCGA_AI_table_final$RNA_VAF <- TCGA_AI_table_final$RNA_alt  / (TCGA_AI_table_final$RNA_alt + 
        TCGA_AI_table_final$RNA_ref)

# Get dispersion parameter from model
model_coefs <- coef(RNA_ASCA_bb_model1)
# Use plogis for dispersion as it seems to be the approach you're using
model_dispersion <- plogis(model_coefs["(Intercept):2"])

# Calculate p-values for model (both tails)
alt_counts <- TCGA_AI_table_final$RNA_alt_ASCA
total_counts <- TCGA_AI_table_final$RNA_alt_ASCA + TCGA_AI_table_final$RNA_ref_ASCA

# Lower tail: observed is lower than expected (favors reference)
TCGA_AI_table_final$RNA_ASCA_p_lower_bb_model1 <- VGAM::pbetabinom(
  alt_counts, size = total_counts, prob = expected_prop, rho = model_dispersion
)

# Upper tail: observed is higher than expected (favors alternate)
TCGA_AI_table_final$RNA_ASCA_p_upper_bb_model1 <- 1 - VGAM::pbetabinom(
  alt_counts - 1, size = total_counts, prob = expected_prop, rho = model_dispersion
)

p_value <- 0.05
TCGA_AI_table_final <- TCGA_AI_table_final %>%
  dplyr::mutate(
    pos_reg_AI_ASCA = ( RNA_ASCA_p_upper_bb_model1 < p_value &
            DNA_p_upper_null_model >= p_value),
    neg_reg_AI_ASCA = ( RNA_ASCA_p_lower_bb_model1 < p_value &
            DNA_p_lower_null_model >= p_value),
    no_reg_AI_ASCA = !(pos_reg_AI_ASCA | neg_reg_AI_ASCA),
    all_reg_AI_ASCA = (pos_reg_AI_ASCA | neg_reg_AI_ASCA)
  ) %>%
  # define the three reg-AI classes
  mutate(
    reg_AI_ASCA_type = case_when(
      pos_reg_AI_ASCA ~ "Pos reg-AI ASCA",
      neg_reg_AI_ASCA ~ "Neg reg-AI ASCA",
      TRUE       ~ "No reg-AI ASCA"
    ),
    # make sure it's a factor in the order you want
    reg_AI_ASCA_type = factor(
      reg_AI_ASCA_type,
      levels = c( "Pos reg-AI ASCA", "No reg-AI ASCA", "Neg reg-AI ASCA"),
      labels = c("Pos reg-AI ASCA", "No reg-AI ASCA", "Neg reg-AI ASCA")
    )
  )
table(TCGA_AI_table_final$reg_AI_ASCA_type)

# Coverage threshold
coverage_thres <- 20
TCGA_AI_table_final_filt <- TCGA_AI_table_final[(TCGA_AI_table_final$RNA_ref_ASCA + TCGA_AI_table_final$RNA_alt_ASCA > coverage_thres)
             & (TCGA_AI_table_final$RNA_alt + TCGA_AI_table_final$RNA_ref > coverage_thres),]
table(TCGA_AI_table_final_filt$reg_AI_ASCA_type)
table(TCGA_AI_table_final_filt$reg_AI_ASCA_type,TCGA_AI_table_final_filt$reg_AI_type)

# (1) Define the new 3‐level variable:
TCGA_AI_table_final_filt$combined_type <- with(
  TCGA_AI_table_final_filt,
  ifelse(
    reg_AI_ASCA_type == "Pos reg-AI ASCA" & reg_AI_type == "Pos reg-AI",
    "Pos reg-AI RNA and ASCA",
    ifelse(
      reg_AI_ASCA_type == "Neg reg-AI ASCA" & reg_AI_type == "Neg reg-AI",
      "Neg reg-AI RNA and ASCA",
      "other"
    )
  )
)

# (2) Turn it into a factor (optional, but useful for ordering/plotting):
TCGA_AI_table_final_filt$combined_type <- factor(
  TCGA_AI_table_final_filt$combined_type,
  levels = c("Pos reg-AI RNA and ASCA", "Neg reg-AI RNA and ASCA", "other")
)
TCGA_AI_table_final_filt <- TCGA_AI_table_final_filt %>%
  filter(!is.na(TCGA_AI_table_final_filt$RNA_ref_ASCA)) %>%
  # filter(!gene %in% c("essential","random")) %>%
  # filter(gene %in% c("random")) %>%
  mutate(
    SNV_varity = case_when(
      SNV_varity == "effectively_syn"     ~ "Eff. synonymous",
      SNV_varity == "nonsynonymous SNV"   ~ "Missense",
      SNV_varity == "stopgain"            ~ "Nonsense",
      # for everything else, keep the original string:
      TRUE                                 ~ as.character(SNV_varity)
    )
  )
        
dim(TCGA_AI_table_final_filt)
table(TCGA_AI_table_final_filt$SNV_varity)
# (3) Check the counts:
table(TCGA_AI_table_final_filt$combined_type)

get_top_pct <- function(x, pct) {
  stopifnot(is.numeric(x), pct >= 0, pct <= 100)
  # 1 - pct/100 is the quantile for the bottom (1 - pct) fraction
  q <- quantile(x, probs = 1 - pct/100, na.rm = TRUE)
  unname(q)
}
# # A) Top 5%
top_thres <- 5
sei_thres <- get_top_pct(TCGA_AI_table_final$Sei_MaxDelta, top_thres)
sei_thres

TCGA_AI_table_final_filt2 <- TCGA_AI_table_final_filt %>%
        filter(Sei_MaxDelta <= sei_thres)

plot <- ggplot(
    data = TCGA_AI_table_final_filt, # or if Sei-removal use this --> TCGA_AI_table_final_filt2
    mapping = aes(
      # x = RNA_VAF_ASCA,
      # y = RNA_VAF,
      x = RNA_VAF_ASCA - RNA_ASCA_expected_ratio,
      y = RNA_VAF - RNA_expected_ratio_bb_model1,
      color = combined_type,
      shape = combined_type
    )
    ) +
    # Scatter points, colored by reg_AI_type
    geom_point(size = 2, alpha = 0.8) +
    
    # Single regression line across all points (ignore color groups)
    geom_smooth(
      aes(group = 1),
      method = "lm",
      formula = y ~ x,
      color = "black",
      se = FALSE
    ) +
    
    # Pearson r and p-value (computed on the entire dataset, not by color)
    stat_cor(
      aes(
        # x = RNA_VAF_ASCA,
        # y = RNA_VAF,
        group = 1,
        x = RNA_VAF_ASCA - RNA_ASCA_expected_ratio,
        y = RNA_VAF - RNA_expected_ratio_bb_model1,
        label = paste(
          after_stat(r.label),
          after_stat(p.label),
          sep = "~`,`~"
        )
      ),
      method = "pearson",
      label.x.npc = "left",   # place near left
      label.y.npc = "top",    # place near top
      size = 4,
      color = "black",
    ) +
    facet_wrap(SNV_varity ~ . ) +
    
    # (Optional) If you wanted separate panels per reg_AI_type, uncomment:
    # facet_wrap(~ reg_AI_type) +
    
    # Final labels and theme
    labs(
      x     = expression(italic(ASCA)~VAF[obs] - italic(ASCA)~VAF[pred]),
      y     = expression(italic(RNA)~VAF[obs] - italic(RNA)~VAF[pred]),
      color = "reg-AI",
      title = "",
      subtitle = ""
    ) +
    # — now override colors and shapes by *the exact* level names —
    scale_color_manual(
      name   = "sig AI",
      values = c(
        "Pos reg-AI RNA and ASCA" = "#4DBBD5FF",
        "other"                   = "#00A087FF",
        "Neg reg-AI RNA and ASCA" = "#3C5488FF"
      )
    ) +
    scale_shape_manual(
      name   = "sig AI",
      values = c(
        "Neg reg-AI RNA and ASCA" = 25,
        "other"                   = 4,
        "Pos reg-AI RNA and ASCA" = 19
      )
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12)
    ) +
    guides(color = "none")

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, plot, width = 325, height = 150, units = "mm") 

write.table(TCGA_AI_table_final_filt, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig8/panel_A.txt", 
                sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(TCGA_AI_table_final_filt, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig8/panel_A.RData")

########################################################################################################
##################################### METHYLATION ANALYSIS #############################################
########################################################################################################

#######################################################
################### SUPP. FIG S8B #####################
#######################################################

plot_df <- TCGA_AI_table_final
plot_df$reg_AI_type <- gsub(" reg-AI", "",plot_df$reg_AI_type)
gene_cognate_types <- names(table(TCGA_AI_table_final$gene_cognate))
gene_cognate_types <- gene_cognate_types[!gene_cognate_types %in% c("both","both_cognate")]
SNV_types <- names(table(TCGA_AI_table_final$SNV_varity))
SNV_types <- c(SNV_types,"Low Impact Missense","Synonymous","All SNVs")
SNV_types <- SNV_types[SNV_types != "stopgain"]
reg_AI_cat <- c("Pos","Neg")
print(gene_cognate_types)
print(SNV_types)
print(reg_AI_cat)
# Methylation_index <- grep("category",colnames(TCGA_AI_table_final))
variables_columns <- c("M_values_CpG_TSS1500_category","M_values_CpG_TSS200_category","M_values_CpG_TSS1500_TSS200_category",
      "M_values_CpG_TSS200_1500_body_first_2kb_category","M_values_CpG_body_exclude_first_2kb_category",
      "M_values_CpG_body_first_2kb_category")

all_df_fisher_test_res <- c()

for (SNV_type_char in SNV_types) {
  print(SNV_type_char)
  df_test1 <- if (SNV_type_char %in% c("effectively_syn","nonsynonymous SNV")) {
    filter(plot_df, SNV_varity == SNV_type_char)
  } else if (SNV_type_char %in% c("Low Impact Missense","Synonymous")) {
    filter(plot_df, SNV_varity_new == SNV_type_char)
  } else if (SNV_type_char == "All SNVs") {
    plot_df
  }

    for (gene_type_char in gene_cognate_types) {
        print(gene_type_char)
        # Gene type
        df_test2 <- df_test1 %>%
            filter(gene_cognate %in% c(gene_type_char))
        
        for (reg_AI_type_char in reg_AI_cat) {
            print(reg_AI_type_char)
            # reg-AI type
            if (reg_AI_type_char == "Pos") {
                df_test3 <- df_test2 %>%
                    filter(reg_AI_type %in% c("Pos","No"))
            } else if (reg_AI_type_char == "Neg") {
                df_test3 <- df_test2 %>%
                    filter(reg_AI_type %in% c("Neg","No"))
            }

            # print(tail(sort(table(df_test3$Gene.Symbol)),3))
            if (nrow(df_test3) == 0) {next}
            # next

            for (col_var in variables_columns) {
                # print(paste0("<--------",col_var,"-------->"))
                df_fisher_test_res <- data.frame(OR = NA, p_value = NA, CI_95_low = NA, CI_95_high = NA,
                                            CI_90_low = NA, CI_90_high = NA, CI_80_low = NA, CI_80_high = NA,
                                                SNV_type = SNV_type_char, gene_type = gene_type_char,
                                                reg_AI_cat = reg_AI_type_char, col_var = col_var)
                # Column variable type
                if (col_var %in% c("M_values_CpG_TSS1500_category","M_values_CpG_TSS200_category",
                          "M_values_CpG_TSS1500_TSS200_category","M_values_CpG_TSS200_1500_body_first_2kb_category",
                          "M_values_CpG_body_first_2kb_category","M_values_CpG_body_exclude_first_2kb_category")) {
                    # Methylation
                    df_test3$variant_class <- NA
                    df_test3[which(df_test3[,col_var] == "Partially methylated"),"variant_class"] <- paste0("High")    
                    df_test3[which(df_test3[,col_var] == "Unmethylated"),"variant_class"] <- paste0("Low")    
                } 
                # else if (col_var %in% c("M_values_CpG_body_exclude_first_2kb_category")) {
                #     # Methylation
                #     df_test3$variant_class <- NA
                #     df_test3[which(df_test3[,col_var] == "Partially methylated"),"variant_class"] <- paste0("High")    
                #     df_test3[which(df_test3[,col_var] == "Methylated"),"variant_class"] <- paste0("Low")   
                # }
                df_test3$reg_AI_type <- factor(df_test3$reg_AI_type, levels = c("No",reg_AI_type_char))
                df_test3$variant_class <- factor(df_test3$variant_class, levels = c("Low","High"))
                # Fisher Test comparing proportions (pos/neg vs no cis-ASE) in X gene type
                # Test is: ([X] variant class: High / Low) / (cis-ASE: pos/neg / no)
                df_test_table <- table(df_test3$variant_class, df_test3$reg_AI_type)
                # Use pseudocount of 1 (NO)
                df_test_table <- (df_test_table)
                fisher_test_res <- fisher.test(df_test_table , conf.int = TRUE, conf.level = 0.95,
                                            alternative = "two.sided")
                df_fisher_test_res$OR <- as.numeric(fisher_test_res$estimate)
                df_fisher_test_res$p_value <- as.numeric(fisher_test_res$p.value)
                # CI
                df_fisher_test_res$CI_95_low <- as.numeric(fisher_test_res$conf.int)[1]
                df_fisher_test_res$CI_95_high <- as.numeric(fisher_test_res$conf.int)[2]
                fisher_test_res <- fisher.test(df_test_table , conf.int = TRUE, conf.level = 0.90,alternative = "two.sided")
                df_fisher_test_res$CI_90_low <- as.numeric(fisher_test_res$conf.int)[1]
                df_fisher_test_res$CI_90_high <- as.numeric(fisher_test_res$conf.int)[2]
                fisher_test_res <- fisher.test(df_test_table , conf.int = TRUE, conf.level = 0.80,alternative = "two.sided")
                df_fisher_test_res$CI_80_low <- as.numeric(fisher_test_res$conf.int)[1]
                df_fisher_test_res$CI_80_high <- as.numeric(fisher_test_res$conf.int)[2]
                # Sample size
                df_fisher_test_res$n_yes_reg_AI_High_variable <- df_test_table["High",reg_AI_type_char]
                df_fisher_test_res$n_yes_reg_AI_Low_variable <- df_test_table["Low",reg_AI_type_char]
                df_fisher_test_res$n_no_reg_AI_High_variable <- df_test_table["High","No"]
                df_fisher_test_res$n_no_reg_AI_Low_variable <- df_test_table["Low","No"]
                # Save
                all_df_fisher_test_res <- rbind(all_df_fisher_test_res, df_fisher_test_res)
            }
        }
    }
}

head(all_df_fisher_test_res)

change_names <- function(df) {
  df %>%
    # 1) rename columns
    rename(
      SNV_varity   = SNV_type,
      gene_cognate = gene_type
    ) %>%
    # 2) recode everything in one go
    mutate(
      # gene_cognate relabel & factor‐level order
      gene_cognate = recode(
        gene_cognate,
        random       = "Passenger",
        essential    = "Essential",
        og_cognate   = "Cognate OG",
        oncogene     = "Noncognate OG",
        tsg          = "Noncognate TSG",
        tsg_cognate  = "Cognate TSG",
        .default     = as.character(gene_cognate)
      ),
      gene_cognate = factor(
        gene_cognate,
        levels = c("Passenger","Essential",
                   "Noncognate OG","Cognate OG",
                   "Noncognate TSG","Cognate TSG",
                   "TP53")
      ),

      # SNV_varity relabel
      SNV_varity = recode(
        SNV_varity,
        effectively_syn       = "ES",
        `nonsynonymous SNV`   = "High aa-impact Miss.",
        `Low Impact Missense` = "Low aa-impact Miss.",
        # stopgain              = "Nonsense",
        .default              = as.character(SNV_varity)
      ),

      # col_var relabel
      col_var = recode(
        col_var,
        M_values_CpG_TSS1500_category                           = "TSS1500",
        M_values_CpG_TSS200_category                            = "TSS200",
        M_values_CpG_TSS1500_TSS200_category                    = "TSS1500-200",
        M_values_CpG_TSS200_1500_body_first_2kb_category        = "TSS1500-200_body_2kb",
        M_values_CpG_body_first_2kb_category                    = "Body_2kb",
        M_values_CpG_body_exclude_first_2kb_category             = "Body_no_2kb",
        .default         = as.character(col_var)
      ),
      col_var = factor(
        col_var,
        levels = c("TSS1500","TSS200","TSS1500-200","TSS1500-200_body_2kb","Body_2kb","Body_no_2kb")
      ),
      # reg_AI_cat relabel & factor-levels
      reg_AI_cat = recode(
        reg_AI_cat,
        Pos = "pos reg-AI",
        Neg = "neg reg-AI",
        .default = as.character(reg_AI_cat)
      ),
      reg_AI_cat = factor(
        reg_AI_cat,
        levels = c("pos reg-AI","neg reg-AI")
      )
    )
}

# Change name of variables
all_df_fisher_test_res_filt <- change_names(all_df_fisher_test_res)

# Remove Inf values
all_df_fisher_test_res_filt <- all_df_fisher_test_res_filt %>% 
  filter(OR != Inf & OR != -Inf & OR != 0)

# ORs FDR correction
all_df_fisher_test_res_filt <- all_df_fisher_test_res_filt %>%
  group_by(reg_AI_cat,SNV_varity) %>%
  mutate(p_value_FDR_adjusted = p.adjust(p_value, method = "fdr"))

# Add labels with conditional coloring
all_df_fisher_test_res_filt <- all_df_fisher_test_res_filt %>%
    mutate(
        # < 0.15 --> "*"; < 0.05 --> "**"; < 0.01 --> "***"
        label = case_when(
            p_value_FDR_adjusted > 0.15 ~ "",
            p_value_FDR_adjusted > 0.05 ~ "*",
            p_value_FDR_adjusted > 0.01 ~ "**",
            !is.na(p_value_FDR_adjusted) ~ "***",
            TRUE ~ NA_character_
        )
    )

plot_df <- all_df_fisher_test_res_filt %>%
        filter(!gene_cognate %in% c("Essential"))

plot <- ggplot(plot_df, aes(x = col_var, y = gene_cognate)) + 
    geom_tile(aes(fill = log2(OR))) +
    scale_fill_gradient2(low = "#0D5EAF", mid = "white", high = "#AF540D", midpoint = 0,
                        #  limits = c(-2.2, 3.65), breaks = seq(-2.2, 3.65, by = 1), 
                        name = expression(log[2](`OR`))) +
    geom_text(aes(label = label), nudge_y = 0, size = 3) +
    scale_color_identity() +  # Ensures colors are taken directly from 'label_color'
    facet_grid(reg_AI_cat ~ SNV_varity) +
    labs(x = "", y = "") +
    theme_classic()+
    theme(
        legend.position = "right",
        axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(size = 9),
        plot.margin = unit(c(0, 1, 0, 1), "cm")
    )

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, plot, width = 325, height = 150, units = "mm") 

# Enrichment analysis --> only  TSS1500-200_body_2kb

calculate_bd_p <- function(df) {
  # must have exactly Passenger + one other
  if(nrow(df) != 2) return(NA_real_)
  other <- setdiff(df$gene_cognate, "Passenger") 
  # reorder Passenger first
  df <- df %>% arrange(match(gene_cognate, c("Passenger", other)))
  # pull out & coerce to numeric
  h_yes <- as.numeric(df$n_yes_reg_AI_High_variable[1:2])
  h_no  <- as.numeric(df$n_no_reg_AI_High_variable[1:2])
  l_yes <- as.numeric(df$n_yes_reg_AI_Low_variable[1:2])
  l_no  <- as.numeric(df$n_no_reg_AI_Low_variable[1:2])
  # bail on any missing
  if(any(is.na(c(h_yes, h_no, l_yes, l_no)))) return(NA_real_)
  # build 2×2×2 array as numeric
  vals <- c(
    # slice 1 = Passenger
    h_yes[1], h_no[1], l_yes[1], l_no[1],
    # slice 2 = other class
    h_yes[2], h_no[2], l_yes[2], l_no[2]
  )
  tbl3 <- array(
    vals,
    dim   = c(2, 2, 2),
    dimnames = list(
      SplicingImpact   = c("High","Low"),
      AlleleImbalance  = c("Yes","No"),
      GeneSet          = c("Passenger", other)
    )
  )
  # run test quietly, catch errors
  pval <- tryCatch({
    suppressWarnings(
      BreslowDayTest(tbl3, correct = FALSE)$p.value
    )
  }, error = function(e) NA_real_)
  pval
}

all_df_bd_res <- map_dfr(
  c("Cognate OG","Noncognate OG","Cognate TSG","Noncognate TSG","Essential","TP53"),
  function(gene_class){
    all_df_fisher_test_res_filt %>%
      filter(gene_cognate %in% c(gene_class, "Passenger")) %>%
      group_by(col_var, SNV_varity, reg_AI_cat) %>%
      nest() %>%
      mutate(
        breslow_day_p = map_dbl(data, calculate_bd_p),
        compare_to    = gene_class
      ) %>%
      select(-data)
  }
)

df_plot <- all_df_fisher_test_res_filt %>%
  # keep exactly the two MaxDelta tests you want
  filter(
    (col_var == "TSS1500-200_body_2kb"),
    (reg_AI_cat == "pos reg-AI")
  ) %>%
  # bring in your Breslow‐Day p
  left_join(
    all_df_bd_res %>%
      select(SNV_varity, reg_AI_cat, col_var,
             compare_to, breslow_day_p),
    by = c(
      "SNV_varity", "reg_AI_cat", "col_var",
      "gene_cognate" = "compare_to"
    )
  ) %>%
  # for each gene+facet, compute Fisher’s combined p
  group_by(gene_cognate, SNV_varity, reg_AI_cat) %>%
  mutate(
    # only two rows per group, so this sums exactly those two p’s
    fisher_X2 = -2 * sum(log(breslow_day_p)),
    fisher_df = 2 * n(),  
    fisher_p  = pchisq(fisher_X2, df = fisher_df, lower.tail = FALSE)
  ) %>%
  ungroup()# %>%
  # drop the fisher_p from the PG rows so only SP shows it
  # mutate(
  #   fisher_p = ifelse(col_var == "SP_MaxDelta", fisher_p, NA_real_)
  # )

df_plot$gene_cognate <- factor(df_plot$gene_cognate, 
      levels = levels(all_df_fisher_test_res_filt$gene_cognate))

plot_levels <- levels(df_plot$gene_cognate)

# 2) build annot_df with a dynamic y based on rank
annot_df <- df_plot %>%
  filter(gene_cognate != "Passenger") %>%
  distinct(reg_AI_cat, SNV_varity, gene_cognate, fisher_p) %>%
  mutate(
    gene_cognate = factor(gene_cognate, levels = plot_levels),
    xstart       = as.numeric(factor("Passenger", levels = plot_levels)),
    xend         = as.numeric(gene_cognate)
  ) %>%
  group_by(reg_AI_cat, SNV_varity) %>%
  arrange(xend, .by_group = TRUE) %>%
  mutate(
    rank = row_number(),
    y0   = 2 + (rank - 1) * 0.25,
    # only label the significant ones
    label = ifelse(!is.na(fisher_p) & fisher_p < 0.05,
                   paste0("p=", format(fisher_p, scientific = TRUE, digits = 2)),
                   NA_character_),
    # only keep y where you have a label
    y = ifelse(!is.na(label), y0, NA_real_)
  ) %>%
  ungroup()

# 3) re‐plot with these staggered y’s
plot <- ggplot(df_plot,
            aes(x     = gene_cognate,
                y     = log2(OR),
                color = gene_cognate,
                group = gene_cognate)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = log2(CI_95_low),
                    ymax = log2(CI_95_high)),
                width    = 0.3,
                position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(reg_AI_cat ~ SNV_varity, scales = "fixed") +
  coord_cartesian(ylim = c(-2, 4)) +   # extend top margin
  # coord_cartesian(ylim = c(-3, max(annot_df$y) + 0.5)) +   # extend top margin
  scale_color_manual(values = c(
    "Passenger"   = "#83c4be",
    "Cognate OG"  = "#37d576",
    "Cognate TSG" = "#7083e5",
    "TP53"        = "#7cd8ef",
    "Essential"   = "#db4538"
  )) +
  # scale_shape_manual(values = c(
  #   "PG_MaxDelta" = 16,
  #   "SP_MaxDelta" = 17
  # )) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    strip.text      = element_text(size = 7.5),
    plot.title      = element_text(size = 10)
  ) +
  labs(
    x     = NULL,
    y     = expression(log[2](OR)),
    shape = "",
    title = ""
  ) +
  geom_segment(data        = annot_df,
               aes(x    = xstart,
                   xend = xend,
                   y    = y-0.2,
                   yend = y-0.2),
               inherit.aes = FALSE) +
  geom_text(data        = annot_df,
            aes(x     = (xstart + xend)/2,
                y     = y + 0.1,
                label = label),
            inherit.aes = FALSE,
            size = 3)

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, plot, width = 275, height = 100, units = "mm") 

plot_list <- list(df_plot = df_plot, annot_df = annot_df)

# Save
# write.table(TCGA_AI_table_final_filt, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig8/panel_B.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(plot_list, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig8/panel_B.RData")

# For the Manuscript
df_plot %>%
  filter(gene_cognate == "Cognate TSG") %>%
  arrange(OR) %>%
  data.frame()

#######################################################
################### SUPP. FIG S9 ######################
#######################################################

##### Reviewer #2 Question 4: Subclonal variants #####

# Number of variants for cancer genes
df <- TCGA_AI_table_final %>%
  # filter(SNV_varity != "stopgain") %>%
  filter(!gene %in% c("essential","random"))
dim(df)
638/nrow(df)

df <- TCGA_AI_table_final %>%
  filter(SNV_varity != "stopgain") %>%
  filter(!gene %in% c("essential","random")) %>%
  filter((( DNA_alt/ (DNA_alt + DNA_ref) ) / (purity)) < 0.2) %>%
  filter(!MUT_CNA2 %in% c(0)) %>%
  filter(all_mRNA_AI)
dim(df)
# 1) % pie chart plot
# Calculate proportions
reg_AI_props <- prop.table(table(df$reg_AI_type)) * 100
CNA_AI_props <- prop.table(table(df$CNA_AI_type)) * 100
table(df$reg_AI_type)
table(df$CNA_AI_type) # 458 Neg CNA-AI
num_neg_CNA_AI_subclonal <- as.numeric(table(df$CNA_AI_type)["Neg CNA-AI"])
num_neg_reg_AI_subclonal <- as.numeric(table(df$reg_AI_type)["Neg reg-AI"])
reg_AI_props
CNA_AI_props
# Regulatory AI pie chart
reg_data <- data.frame(
  Category = names(reg_AI_props),
  Proportion = as.numeric(reg_AI_props)
)

p1 <- ggplot(reg_data, aes(x = "", y = Proportion, fill = Category)) +
  geom_col(width = 1, color = "white", size = 0.8) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c(
    "Pos reg-AI" = "#FF6B6B", 
    "No reg-AI" = "#4ECDC4", 
    "Neg reg-AI" = "#45B7D1"
  )) +
  labs(title = "reg-AI distribution",
       fill = "Reg-AI") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    plot.margin = margin(0, 0, 0, 0, "cm")  # Remove plot margins
  ) +
  geom_text(aes(label = paste0(round(Proportion, 1), "%")), 
            position = position_stack(vjust = 0.5),
            color = "black", fontface = "bold", size = 4)

# CNA AI pie chart
cna_data <- data.frame(
  Category = names(CNA_AI_props),
  Proportion = as.numeric(CNA_AI_props)
)

p2 <- ggplot(cna_data, aes(x = "", y = Proportion, fill = Category)) +
  geom_col(width = 1, color = "white", size = 0.8) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c(
    "Pos CNA-AI" = "#FF6B6B", 
    "No CNA-AI" = "#4ECDC4", 
    "Neg CNA-AI" = "#45B7D1"
  )) +
  labs(title = "CNA-AI distribution",
       fill = "CNA-AI") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    plot.margin = margin(0, 0, 0, 0, "cm")  # Remove plot margins
  ) +
  geom_text(aes(label = paste0(round(Proportion, 1), "%")), 
            position = position_stack(vjust = 0.5),
            color = "black", fontface = "bold", size = 4)

plot <- p1 + p2 + 
  plot_annotation(title = "reg-AI vs CNA-AI % - Subclonal variants",
                  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))

final_figure_path <- "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/test.png"
ggsave(final_figure_path, plot, width = 200, height = 200, units = "mm") 
            
# 2) % of significant CNA-AI due to subclonal variants
# df <- TCGA_AI_table_final %>%
#   filter(SNV_varity != "stopgain") %>%
#   filter(!gene %in% c("essential","random")) %>%  
#   filter(!MUT_CNA2 %in% c(0)) %>%
#   # filter(MUT_CNA2 == 1 & WT_CNA == 1) %>%
#   filter(all_mRNA_AI)
# dim(df)
# table(df$CNA_AI_type)
# table(df$reg_AI_type)
# df$subclonal <- ( df$DNA_alt/ (df$DNA_alt + df$DNA_ref) )/ (df$purity) < 0.2
# table(df$subclonal,df$CNA_AI_type)
# table(df$subclonal,df$reg_AI_type)

# num_neg_CNA_AI_subclonal/table(df$CNA_AI_type)["Neg CNA-AI"] * 100
# num_neg_CNA_AI_subclonal/sum(table(df$CNA_AI_type)[c("Pos CNA-AI","Neg CNA-AI")]) * 100

plot_list <- list(plot_df1 = reg_data, plot_df2 = cna_data)

# Save
# write.table(plot_list, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig9/panel_A.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
saveRDS(plot_list, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/SuppFig/SuppFig9/panel_A.RData")

###################################################################
################### CORREIA et. al. ANALYSIS ######################
###################################################################


# % OF AI VARIABILITY COMPARING TO Correia et. al.

# 1) Our original models (correcting by purity)

# df_processed <- TCGA_AI_table_final %>%
#   mutate(cancer = gsub("TCGA-", "", cancer))
# # Calculate pan-cancer proportions
# pancan_summary_df <- df_processed %>%
#   filter(all_mRNA_AI) %>%
#   mutate(
#     AI_source = case_when( 
#       all_reg_AI & !all_CNA_AI  ~ "reg-AI",
#       !all_reg_AI & all_CNA_AI  ~ "CNA-AI"
#     )
#   ) %>%
#   filter(!is.na(AI_source)) %>%
#   group_by(AI_source) %>%
#   summarise( n = n(), .groups = "drop" ) %>%
#   mutate(
#     total_exclusive_snvs_pancan = sum(n),
#     percent = (n / total_exclusive_snvs_pancan) * 100
#   )

# # a) Our analysis --> original

# total_reg_AI <- paste0(round(pancan_summary_df[2,"percent"],1),"%")
# total_CNA_AI <- paste0(round(pancan_summary_df[1,"percent"],1),"%")

# df_a <- data.frame(var_explained = c(total_reg_AI,total_CNA_AI), 
#                     ASE_type = c("reg-AI","CNA-ASE"), 
#                     paper = c("This work"))

# # b) Our analysis --> without correcting by purity + only Breast + PIK3CA
# # we skip --> nonsynonymous + stringent filterings
# pancan_summary_df <- df_processed %>%
#   filter(Gene.Symbol == "PIK3CA") %>%
#   filter(SNV_varity_new %in% c("High Impact Missense","Low Impact Missense")) %>%
#   filter(cancer == "BRCA") %>%
#   # filter(SNV_varity == "nonsynonymous SNV") %>%
#   # filter((DNA_alt + DNA_ref) >= 30) %>%
#   # filter((RNA_alt + RNA_ref) >= 30) %>%
#   filter(all_mRNA_AI) %>%
#   mutate(
#     AI_source = case_when( 
#       all_reg_AI & !all_CNA_AI  ~ "reg-AI",
#       !all_reg_AI & all_CNA_AI  ~ "CNA-AI"
#     )
#   ) %>%
#   filter(!is.na(AI_source)) %>%
#   group_by(AI_source) %>%
#   summarise( n = n(), .groups = "drop" ) %>%
#   mutate(
#     total_exclusive_snvs_pancan = sum(n),
#     percent = (n / total_exclusive_snvs_pancan) * 100
#   )
# pancan_summary_df
# # Save
# total_reg_AI <- paste0(round(pancan_summary_df[2,"percent"],1),"%")
# total_CNA_AI <- paste0(round(pancan_summary_df[1,"percent"],1),"%")

# df_f <- data.frame(var_explained = c(total_reg_AI,total_CNA_AI), 
#                     ASE_type = c("reg-AI","CNA-AI"), 
#                     paper = c("This work (5): (4) + only PIK3CA"))
# # g) Paper PIK3CA
# # PIK3CA_paper_METABRIC_breast_cis_ASE <- 0.206
# # PIK3CA_paper_METABRIC_breast_CNA_ASE <- 1 - PIK3CA_paper_METABRIC_breast_cis_ASE
# # PIK3CA_paper_TCGA_breast_cis_ASE <- 0.144
# # PIK3CA_paper_TCGA_breast_CNA_ASE <- 1 - PIK3CA_paper_TCGA_breast_cis_ASE
# PIK3CA_paper_cis_ASE <- 0.16
# PIK3CA_paper_CNA_ASE <- 1 - PIK3CA_paper_cis_ASE
# # Save
# df_g <- data.frame(var_explained = c(PIK3CA_paper_cis_ASE,PIK3CA_paper_CNA_ASE), 
#                     ASE_type = c("cis_ASE","CNA_ASE"), 
#                     paper = c("Correia et al. paper"))

# # FINAL PLOT
# # Merge
# final_df <- rbind(df_a,df_b,df_c,df_d,df_e,df_f,df_g)
# final_df <- ddply(final_df, "paper",
#                    transform, label_ypos = cumsum(var_explained))
# final_df$paper_category <- "Correia et. al 2022"

# plots <- list()
# i <- 1
# for (paper in c("Correia et. al 2022")) {

#     final_df_tmp <- final_df[final_df$paper_category == paper,]
#     final_df_tmp$paper <- gsub("ASE \\+ purity","ASE + purity\n",final_df_tmp$paper)
#     final_df_tmp$paper <- gsub("BRCA","BRCA\n",final_df_tmp$paper)
#     final_df_tmp$paper <- factor(final_df_tmp$paper, levels = unique(as.character(final_df_tmp$paper)))

#     plot <- ggplot(data= final_df_tmp, aes(x=paper, y=var_explained, fill=ASE_type)) +
#         # geom_bar(stat="identity")+
#         geom_col() + 
#         geom_text(aes(label = scales::percent(var_explained, accuracy = 0.01)), position = position_stack(.5), size = 4) +
#         coord_flip() + 
#         labs(fill = "", title = "Explained cis vs CNA ASE variability") +
#         # geom_text(aes(y=round(var_explained,2), label=round(var_explained,2)), vjust= -0.5, 
#         #         color="white", size=3.5)+
#         scale_fill_brewer(palette="Accent")+
#         facet_wrap(paper_category ~ .) +
#         # scale_fill_manual(values=c('#999999','#E69F00')) +
#         theme_classic() +
#         theme(plot.title = element_text(hjust = 0.5, size = 18),
#             axis.text.y = element_text(size = 10),
#             legend.position = "top")
#     plots[[i]] <- plot
#     i <- i + 1
# }

# final_plot <- plot_grid(plots[[1]],plots[[2]], nrow = 2, align = "v")

# final_figure_path <- "/g/strcombio/fsupek_home/gpalou/analysis_results/ASE_project/variance_explained/comparison_between_papers.png"
# ggsave(final_figure_path, final_plot, width = 250, height = 300, units = "mm") 

# # Save
# write.table(final_df, file = "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/SuppFig1/SuppFig1A.txt", 
#                 sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
# saveRDS(final_df, "/g/strcombio/fsupek_home/gpalou/Manuscript/ASE/figures/SuppFig1/SuppFig1A.RData")

# 2) Updated models NOT correcting by purity

# A) Create bb-models

# Dataframe
df <- TCGA_AI_table_final %>%
  filter(cancer == "TCGA-BRCA")
dim(df)
# Compute Estimated Mutant Copy Number (MUT_CNA) ---- Original formula
# df$estimate_MUT_CNA <- (df$DNA_ASE / df$purity) * 
#                     ((df$purity * (df$Copy_Number)) + 
#                     (2 * (1 - df$purity)))
# Compute Estimated Mutant Copy Number (MUT_CNA) ---- Not correcting by purity formula
df$estimate_MUT_CNA <- (df$DNA_ASE) * (df$Copy_Number)

# Compute the fraction of mutant copy number relative to total CN, capped at 1
df$MUT_CNA_fraction <- pmin(df$estimate_MUT_CNA / (df$Copy_Number), 1)
df$DNA_alt <- as.integer(df$DNA_alt)
df$DNA_ref <- as.integer(df$DNA_ref)

RNA_bb_model1 <- VGAM::vglm(
  cbind(RNA_alt, RNA_ref) ~ MUT_CNA_fraction + sqrt(gene_exp),
  family = VGAM::betabinomial,
  na.action = na.exclude,
  data = df
)

RNA_null_model <- VGAM::vglm(
  cbind(RNA_alt, RNA_ref) ~ sqrt(gene_exp),
  family = VGAM::betabinomial,
  na.action = na.exclude,
  data = df
)

DNA_bb_model1 <- VGAM::vglm(
  cbind(DNA_alt, DNA_ref) ~ MUT_CNA_fraction,
  family = VGAM::betabinomial,
  na.action = na.exclude,
  data = df
)

DNA_null_model <- VGAM::vglm(
  cbind(DNA_alt, DNA_ref) ~ 1,
  family = VGAM::betabinomial,
  na.action = na.exclude,
  data = df
)

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
  return(results_df)
}

analyze_all_allelic_imbalance_models <- function(df_merged) {
  # Make a copy of the input dataframe
  results_df <- df_merged
  # List of all DNA models with their model_type identifiers
  dna_models <- list(
    list(model = DNA_bb_model1, type = "DNA_bb_model1"),
    list(model = DNA_null_model, type = "DNA_null_model")
  )
  # List of all RNA models with their model_type identifiers
  rna_models <- list(
    list(model = RNA_bb_model1, type = "RNA_bb_model1"),
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
  return(results_df)
}
# Process all models and get comprehensive results
TCGA_AI_table_final <- analyze_all_allelic_imbalance_models(df_merged = df)

dim(TCGA_AI_table_final)

p_value <- 0.05
TCGA_AI_table_final <- TCGA_AI_table_final %>%
  dplyr::mutate(
    pos_mRNA_AI = ( RNA_p_upper_null_model < p_value ),
    neg_mRNA_AI = ( RNA_p_lower_null_model < p_value ),
    no_mRNA_AI = !(pos_mRNA_AI | neg_mRNA_AI),
    pos_CNA_AI = (DNA_p_upper_null_model < p_value &
            DNA_p_upper_bb_model1 >= p_value),
    neg_CNA_AI = (DNA_p_lower_null_model < p_value &
            DNA_p_lower_bb_model1 >= p_value),
    no_CNA_AI = !(pos_CNA_AI | neg_CNA_AI),
    pos_reg_AI = ( RNA_p_upper_bb_model1 < p_value &
            DNA_p_upper_null_model >= p_value),
    neg_reg_AI = ( RNA_p_lower_bb_model1 < p_value &
            DNA_p_lower_null_model >= p_value),
    no_reg_AI = !(pos_reg_AI | neg_reg_AI),
    all_mRNA_AI = (pos_mRNA_AI | neg_mRNA_AI),
    all_CNA_AI = (pos_CNA_AI | neg_CNA_AI),
    all_reg_AI = (pos_reg_AI | neg_reg_AI)
  ) %>%
  filter(!is.na(pos_reg_AI)) %>%
  # define the three reg-AI classes
  mutate(
    reg_AI_type = case_when(
      pos_reg_AI ~ "Pos reg-AI",
      neg_reg_AI ~ "Neg reg-AI",
      TRUE       ~ "No reg-AI"
    ),
    # make sure it's a factor in the order you want
    reg_AI_type = factor(
      reg_AI_type,
      levels = c( "Pos reg-AI", "No reg-AI", "Neg reg-AI"),
      labels = c("Pos reg-AI", "No reg-AI", "Neg reg-AI")
    )
  ) %>%
  # define the three CNA-AI classes
  mutate(
    CNA_AI_type = case_when(
      pos_CNA_AI ~ "Pos CNA-AI",
      neg_CNA_AI ~ "Neg CNA-AI",
      TRUE       ~ "No CNA-AI"
    ),
    # make sure it's a factor in the order you want
    CNA_AI_type = factor(
      CNA_AI_type,
      levels = c( "Pos CNA-AI", "No CNA-AI", "Neg CNA-AI"),
      labels = c("Pos CNA-AI", "No CNA-AI", "Neg CNA-AI")
    )
  ) %>%
  # define the three mRNA-AI classes
  mutate(
    mRNA_AI_type = case_when(
      pos_mRNA_AI ~ "Pos mRNA-AI",
      neg_mRNA_AI ~ "Neg mRNA-AI",
      TRUE       ~ "No mRNA-AI"
    ),
    # make sure it's a factor in the order you want
    mRNA_AI_type = factor(
      mRNA_AI_type,
      levels = c( "Pos mRNA-AI", "No mRNA-AI", "Neg mRNA-AI"),
      labels = c("Pos mRNA-AI", "No mRNA-AI", "Neg mRNA-AI")
    )
  )

# B) Calculate % reg-AI vs CNA-AI
# Filterings:
# Our analysis --> without correcting by purity + only BRCA + PIK3CA
# we skip --> nonsynonymous

library(dplyr)
library(stats)   # for qnorm()

pancan_summary_df <- TCGA_AI_table_final %>%
  # filter(cancer == "TCGA-BRCA")                         # already BRCA
  filter(Gene.Symbol == "PIK3CA") %>%
  filter(SNV_varity_new %in% c("High Impact Missense",
                               "Low Impact Missense")) %>%
  filter((DNA_alt + DNA_ref) >= 30,
         (RNA_alt + RNA_ref) >= 30,
         all_mRNA_AI) %>%
  mutate(
    AI_source = case_when(
      all_reg_AI & !all_CNA_AI ~ "reg-AI",
      !all_reg_AI & all_CNA_AI ~ "CNA-AI"
    )
  ) %>%
  filter(!is.na(AI_source)) %>%
  group_by(AI_source) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    total_exclusive_snvs_pancan = sum(n),
    prop   = n / total_exclusive_snvs_pancan,            # proportion (0–1)
    se     = sqrt(prop * (1 - prop) /
                  total_exclusive_snvs_pancan),          # Wald standard error
    z      = qnorm(0.975),                               # 1.96 for 95 % CI
    ci_low = pmax(0, prop - z * se),                     # lower bound, clipped
    ci_high= pmin(1, prop + z * se),                     # upper bound, clipped
    percent      = prop      * 100,                      # back to %
    percent_low  = ci_low   * 100,
    percent_high = ci_high  * 100
  ) %>%
  select(AI_source, n, total_exclusive_snvs_pancan,
         percent, percent_low, percent_high)

pancan_summary_df

# # Save
# total_reg_AI <- paste0(round(pancan_summary_df[2,"percent"],1),"%")
# total_CNA_AI <- paste0(round(pancan_summary_df[1,"percent"],1),"%")