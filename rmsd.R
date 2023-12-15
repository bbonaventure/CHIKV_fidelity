library(VariantAnnotation)
library(Biostrings)
library(VariantAnnotation)
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(reshape2)
library(stringr)
library(rtracklayer)


#RMSD calculated from variant calling file

# Get a list of BAM files matching the pattern "deduplicated*.bam"
bam_files <- list.files(pattern = "deduplicated.*\\.bam$")
#set thresholds
thresholds <- c(0.01, 0.005, 0.001, 0.0001)
#initialize a data frame to store the average RMSD and SD for each condition
rmsd_summary <- data.frame(condition = character(), average_RMSD = numeric(), SD_RMSD = numeric(), SEM_RMSD = numeric(), Threshold = character(), stringsAsFactors = FALSE)
 
 # Loop over each condition in the condition_list
 for (bam_file in bam_files){
   RMSD_condition <-  results_filtered[which(results_filtered$bam_file == bam_file), ]
   
   for (threshold in thresholds) {
     # Filter the tsv dataframe based on the current threshold
     RMSD_df <- RMSD_condition[which(RMSD_condition$ALT_FREQ > threshold), ]

   # Transforming the data
   transformed_data <- RMSD_df %>%
     mutate(
       A = ifelse(REF == "A", REF_DP, ifelse(ALT == "A", ALT_DP, 0)),
       T = ifelse(REF == "T", REF_DP, ifelse(ALT == "T", ALT_DP, 0)),
       C = ifelse(REF == "C", REF_DP, ifelse(ALT == "C", ALT_DP, 0)),
       G = ifelse(REF == "G", REF_DP, ifelse(ALT == "G", ALT_DP, 0)),
       DEL = ifelse(grepl("^-+", ALT), ALT_DP, 0),
       INS = ifelse(grepl("^\\++", ALT), ALT_DP, 0)) %>%
     group_by(REGION, POS, Threshold, condition) %>%
     summarise(
       A = sum(A, na.rm = TRUE),
       T = sum(T, na.rm = TRUE),
       C = sum(C, na.rm = TRUE),
       G = sum(G, na.rm = TRUE),
       DEL = sum(DEL, na.rm = TRUE),
       INS = sum(INS, na.rm = TRUE),
       .groups = 'drop'
     ) %>%
     mutate(TOTAL_DP = A + T + C + G + DEL + INS) 
   
   # Calculate probabilities for each nucleotide, including INS and DEL
   transformed_data <- transformed_data %>%
     mutate(prob_A = A / TOTAL_DP,
            prob_T = T / TOTAL_DP,
            prob_C = C / TOTAL_DP,
            prob_G = G / TOTAL_DP,
            prob_INS = INS / TOTAL_DP,
            prob_DEL = DEL / TOTAL_DP,
            prob_GAP = prob_INS + prob_DEL  # Combine insertion and deletion probabilities
     )
    # Calculate the prob_REF which is the max probability among A, T, C, G, INS, DEL for each position
   transformed_data <- transformed_data %>%
     mutate(prob_REF = pmax(prob_A, prob_T, prob_C, prob_G, prob_INS, prob_DEL, na.rm = TRUE))
   
   # Identify the nucleotide with the max probability at each position and calculate the RMSD components
   transformed_data <- transformed_data %>%
     rowwise() %>% # to perform operations row by row
       mutate(
         # Calculate squared difference for the reference nucleotide
         sq_diff_ref = (1 - prob_REF)^2,
         # Calculate squared differences for all nucleotides
         sq_diff_A = ifelse(prob_A != prob_REF, prob_A^2, 0),
         sq_diff_T = ifelse(prob_T != prob_REF, prob_T^2, 0),
         sq_diff_C = ifelse(prob_C != prob_REF, prob_C^2, 0),
         sq_diff_G = ifelse(prob_G != prob_REF, prob_G^2, 0),
         sq_diff_GAP = ifelse( prob_GAP != prob_REF, prob_GAP^2, 0),
         # Sum the squared differences excluding the ref nucleotide difference
         sq_diff_sum = sq_diff_A + sq_diff_T + sq_diff_C + sq_diff_G + sq_diff_GAP,
         # Calculate RMSD for the position
         rmsd = sqrt((sq_diff_ref + sq_diff_sum) / 5) # divide by 5 for the number of nucleotide types including GAP
       ) %>%
     ungroup() 
  
    # Vector to store RMSD values for each position
   rmsd_values <- numeric(nrow(transformed_data))
   
   #add rmsd_values of non mutated position. non mutation = rmsd = 0
   #get number of non-mutated position
   non_mutated_position_number <- 11852 - length(transformed_data$POS)
   rmsd_other_position_genome <- rep(0,non_mutated_position_number)
   condition_rmsd <- c(transformed_data$rmsd,rmsd_other_position_genome)
   
   
 # Calculate the average RMSD and its standard deviation across all positions
 average_rmsd <- mean(condition_rmsd)
 sd_rmsd <- sd(condition_rmsd)
 sem_rmsd <- sd_rmsd / sqrt(length(condition_rmsd))

 sample <- unique(transformed_data$condition)

 
 # Store the results
rmsd_summary <- rbind(condition_rmsd_summary, data.frame(condition = sample, average_RMSD = average_rmsd, SD_RMSD = sd_rmsd, SEM_RMSD = sem_rmsd, Threshold =threshold ))
 }
 }

write.table(condition_rmsd_summary, file = "RMSD_table.txt", sep = "\t", row.names = FALSE)
 
