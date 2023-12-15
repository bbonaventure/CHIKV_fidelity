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


#Analysis using the result table from variant calling script

# list of all mutation types from all data
all_mutation_types <- unique(results_filtered %>%
                               mutate(mutation_type = paste(REF, "to", ALT)) %>%
                               pull(mutation_type))

# Create an empty dataframe to store the results 
results_df <- data.frame( )
#set different threshold values
thresholds <- c(0.01, 0.005, 0.001, 0.0001)

##get the number of total variants, the one in coding sequences, and the one inducing nonsynonymous mutation
for (bam_file in bam_files){
  table_condition <-  results_filtered[which(results_filtered$bam_file == bam_file), ]
  
  for (threshold in thresholds) {
    # Filter the tsv dataframe based on the current threshold
    filtered_df <- table_condition[which(table_condition$ALT_FREQ > threshold), ]
    
    # Calculate the number of variants for the current threshold
    num_variants <- as.numeric(sum(nrow(filtered_df)))
    
    # Calculate the number of coding variants for the current threshold
    num_variants_coding_sequences <- nrow(subset(filtered_df, filtered_df$GFF_FEATURE != "" ))
    
    # Calculate the number of coding variants for the current threshold wihich change the Amino acid
    num_variants_nonsynonymous_mutation <- nrow(subset(filtered_df, filtered_df$GFF_FEATURE != "" & filtered_df$REF_AA != filtered_df$ALT_AA ))
    
    # get mutation type
    filtered_df <- filtered_df %>%
      mutate(mutation_type = paste(REF, "to", ALT))
    
    # Count occurrences of each mutation type
    mutation_counts <- filtered_df %>%
      group_by(mutation_type) %>%
      summarize(count = n(), .groups = 'drop')
    
    # Ensure all mutation types are represented
    all_mutation_df <- data.frame(mutation_type = all_mutation_types)
    mutation_counts <- full_join(all_mutation_df, mutation_counts, by = "mutation_type")
    
    # Replace NA with 0 if needed
    mutation_counts[is.na(mutation_counts)] <- 0
    
    # Reshape mutation count data to wide format
    wide_mutation_counts <- pivot_wider(mutation_counts, names_from = mutation_type, values_from = count, values_fill = list(count = 0))
    
    mutation_df <- data.frame(Threshold = threshold, 
                               Condition = bam_file, 
                               NumVariants = num_variants, 
                               NumVariantsCodingSequence = num_variants_coding_sequences, 
                               NumVariantsInducingMutation = num_variants_nonsynonymous_mutation)
  
  # Merge with mutation counts
  merged_df <- merge(mutation_df, wide_mutation_counts, by = NULL)
  results_df <- rbind(results_df, merged_df)
  }
}

write.table(results_df, file = "Number_of_variants.txt", sep = "\t", row.names = FALSE)
