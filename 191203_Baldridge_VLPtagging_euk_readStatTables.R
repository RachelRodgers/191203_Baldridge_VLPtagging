# epg5_hecatomb_euk_readStatsTable.R

# Script will produce:
#   1.) A data frame of individual reads and their alignment stats

library("phyloseq")
library("data.table")
library("tidyverse")
library("speedyseq")

load("../data/RData_Objects/191203_Baldridge_VLPtagging_euk_stdCountTables.R")

# UniProt target descriptions
descriptions <- read.delim("../data/uniprot.description.db", header = FALSE,
                           col.names = c("target", "protein"),
                           stringsAsFactors = FALSE)

# Merge eukaryotic virus taxonomy table with alignment data and count table
mergedData <- Reduce(f = function(df1, df2) 
  {base::merge(df1, df2, by = "id", all.x = TRUE)},
  x = list(eukTaxTable, allAln, cleanCountTable)) # 532 x 85

#---- Generate Raw Values Phyloseq Object -----# 
# A phyloseq object without merged counts or summarized alignment stats

rawPSCounts <- mergedData %>%
  select(sample_names(physeqRaw)) %>% 
  otu_table(taxa_are_rows = TRUE) # 532 sequences x 65 samples

rawPSTax <- mergedData %>% 
  select(all_of(c(ranks, allAlnCols))) %>% 
  mutate(across(where(is.numeric), as.character)) %>% 
  as.matrix() %>% 
  tax_table() # 532 sequences x 20 variables (ranks + alignment stats + id column)

physeqAlign <- phyloseq(rawPSCounts, rawPSTax) # 532 taxa x 65 samples

#----- Generate Table of Individual Read Alignment Stats -----#

readAlignmentData <- physeqAlign %>%
  speedyseq::psmelt() %>% 
  tibble::as_tibble() %>% 
  mutate(across(all_of(allAlnCols), as.character),
         across(all_of(numStatVars), as.numeric),
         across(where(is.character), as.factor)) %>% 
  merge(descriptions, by = "target", all = TRUE) %>% 
  filter(!is.na(OTU)) %>% 
  merge(librarySize, by = "Sample", all = TRUE) %>% 
  mutate(proportional_abundance = Abundance/library_size,
         scaled_abundance_min = (min(librarySize$library_size)*proportional_abundance),
         scaled_abunance_mean = (mean(librarySize$library_size)*proportional_abundance),
         scaled_abundance_median = (median(librarySize$library_size)*proportional_abundance)) %>% 
  select(-target) %>% 
  mutate(alignment_length_adjusted = base::ifelse(query_type == "aa",
                                                  yes = alignment_length*3,
                                                  no = alignment_length*1)) # 34,580 x 29

# Generate a filterd table to remove rows with 0 abundance
readAlignmentDataFiltered <- readAlignmentData %>% 
  filter(Abundance > 0) # 533 x 29

#----- Save Data  -----#

saveRDS(readAlignmentData, "../data/RData_Objects/readAlignmentData.RDS")
saveRDS(readAlignmentDataFiltered, "../data/RData_Objects/readAlignmentDataFiltered.RDS")

save.image("../data/RData_Objects/191203_Baldridge_VLPtagging_euk_readStatsTable.R")

writeLines(capture.output(sessionInfo()),
           "191203_Baldridge_VLPtagging_euk_readStatsTable_session_info.txt")
Sys.Date()
getwd()
sessionInfo()
