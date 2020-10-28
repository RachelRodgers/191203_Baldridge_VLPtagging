# epg5_hecatomb_euk_stdCountTables.R

# Script will produce:
#   1.) A standardized count table at the species level
#   2.) A standardized count table at the genus level

library("phyloseq")
library("data.table")
library("tidyverse")
library("speedyseq")

load("../data/RData_Objects/191203_Baldridge_VLPtagging_euk_physeqRaw.R")

# Baltimore virus classifications
baltimore <- read.delim("../data/2020_07_27_Viral_classification_table_ICTV2019.txt")

#----- Function Definitions -----#

StandardizeCounts <- function(physeqObject, viralClassTable, statColNames,
                              librarySizeTable) {
  physeqObject %>% 
    speedyseq::psmelt() %>% 
    tibble::as_tibble() %>% 
    mutate(across(all_of(statColNames), as.character),
           across(all_of(statColNames), as.numeric),
           across(where(is.character), as.factor)) %>% 
    merge(viralClassTable, by = "Family", all = TRUE) %>% 
    filter(!is.na(OTU)) %>% 
    merge(librarySizeTable, by = "Sample", all = TRUE) %>% 
    mutate(proportional_abundance = Abundance/library_size,
           scaled_abundance_min = (min(librarySizeTable$library_size)*proportional_abundance),
           scaled_abunance_mean = (mean(librarySizeTable$library_size)*proportional_abundance),
           scaled_abundance_median = (median(librarySizeTable$library_size)*proportional_abundance))
  
}

#----- Standardize Counts in cleanCountTable to Per-Sample Library Size -----#

librarySize <- enframe(colSums(cleanCountTable %>% select(-"id")),
                       name = "Sample", value = "library_size") %>% 
  as.data.frame()

#----- Generate Standardized Counts Objects: Species & Genus -----#

# Melt physeqRaw to generate a table where every row is a unique sample-taxon
#   pair, and columns contain sample metadata and alignment data.
#   Once melted, convert the numeric alignment stat columns back (stored
#   in the numStatVars variable) to numeric, and convert character columns to 
#   factors. Add Baltimore virus classification information, and library
#   size information. Use library size to generate scaled abundance information.

physeqMeltSpecies <- StandardizeCounts(physeqObject = physeqRaw,
                                       viralClassTable = baltimore,
                                       statColNames = numStatVars,
                                       librarySizeTable = librarySize) # 13,390 taxa x 27 vars

physeqRawGenus <- speedyseq::tax_glom(physeqRaw, "Genus", NArm = FALSE) # 98 x 65
physeqMeltGenus <- StandardizeCounts(physeqObject = physeqRawGenus,
                                     viralClassTable = baltimore,
                                     statColNames = numStatVars,
                                     librarySizeTable = librarySize) # 6,370 taxa x 27 vars


#----- Save Data  -----#

saveRDS(physeqMeltSpecies, "../data/RData_Objects/physeqMeltSpecies.RDS")
saveRDS(physeqMeltGenus, "../data/RData_Objects/physeMeltGenus.RDS")

save.image("../data/RData_Objects/191203_Baldridge_VLPtagging_euk_stdCountTables.R")

writeLines(capture.output(sessionInfo()),
           "191203_Baldridge_VLPtagging_euk_stdCountTables_session_info.txt")
Sys.Date()
getwd()
sessionInfo()

