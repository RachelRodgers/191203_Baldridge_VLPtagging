# 191203_Baldridge_VLPtagging_euk_physeqRaw.R

# Script will produce:
#   1.) A raw count physeq object

library("phyloseq")
library("data.table")
library("tidyverse")

#----- Metadata (Phage & Eukaryotic) -----#

#----- Eukaryotic Taxonomy Table -----#

# Load eukaryotic virus taxonomy table and remove any Bacteria sequences.
#   Replace any spaces and underscores which appear within a column with
#   hyphens "-" for easier parsing when summing counts across taxa of the
#   same name.
eukTaxTable <- fread(file = "../data/viruses_tax_table.tsv", header = TRUE,
                     stringsAsFactors = TRUE) %>% 
  filter(Kingdom != "Bacteria") %>% 
  mutate(across(where(is.factor), ~ str_replace_all(., "\\s|_", "-"))) %>% 
  as.data.frame() # 532 eukaryotic viruses
eukTaxTable$id <- as.character(eukTaxTable$id)

#----- Full Count Table (Phage & Eukaryotic) -----#

# Contains all sequences - generated from all samples;
#   Remove sequence column and replace NAs with 0.
fullCountTable <- fread(file = "../data/seqtable.all", 
                        header = TRUE, sep = "\t") %>% 
  as.data.frame() # 270,672 total sequences

cleanCountTable <- fullCountTable %>% 
  select(-sequence) %>% 
  mutate(across(where(is.numeric),~ replace(., is.na(.), 0)))
cleanCountTable$id <- as.character(cleanCountTable$id)

#----- Alignment Statistics -----#

# aa checked alignment:
aaAln <- fread(file = "../data/aa.aln.m8", header = TRUE, sep = "\t") %>% # 27,022 x 13
  rename(id = query) %>% 
  mutate(query_type = "aa") %>% 
  as.data.frame()

# nt checked alignment:
ntAln <- fread(file = "../data/nt.aln.m8", header = TRUE, sep = "\t") %>%  # 40 x 12
  rename(id = query) %>% 
  mutate(query_type = "nt") %>% 
  as.data.frame()

# bind aa and alignment tables by row:
allAln <- rbind(aaAln, ntAln) # 27,062 sequences
allAln$id <- as.character(allAln$id)

# grab the alingment stat variables from allAln (remove the "id" variable),
#   as well as only the numeric stat variables from allAln.
allAlnCols <- colnames(allAln)
allStatVars <- allAlnCols[allAlnCols != "id"]
numStatVars <- allAlnCols[map_lgl(allAln, is.numeric)]

#----- Generate Raw Count Phyloseq Object -----#

ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Merge eukaryotic viral taxonomy table with alignment statistics;
#   summarize alignment statistics across taxa of the same name.
#   The "target" column will necessarily drop off during grouping.
#   Table will likely contain fewer sequences than the eukaryotic
#   taxonomy table because we are collapsing shared taxa. This will be the 
#   tax_table for physeq (ranks and alignment statistics retained).
taxAlnAvg <- merge(eukTaxTable, allAln, by = "id") %>% 
  select(-id) %>% 
  mutate(across(all_of(ranks), as.character)) %>% 
  unite(col = "lineage", all_of(ranks), sep = "_") %>% 
  group_by(lineage) %>% 
  summarize(across(all_of(numStatVars), mean)) %>%
  separate(col = "lineage", into = ranks, sep = "_") %>% 
  mutate(across(where(is.numeric), as.character)) %>% 
  as.matrix() # 206 x 18 (ranks + aln stats)

# Merge eukaryotic viral taxonomy table with count table;
#   sum counts across taxa of the same name.  Table will like contain fewer 
#   sequences than the eukaryotic taxonomy table because we are collapsing 
#   across shared taxa. It should have the same number of rows as taxAlnAvg.
#   This will be the otu_table (count table) for physeq (counts retained, 
#   ranks removed).
taxCountsSum <- merge(eukTaxTable, cleanCountTable, by = "id") %>% 
  select(-id) %>% 
  unite("lineage", all_of(ranks), sep = "_") %>% 
  group_by(lineage) %>% 
  summarize(across(where(is.numeric), sum)) %>% 
  separate("lineage", all_of(ranks), sep = "_") %>% 
  select(-all_of(ranks)) %>% 
  as.matrix() # 206 x 65 (samples w/counts summed by lineage)

# Generate physeq object
physeqRaw <- phyloseq(otu_table(taxCountsSum, taxa_are_rows = TRUE),
                      tax_table(taxAlnAvg)) # 206 taxa x 65 samples

#----- Save Data -----#

dir.create("../data/RData_Objects")

saveRDS(physeqRaw, "../data/RData_Objects/physeqRaw.RDS")
save.image("../data/RData_Objects/191203_Baldridge_VLPtagging_euk_physeqRaw.R")

writeLines(capture.output(sessionInfo()),
           "191203_Baldridge_VLPtagging_euk_physeqRaw_session_info.txt")
Sys.Date()
getwd()
sessionInfo()





