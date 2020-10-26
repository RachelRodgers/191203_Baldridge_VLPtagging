# 191203_Baldridge_VLPtagging_phage_analysis.R

library("tidyverse")

# For each sample, read in the CAT summary file
fileList <- list.files("../data/summary_files/")
sampleNames <- map_chr(.x = fileList, .f = ~ str_remove(.x, "_CAT.summary.txt"))
filePaths <- paste0("../data/summary_files/", fileList)
names(filePaths) <- sampleNames

catSummaries <- map(.x = filePaths,
                    .f = read.delim,
                    sep = "\t", stringsAsFactors = FALSE,
                    blank.lines.skip = TRUE, comment.char = "#",
                    header = FALSE,
                    col.names = c("rank", "clade", "num_contigs",
                                  "num_ORFs", "num_pos"))

superkingdomList <- vector(mode = "list", length = length(catSummaries))

for (i in 1:length(catSummaries)) {
  currentSummary <- catSummaries[[i]]
  currentName <- names(catSummaries[i])
  pivotSummary <- currentSummary %>% 
    select(rank, clade, num_contigs) %>% 
    filter(rank == "superkingdom") %>% 
    pivot_wider(names_from = rank, values_from = num_contigs) %>% 
    mutate("sample" = currentName) %>% 
    rename(num_contigs = superkingdom)
  superkingdomList[[i]] <- pivotSummary
  names(superkingdomList)[i] <- currentName
}

# bind the DFs together
superkingdomsDF <- Reduce(f = function(df1, df2) {rbind(x = df1, y = df2,
                                                        make.row.names = FALSE)},
                          x = superkingdomList)

# add total numbers of contigs so we can make percentages for different categories
superkingdomsSum <- superkingdomsDF %>% 
  group_by(sample) %>% 
  mutate("total_contigs" = sum(num_contigs),
         "percent" = num_contigs/total_contigs)

superkingdomsSum$clade <- factor(superkingdomsSum$clade,
                                    levels = c("Bacteria", "Eukaryota",
                                               "Archaea",
                                               "not classified", 
                                               "Viruses"))


# Make named vector of samples, named by their total contigs
sampleLabelsVec <- superkingdomsSum %>% 
  select(sample, total_contigs) %>% 
  unique() %>% 
  arrange(total_contigs) %>% 
  deframe()

sampleOrder <- names(sampleLabelsVec)

superkingdomsSum$sample <- factor(superkingdomsSum$sample,
                                  levels = sampleOrder)

cbPaletteGrey <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                   "#0072B2", "#D55E00", "#CC79A7")

ggplot(superkingdomsSum, aes(x = sample, y = percent, fill = clade)) +
  geom_bar(stat = "identity", color = "black") +
  ylab("Percent of Contigs > 1000bp") +
  xlab("Sample (labeled by total number of contigs)") +
  scale_fill_manual(values = cbPaletteGrey[c(3,4,5,1,7)]) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(labels = sampleLabelsVec)


samplesWithViruses <- superkingdomsSum %>% 
  filter(clade == "Viruses") %>% 
  select(sample)

View(catSummaries$Baldridge_CF912K06_9F_SIC_index_0537_SIC_index_0551_AGCAGATAC_CCACGGTCT_S122_L001)

