library(dplyr)

###################
#data

#Read the GR data
gr_stat.df = read.table(file = "data/summary_tables/growth_rate_stats.txt", header = T, sep = "\t")
isonames = gr_stat.df$iso_name %>% unique() 

# join to sample level metadata. Just the species in this case
# sample metadata
sample_ID_map = read.table("data/sample_metadata/sample_ID_mapping_all_samples_05052022.txt", header = T, sep = "\t")
head(sample_ID_map)

write.table(sample_ID_map %>% filter(Spp_based_on_ITS_map == "Nf" & sample %in% isonames) %>% select(Sequence_label), "data/sample_metadata/sequence_IDs.Nf.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(sample_ID_map %>% filter(Spp_based_on_ITS_map == "Nd" & sample %in% isonames) %>% select(Sequence_label), "data/sample_metadata/sequence_IDs.Nd.txt", row.names = F, col.names = F, quote = F, sep = "\t")
