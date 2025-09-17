library(dplyr)
library(ggplot2)
source("growth_curves/make_site_metadata.r")
source("R_scripts/ggplot_theme.txt")

#######################
#Read data

#makes site_metadata
site_metadata
nrow(site_metadata)
site_metadata$Site %>% unique %>% sort

gr_stat.df = read.table(file = "data/summary_tables/growth_rate_stats.txt", header = T, sep = "\t")
nrow(gr_stat.df)
gr_stat.df$Site = sapply(strsplit(gr_stat.df$iso_name, " "), "[[", 1)
gr_stat.df$Site %>% unique %>% sort
gr_stat.df$Site %>% unique %>% length
#join the gr stats to metadata
gr_stat.df.metadata = left_join(
    gr_stat.df,
    site_metadata,
    by = "Site"
)

#check rows on joined and unique sites
nrow(gr_stat.df.metadata)
colnames(gr_stat.df.metadata)
unique(gr_stat.df.metadata$Site)
unique(gr_stat.df.metadata) %>% nrow
#all is well

#check for all na rows
apply(gr_stat.df.metadata, 1, FUN=function(x){ all(is.na(x)) } ) %>% sum

# join to sample level metadata. Just the species in this case
sample_ID_map = read.table("data/sample_metadata/sample_ID_mapping_all_samples_05052022.txt", header = T, sep = "\t")
sample_ID_map
sample_ID_map %>% nrow
colnames(sample_ID_map)[c(3, 5)] = c("iso_name", "spp")
sample_ID_map %>% select(iso_name, spp) %>% unique

gr_stat.df.metadata.sp = left_join(
    gr_stat.df.metadata,
    sample_ID_map %>% select(iso_name, spp) %>% unique,
    by = "iso_name"
)
nrow(gr_stat.df.metadata.sp)
is.na(gr_stat.df.metadata.sp$spp) %>% sum

gr_stat.df.metadata.sp %>% filter(is.na(spp)) %>% select(iso_name) %>% unique
######################
#

######################
#Filter to 30 and 35 C
#and calc means

gr.means.30_35 = gr_stat.df.metadata.sp %>%
    na.omit() %>%
    filter(temp == 30 | temp == 35) %>%
    group_by(spp, temp, iso_name) %>%
    summarize(gr.mean = mean(gr.est))

p1 = ggplot(gr.means.30_35 %>% filter(spp == "Nf"), aes(x = as.factor(temp), y = gr.mean)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.1)) +
    labs(x = expression(paste("Temperature ("~degree~C,")")), 
         y = expression(paste("Grwoth rate (mm"^-1, " day"^-1,")"))
    ) +
    my_gg_theme
p1

pdf("figures/Nf_30_35C_gr.pdf")
p1
dev.off()

ggplot(gr.means.30_35 %>% filter(spp == "Nd"), aes(x = as.factor(temp), y = gr.mean)) +
    geom_boxplot() +
    geom_point(position = position_jitter(width = 0.1)) +
    my_gg_theme

#distance cluster of 30C temps

#Nf
Nf.30 = gr.means.30_35 %>%
    filter(spp == "Nf" & temp == "30")
Nf.30


Nf.30.d = dist(Nf.30$gr.mean)
Nf.30.clust = hclust(Nf.30.d, method = "ward.D")
plot(Nf.30.clust)

pdf("figures/Nf.30C.clust.pdf")
plot(Nf.30.clust)
dev.off()

Nf.35 = gr.means.30_35 %>%
    filter(spp == "Nf" & temp == "35")
Nf.35


Nf.35.d = dist(Nf.35$gr.mean)
Nf.35.clust = hclust(Nf.35.d, method = "ward.D")
plot(Nf.35.clust)

pdf("figures/Nf.35C.clust.pdf")
plot(Nf.35.clust)
dev.off()


#Nd
Nd.30 = gr.means.30_35 %>%
    filter(spp == "Nd" & temp == "30")
Nd.30


Nd.30.d = dist(Nd.30$gr.mean)
Nd.30.clust = hclust(Nd.30.d)
plot(Nd.30.clust)
